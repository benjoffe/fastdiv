// div_report -- what fraction of divisors land in each BenJoffe::params() kind:
// unit (|d|=1), power-of-two, sign-selected, or biased round-up/round-down.
//
//   ./div_report
//   ./div_report -min=N -max=N
//
// -min/-max (either or both) restrict every width to [max(min, type_min),
// min(max, type_max)]; a width whose intersection with that range is empty is
// skipped. 64-bit sweeps that range exhaustively if it's <= 2^32 divisors,
// otherwise falls back to a random sample of 2^32 drawn from within it.
//
// Runs smallest to largest, each stage live-updating one progress line while
// it runs (same CLEAR_LINE redraw div_test uses), then leaving the final
// counts in place before moving to the next width.
//
// 8/16-bit: exhaustive, single-threaded (trivially small either way). 32-bit:
// exhaustive, split across hardware threads. 64-bit: see above.
#ifdef _WIN32
#include <io.h>
#define isatty _isatty
#define STDOUT_FILENO _fileno(stdout)
#else
#include <unistd.h>
#endif

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "algorithms/benjoffe.hpp"

using namespace fastdiv;
using Clock = std::chrono::steady_clock;

namespace {

const char* const CLEAR_LINE = "\r\033[K";

std::atomic<uint64_t> g_progress{0};

// Workers call this every FLUSH_EVERY iterations (and once more with the
// remainder) instead of every iteration, so the atomic isn't a bottleneck.
constexpr uint64_t FLUSH_EVERY = 1 << 16;

struct Counts {
    uint64_t unit = 0, pow2 = 0, sign_selected = 0, round_up = 0, round_down = 0;
    // Sub-counts of the sign-selected family: divisors where BOTH rounding
    // directions have a valid shift, and exactly one of them lands cleanly
    // at S=W. These are sign-selected candidates that the preference logic
    // in benjoffe.hpp demotes to round-up/round-down (so they're already
    // counted in round_up/round_down above -- the sub-rows just surface how
    // often the demotion fires).
    uint64_t ss_wbit_down_only = 0, ss_wbit_up_only = 0;
    uint64_t total() const { return unit + pow2 + sign_selected + round_up + round_down; }
    Counts& operator+=(const Counts& o) {
        unit += o.unit;
        pow2 += o.pow2;
        sign_selected += o.sign_selected;
        round_up += o.round_up;
        round_down += o.round_down;
        ss_wbit_down_only += o.ss_wbit_down_only;
        ss_wbit_up_only += o.ss_wbit_up_only;
        return *this;
    }
};

template <class T>
inline void classify(T d, Counts& c) {
    constexpr int W = traits<T>::bits;
    const auto p = BenJoffe::params<T>(d);
    if (p.mode == BenJoffe::ModeTrivial) {
        if (p.trivial.kind == Trivial::Kind::One || p.trivial.kind == Trivial::Kind::NegOne) c.unit++;
        else c.pow2++;
        return;
    }
    if (p.mode == BenJoffe::ModeSignSelected) c.sign_selected++;
    else if (p.mode == BenJoffe::ModeRoundUp) c.round_up++;
    else c.round_down++;
    if (p.upShift && p.downShift) {
        if (p.downShift == W && p.upShift != W) c.ss_wbit_down_only++;
        else if (p.upShift == W && p.downShift != W) c.ss_wbit_up_only++;
    }
}

// Fast, independent per-thread streams (not cryptographic, just decorrelated).
inline uint64_t splitmix64(uint64_t& state) {
    uint64_t z = (state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

unsigned hw_threads() {
    const unsigned n = std::thread::hardware_concurrency();
    return n ? n : 1;
}

template <class F>
Counts parallel_for(unsigned nt, F&& per_thread) {
    std::vector<Counts> parts(nt);
    std::vector<std::thread> threads;
    for (unsigned t = 0; t < nt; t++) threads.emplace_back([&parts, t, &per_thread] { parts[t] = per_thread(t); });
    for (auto& th : threads) th.join();
    Counts total;
    for (auto& p : parts) total += p;
    return total;
}

// Exhaustive sweep of [lo, hi] (inclusive), split across threads. Works for
// int32_t and int64_t alike -- offsets are unsigned so this never overflows
// even when lo/hi sit at the type's extremes.
template <class T>
Counts sweep_range_parallel(int64_t lo, int64_t hi) {
    const unsigned nt = hw_threads();
    const uint64_t total = uint64_t(hi) - uint64_t(lo) + 1;
    const uint64_t chunk = total / nt;
    return parallel_for(nt, [&](unsigned t) {
        const uint64_t start = chunk * t;
        const uint64_t end = (t == nt - 1) ? total : (start + chunk);
        Counts c;
        uint64_t since_flush = 0;
        for (uint64_t off = start; off < end; off++) {
            const int64_t v = int64_t(uint64_t(lo) + off);
            if (v != 0) classify<T>(T(v), c);
            if (++since_flush == FLUSH_EVERY) {
                g_progress.fetch_add(since_flush, std::memory_order_relaxed);
                since_flush = 0;
            }
        }
        if (since_flush) g_progress.fetch_add(since_flush, std::memory_order_relaxed);
        return c;
    });
}

// 8/16-bit: exhaustive, single-threaded (trivially fast regardless of range).
template <class T>
Counts sweep_small_range(int64_t lo, int64_t hi) {
    Counts c;
    for (int64_t v = lo; v <= hi; v++) {
        if (v != 0) classify<T>(T(v), c);
    }
    g_progress.store(uint64_t(hi) - uint64_t(lo) + 1, std::memory_order_relaxed);
    return c;
}

// 64-bit, range too large to sweep: random sample of `count` nonzero
// divisors drawn uniformly from [lo, hi]. The full-range case (the default,
// no -min/-max) skips the modulo and uses the raw bits directly -- no bias
// since the range already spans the whole type.
Counts sample_range64(int64_t lo, int64_t hi, uint64_t count, uint64_t seed) {
    const bool full = (lo == INT64_MIN && hi == INT64_MAX);
    const uint64_t span = full ? 0 : (uint64_t(hi) - uint64_t(lo) + 1);
    const uint64_t reject_at = full ? 0 : (UINT64_MAX - UINT64_MAX % span);  // avoid modulo bias
    const unsigned nt = hw_threads();
    const uint64_t per = count / nt;
    return parallel_for(nt, [&](unsigned t) {
        const uint64_t n = (t == nt - 1) ? (count - per * (nt - 1)) : per;
        uint64_t state = seed ^ (0x9E3779B97F4A7C15ull * (t + 1));
        Counts c;
        uint64_t since_flush = 0;
        for (uint64_t i = 0; i < n; i++) {
            int64_t d;
            do {
                uint64_t r;
                do { r = splitmix64(state); } while (!full && r >= reject_at);
                d = full ? int64_t(r) : int64_t(uint64_t(lo) + r % span);
            } while (d == 0);
            classify<int64_t>(d, c);
            if (++since_flush == FLUSH_EVERY) {
                g_progress.fetch_add(since_flush, std::memory_order_relaxed);
                since_flush = 0;
            }
        }
        if (since_flush) g_progress.fetch_add(since_flush, std::memory_order_relaxed);
        return c;
    });
}

void report(const char* label, int W, const Counts& c, double secs) {
    const uint64_t n = c.total();
    auto pct = [&](uint64_t v) { return 100.0 * double(v) / double(n); };
    // Pad `s` so that any UTF-8 continuation bytes don't push numeric columns
    // left (printf's %-Ns pads to byte count, not display columns).
    auto row = [&](const char* s, uint64_t v) {
        int extra = 0;
        for (const char* p = s; *p; p++)
            if ((unsigned char)*p >= 0x80 && (unsigned char)*p < 0xC0) extra++;
        std::printf("  %-*s %14llu   %6.3f%%\n", 42 + extra, s, (unsigned long long)v, pct(v));
    };
    std::printf("%s%s   (n = %llu, %.1fs)\n", CLEAR_LINE, label, (unsigned long long)n, secs);
    row("unit / power-of-two", c.unit + c.pow2);
    row("round-up only", c.round_up);
    row("round-down only", c.round_down);
    row("sign-selected", c.sign_selected);
    char sub[96];
    std::snprintf(sub, sizeof sub, "  ├ with %d-bit shift for round-down only", W);
    row(sub, c.ss_wbit_down_only);
    std::snprintf(sub, sizeof sub, "  └ with %d-bit shift for round-up only", W);
    row(sub, c.ss_wbit_up_only);
    std::printf("\n");
}

// Runs `fn` on a background thread, redrawing a single progress line from
// g_progress (which `fn`'s workers update) until it finishes.
template <class F>
void run_stage(const std::string& label, int W, uint64_t total, F&& fn) {
    const bool live = isatty(STDOUT_FILENO);
    g_progress.store(0, std::memory_order_relaxed);
    const auto t0 = Clock::now();
    std::future<Counts> fut = std::async(std::launch::async, std::forward<F>(fn));
    while (fut.wait_for(std::chrono::milliseconds(100)) != std::future_status::ready) {
        if (!live) continue;
        const uint64_t done = g_progress.load(std::memory_order_relaxed);
        const double pct = 100.0 * double(done) / double(total);
        const double secs = std::chrono::duration<double>(Clock::now() - t0).count();
        std::printf("%s%s: %llu / %llu  (%.1f%%)  %.1fs", CLEAR_LINE, label.c_str(), (unsigned long long)done,
                    (unsigned long long)total, pct, secs);
        std::fflush(stdout);
    }
    const Counts c = fut.get();
    report(label.c_str(), W, c, std::chrono::duration<double>(Clock::now() - t0).count());
}

struct Options {
    bool has_min = false, has_max = false;
    int64_t min_v = 0, max_v = 0;
};

Options parse_args(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; i++) {
        const std::string a = argv[i];
        if (a.rfind("-min=", 0) == 0) { o.has_min = true; o.min_v = std::stoll(a.substr(5)); }
        else if (a.rfind("-max=", 0) == 0) { o.has_max = true; o.max_v = std::stoll(a.substr(5)); }
        else {
            std::fprintf(stderr, "Unknown argument: %s\nUsage: div_report [-min=N] [-max=N]\n", a.c_str());
            std::exit(2);
        }
    }
    if (o.has_min && o.has_max && o.min_v > o.max_v) {
        std::fprintf(stderr, "-min (%lld) must be <= -max (%lld)\n", (long long)o.min_v, (long long)o.max_v);
        std::exit(2);
    }
    return o;
}

// 8/16-bit stage: intersect the requested range with the type's own range.
template <class T>
void run_small_width(const char* name, const Options& o) {
    const int64_t tlo = int64_t(t_min<T>()), thi = int64_t(t_max<T>());
    const int64_t lo = o.has_min ? std::max(o.min_v, tlo) : tlo;
    const int64_t hi = o.has_max ? std::min(o.max_v, thi) : thi;
    if (lo > hi) { std::printf("%s: empty after intersecting with -min/-max, skipped\n\n", name); return; }
    char label[96];
    std::snprintf(label, sizeof label, "%s  (sweep [%lld, %lld])", name, (long long)lo, (long long)hi);
    run_stage(label, traits<T>::bits, uint64_t(hi) - uint64_t(lo) + 1,
              [lo, hi] { return sweep_small_range<T>(lo, hi); });
}

void run_32bit(const Options& o) {
    const int64_t tlo = INT32_MIN, thi = INT32_MAX;
    const int64_t lo = o.has_min ? std::max(o.min_v, tlo) : tlo;
    const int64_t hi = o.has_max ? std::min(o.max_v, thi) : thi;
    if (lo > hi) { std::printf("32-bit: empty after intersecting with -min/-max, skipped\n\n"); return; }
    char label[96];
    std::snprintf(label, sizeof label, "32-bit  (sweep [%lld, %lld])", (long long)lo, (long long)hi);
    run_stage(label, 32, uint64_t(hi) - uint64_t(lo) + 1,
              [lo, hi] { return sweep_range_parallel<int32_t>(lo, hi); });
}

void run_64bit(const Options& o) {
    const int64_t lo = o.has_min ? o.min_v : INT64_MIN;
    const int64_t hi = o.has_max ? o.max_v : INT64_MAX;
    if (lo > hi) { std::printf("64-bit: empty range (-min > -max), skipped\n\n"); return; }

    const bool full = (lo == INT64_MIN && hi == INT64_MAX);
    const uint64_t span = full ? 0 : (uint64_t(hi) - uint64_t(lo) + 1);  // 0 here means "exactly 2^64"
    const bool exceeds_2_32 = full || span > (uint64_t(1) << 32);

    char label[128];
    if (!exceeds_2_32) {
        std::snprintf(label, sizeof label, "64-bit  (sweep [%lld, %lld])", (long long)lo, (long long)hi);
        run_stage(label, 64, span, [lo, hi] { return sweep_range_parallel<int64_t>(lo, hi); });
    } else {
        // Fresh seed per run so successive invocations draw different samples.
        // std::random_device returns 32 bits on most platforms; combine two
        // draws for a full 64-bit seed. XOR with the clock as a backup in case
        // random_device is deterministic (older MinGW, some embedded libstdc++).
        std::random_device rd;
        const uint64_t seed = (uint64_t(rd()) << 32) ^ uint64_t(rd())
                            ^ uint64_t(std::chrono::steady_clock::now().time_since_epoch().count());
        if (full) {
            std::snprintf(label, sizeof label, "64-bit  (random sample of 2^32)");
        } else {
            std::snprintf(label, sizeof label, "64-bit  (random sample of 2^32 from [%lld, %lld])",
                          (long long)lo, (long long)hi);
        }
        run_stage(label, 64, uint64_t(1) << 32,
                  [lo, hi, seed] { return sample_range64(lo, hi, uint64_t(1) << 32, seed); });
    }
}

}  // namespace

int main(int argc, char** argv) {
    const Options o = parse_args(argc, argv);
    std::printf("fastdiv / BenJoffe::params kind distribution  (%u hardware threads)\n\n", hw_threads());

    run_small_width<int8_t>("8-bit", o);
    run_small_width<int16_t>("16-bit", o);
    run_32bit(o);
    run_64bit(o);

    return 0;
}
