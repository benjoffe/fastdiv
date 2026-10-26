// div_test -- exhaustively verify a division algorithm over every divisor of a
// given bit-width, using the runtime (params-struct) form of each algorithm.
// Aborts immediately on the first mismatch.
//
//   ./div_test -bits=16 -algo=benjoffe
//
//   -bits=8|16|32|64   required
//   -algo=NAME|all     required: traditional gm-v1 gm-v2 robison-v1 robison-v2 benjoffe
//   -d=N               test only the divisor N
//   -boundary          only numerators adjacent to multiples of the divisor
//   -quick[=N]         N numerators each: around zero, up to INT_MAX,
//                      down to INT_MIN, and random (default N = 1e9)
//
// Output: one table row per divisor. At 8/16-bit that is all there is. At
// 32/64-bit, when stdout is a terminal, the row in progress is redrawn in
// place with the numerator currently being tested, then replaced by the final
// row. Column widths are derived from the bit-width and the longest PARAMS
// text, and shrunk to the terminal width, so no width has dead space.
#ifdef _WIN32
#define NOMINMAX
#include <io.h>
#include <windows.h>
#define isatty _isatty
#define STDOUT_FILENO _fileno(stdout)
#else
#include <sys/ioctl.h>
#include <unistd.h>
#endif

#include <algorithm>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "algorithms/benjoffe.hpp"
#include "algorithms/traditional.hpp"
#include "algorithms/granlund_montgomery_v1.hpp"
#include "algorithms/granlund_montgomery_v2.hpp"
#include "algorithms/robison.hpp"

using namespace fastdiv;
using Clock = std::chrono::steady_clock;

namespace {

const char* const RESET     = "\033[0m";
const char* const AMBER     = "\033[38;5;216m";
const char* const DIM_AMBER = "\033[38;5;240m";
const char* const YELLOW    = "\033[93m";
const char* const RED       = "\033[38;5;203m";
const char* const GREEN     = "\033[38;5;121m";
const char* const CLEAR_LINE = "\r\033[K";

struct Options {
    int bits = 0;
    std::string algo;
    bool has_d = false;
    long long d = 0;
    bool boundary = false;
    bool quick = false;
    uint64_t quick_n = 1000000000ull;
};

volatile std::sig_atomic_t g_stop = 0;
void on_sigint(int) { g_stop = 1; }

// ---- number / text formatting ---------------------------------------------

std::string fmt_u(u128 v) {
    if (v == 0) return "0";
    std::string s;
    while (v) { s.insert(s.begin(), char('0' + int(v % 10))); v /= 10; }
    return s;
}
std::string fmt_i(i128 v) { return v < 0 ? "-" + fmt_u(u128(-v)) : fmt_u(u128(v)); }
std::string fmt_signed(i128 v) { return v < 0 ? fmt_i(v) : "+" + fmt_i(v); }

std::string with_commas(u128 v) {
    std::string s = fmt_u(v);
    for (int i = int(s.size()) - 3; i > 0; i -= 3) s.insert(size_t(i), ",");
    return s;
}

// 5-7 characters at any magnitude: "0.12s", "99.99s", "12m05s", "123h09m".
std::string fmt_time(double secs) {
    char b[32];
    if (secs < 99.995)        std::snprintf(b, sizeof b, "%.2fs", secs);
    else if (secs < 6000)     std::snprintf(b, sizeof b, "%dm%02ds", int(secs) / 60, int(secs) % 60);
    else                      std::snprintf(b, sizeof b, "%dh%02dm", int(secs) / 3600, int(secs) / 60 % 60);
    return b;
}

// Visible character count: skips ANSI escapes and UTF-8 continuation bytes.
size_t visible_len(const std::string& s) {
    size_t n = 0;
    bool esc = false;
    for (unsigned char c : s) {
        if (esc)               { if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) esc = false; }
        else if (c == '\033')  esc = true;
        else if ((c & 0xC0) != 0x80) n++;
    }
    return n;
}

std::string pad_r(const std::string& s, size_t n) {
    const size_t v = visible_len(s);
    return v >= n ? s : s + std::string(n - v, ' ');
}
std::string pad_l(const std::string& s, size_t n) {
    const size_t v = visible_len(s);
    return v >= n ? s : std::string(n - v, ' ') + s;
}

// Plain ASCII text cut to n columns, ending in an ellipsis.
std::string truncate(const std::string& s, size_t n) {
    return s.size() <= n ? s : s.substr(0, n - 1) + "\xe2\x80\xa6";
}

size_t term_width() {
#ifdef _WIN32
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    if (h != INVALID_HANDLE_VALUE && GetConsoleScreenBufferInfo(h, &csbi)) {
        int w = csbi.srWindow.Right - csbi.srWindow.Left + 1;
        if (w > 0) return size_t(w);
    }
    return 0;
#else
    struct winsize ws;
    return (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0) ? ws.ws_col : 0;
#endif
}

// Last-resort clip of an already-rendered row to the terminal width (keeps
// ANSI codes intact). Only reached if the table cannot be shrunk enough.
std::string clip(const std::string& s, size_t tw) {
    if (tw == 0 || visible_len(s) <= tw) return s;
    std::string out;
    size_t vis = 0;
    bool esc = false;
    for (unsigned char c : s) {
        if (esc) { out += char(c); if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z')) esc = false; }
        else if (c == '\033') { out += char(c); esc = true; }
        else if ((c & 0xC0) == 0x80) out += char(c);
        else { if (vis >= tw) break; out += char(c); vis++; }
    }
    return out + RESET;
}

// ---- the table --------------------------------------------------------------

struct Table {
    size_t w_div = 0, w_num = 0, w_chk = 0, w_time = 0, w_par = 0;  // w_par == 0: no PARAMS column
    size_t tw = 0;                                                    // terminal width, 0 = unlimited

    size_t total() const {
        const size_t ncols = w_par ? 5 : 4;
        return w_div + w_num + w_chk + w_time + w_par + 3 * ncols + 1;
    }
    // Shrink PARAMS (the only column with unbounded text) to fit the terminal.
    void fit() {
        if (!tw || !w_par || total() <= tw) return;
        const size_t over = total() - tw;
        w_par = std::max<size_t>(12, w_par > over ? w_par - over : 0);
    }

    std::string sep() const {
        auto d = [](size_t n) { return std::string(n + 2, '-'); };
        std::string s = std::string(AMBER) + "+" + d(w_div) + "+" + d(w_num) + "+" + d(w_chk) + "+" + d(w_time) + "+";
        if (w_par) s += d(w_par) + "+";
        return clip(s + RESET, tw);
    }
    std::string header() const {
        std::string s = std::string(AMBER) + "| " + pad_l("DIVISOR", w_div) + " | " + pad_l("NUMERATORS", w_num) +
                        " | " + pad_l("CHECKS", w_chk) + " | " + pad_l("TIME", w_time);
        if (w_par) s += " | " + pad_r("PARAMS", w_par);
        return clip(s + " |" + RESET, tw);
    }
    // `done` rows are coral; the in-progress row is uncoloured. A green
    // numerator cell marks a divisor that was checked against every numerator.
    std::string row(bool done, const std::string& div, const std::string& num, bool num_green,
                    uint64_t checks, double secs, const std::string& params) const {
        const char* col = done ? AMBER : RESET;
        std::ostringstream o;
        o << col << "| " << pad_l(div, w_div) << " | ";
        if (num_green) o << GREEN << pad_l(num, w_num) << col;
        else           o << pad_l(num, w_num);
        o << " | " << pad_l(with_commas(checks), w_chk) << " | " << pad_l(fmt_time(secs), w_time);
        if (w_par) {
            o << " | ";
            if (params.empty()) o << DIM_AMBER << pad_r("-", w_par) << col;
            else                o << pad_r(truncate(params, w_par), w_par);
        }
        o << " |" << RESET;
        return clip(o.str(), tw);
    }
};

std::string full_range(int W) {
    switch (W) {
        case 8:  return "-128 \xe2\x86\x92 127";
        case 16: return "-32768 \xe2\x86\x92 32767";
        case 32: return "-2\xc2\xb3\xc2\xb9 \xe2\x86\x92 2\xc2\xb3\xc2\xb9-1";
        default: return "-2\xe2\x81\xb6\xc2\xb3 \xe2\x86\x92 2\xe2\x81\xb6\xc2\xb3-1";
    }
}

// Longest PARAMS text any divisor in this run can produce. 8/16-bit: exact
// (every divisor is described up front). 32/64-bit: probed over the edges,
// powers of two and a random sample -- any outlier is truncated on display.
template <class T, class Algo>
size_t params_width(const Options& o) {
    constexpr int W = traits<T>::bits;
    using U = uns_t<T>;
    size_t m = 0;
    auto probe = [&](T d) {
        if (d == 0) return;
        m = std::max(m, Algo::template describe<T>(Algo::template params<T>(d)).size());
    };
    auto probe_mag = [&](u128 mag) {
        if (mag == 0 || mag > (u128(1) << (W - 1))) return;
        if (mag <= u128(t_max<T>())) probe(T(U(mag)));
        probe(T(U(U(0) - U(mag))));
    };
    if (o.has_d) {
        probe(T(o.d));
    } else if (W <= 16) {
        for (i128 v = i128(t_min<T>()); v <= i128(t_max<T>()); ++v) probe(T(v));
    } else {
        const u128 top = u128(1) << (W - 1);
        for (u128 k = 1; k <= 2048; ++k) { probe_mag(k); probe_mag(top - k + 1); }
        for (int k = 1; k < W; ++k) for (int j = -1; j <= 1; ++j) probe_mag((u128(1) << k) + u128(i128(j)));
        uint64_t s = 0x9E3779B97F4A7C15ull;
        for (int i = 0; i < 4096; ++i) {
            s ^= s << 13; s ^= s >> 7; s ^= s << 17;
            probe_mag(u128(s) % top + 1);
        }
    }
    return m;
}

// Runs one algorithm at one width over every divisor.
template <class T, class Algo>
int run_algo(const Options& o) {
    constexpr int W = traits<T>::bits;
    using U = uns_t<T>;
    const i128 lo = i128(t_min<T>()), hi = i128(t_max<T>());

    // -quick only applies when it actually samples less than the whole range.
    const bool use_quick = o.quick && u128(o.quick_n) * 4 < (u128(1) << W);
    const bool exhaustive = !use_quick && !o.boundary;
    std::string num_desc;
    if (use_quick)      num_desc = "quick 4 \xc3\x97 " + with_commas(o.quick_n);
    else if (o.boundary) num_desc = "boundary points";
    else                 num_desc = full_range(W);

    const bool live = isatty(STDOUT_FILENO) && W >= 32;

    // ---- column widths ----
    const size_t w_mag = fmt_u(u128(1) << (W - 1)).size() + 1;    // "-128", "-2147483648", ...
    const u128 dmin = o.has_d ? u128(o.d < 0 ? -i128(o.d) : i128(o.d)) : 1;
    u128 max_checks;
    if (use_quick)         max_checks = u128(o.quick_n) * 4;
    else if (o.boundary)   max_checks = 4 * ((u128(1) << (W - 1)) / dmin + 2) + 8;
    else                   max_checks = u128(1) << W;

    Table t;
    t.tw = term_width();
    t.w_div  = std::max<size_t>(7, o.has_d ? fmt_signed(i128(o.d)).size() : w_mag);
    t.w_num  = std::max({size_t(10), visible_len(num_desc), live ? w_mag : size_t(0)});
    t.w_chk  = std::max<size_t>(6, with_commas(max_checks).size());
    t.w_time = W <= 16 ? 6 : 7;
    if (const size_t pw = params_width<T, Algo>(o)) t.w_par = std::max<size_t>(6, pw);
    t.fit();

    std::cout << "\n" << AMBER << "Algorithm: " << Algo::name << "   Width: " << W << "-bit" << RESET << "\n\n";
    std::cout << t.sep() << "\n" << t.header() << "\n" << t.sep() << "\n";

    const auto t_start = Clock::now();
    uint64_t total_checks = 0, divisors_done = 0;

    auto finish = [&](bool interrupted) {
        std::cout << CLEAR_LINE << t.sep() << "\n";
        const double secs = std::chrono::duration<double>(Clock::now() - t_start).count();
        std::cout << "\n"
                  << (interrupted ? "Interrupted. " : "") << "Divisors tested: " << with_commas(divisors_done)
                  << "   Checks: " << with_commas(total_checks) << "   Time: " << std::fixed
                  << std::setprecision(1) << secs << "s\n";
        if (!interrupted) std::cout << GREEN << "ALL PASSED" << RESET << "\n";
    };

    auto test_divisor = [&](T d) {
        const auto params = Algo::template params<T>(d);
        const std::string params_str = Algo::template describe<T>(params);
        const std::string div_str = fmt_signed(i128(d));
        const auto t_div = Clock::now();
        auto last_live = t_div;
        uint64_t checks = 0;
        auto secs = [&]() { return std::chrono::duration<double>(Clock::now() - t_div).count(); };

        auto tick = [&](T x) {
            if (g_stop) { finish(true); std::exit(130); }
            if (!live) return;
            const auto now = Clock::now();
            if (now - last_live < std::chrono::milliseconds(100)) return;
            last_live = now;
            std::cout << CLEAR_LINE << t.row(false, div_str, fmt_i(i128(x)), false, checks, secs(), params_str)
                      << std::flush;
        };

        auto check = [&](T x) {
            if (d == T(-1) && x == t_min<T>()) return;  // INT_MIN / -1 overflows; undefined
            const T expect = ref_floordiv<T>(x, d);
            const T got = Algo::template floordiv<T>(x, params);
            if (expect != got) {
                std::cout << CLEAR_LINE << t.row(true, div_str, num_desc, false, checks, secs(), params_str) << "\n"
                          << t.sep() << "\n"
                          << RED << "\nFAIL: divisor " << fmt_i(i128(d)) << ", numerator " << fmt_i(i128(x))
                          << "\n      expected " << fmt_i(i128(expect)) << ", got " << fmt_i(i128(got)) << RESET
                          << "\n";
                std::exit(1);
            }
            ++checks;
            if ((checks & 0xFFFFF) == 0) tick(x);
        };

        // ---- numerators ----
        const u128 dabs = u128(abs_wide<T>(d));
        if (use_quick) {
            const uint64_t n = o.quick_n;
            auto run_from = [&](T start) {
                for (uint64_t c = 0; c < n; ++c) check(T(U(U(start) + U(c))));
            };
            run_from(T(U(U(0) - U(n / 2))));           // around zero
            run_from(T(U(U(t_max<T>()) - U(n - 1))));  // down from INT_MAX
            run_from(t_min<T>());                      // up from INT_MIN
            uint64_t s = 0x9E3779B97F4A7C15ull ^ uint64_t(d);
            for (uint64_t c = 0; c < n; ++c) {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                check(T(U(s)));
            }
        } else if (o.boundary) {
            // Numerators either side of every multiple of the divisor, swept in
            // ascending order from INT_MIN to INT_MAX (so the live row moves
            // steadily through the range). Every emit is strictly above the
            // previous one, which also drops duplicates.
            i128 last = lo - 1;
            auto emit = [&](i128 v) {
                if (v < lo || v > hi || v <= last) return;
                last = v;
                check(T(v));
            };
            const u128 kmax = (u128(hi) + 1) / dabs;
            emit(lo);
            emit(lo + 1);
            for (u128 k = kmax; k >= 1; --k) {          // negative multiples, largest magnitude first
                const i128 m = i128(k * dabs);
                emit(-m - 1);
                emit(-m);
            }
            emit(-1);
            emit(0);
            emit(1);
            for (u128 k = 1; k <= kmax; ++k) {          // positive multiples
                const i128 m = i128(k * dabs);
                emit(m - 1);
                emit(m);
            }
            emit(hi - 1);
            emit(hi);
        } else {
            for (i128 v = lo; v <= hi; ++v) check(T(v));
        }

        ++divisors_done;
        total_checks += checks;
        std::cout << CLEAR_LINE << t.row(true, div_str, num_desc, exhaustive, checks, secs(), params_str) << "\n";
        if (g_stop) { finish(true); std::exit(130); }
    };

    if (o.has_d) {
        test_divisor(T(o.d));
    } else {
        const u128 max_mag = u128(1) << (W - 1);
        for (u128 mag = 1; mag <= max_mag; ++mag) {
            if (mag <= u128(hi)) test_divisor(T(U(mag)));   // +2^(W-1) is not representable
            test_divisor(T(U(U(0) - U(mag))));
        }
    }
    finish(false);
    return 0;
}

template <class T>
int run_width(const Options& o, const std::string& algo) {
    if (algo == "traditional") return run_algo<T, Traditional>(o);
    if (algo == "gm-v1")       return run_algo<T, GranlundMontgomeryV1>(o);
    if (algo == "gm-v2")       return run_algo<T, GranlundMontgomeryV2>(o);
    if (algo == "robison-v1")  return run_algo<T, RobisonV1>(o);
    if (algo == "robison-v2")  return run_algo<T, RobisonV2>(o);
    if (algo == "benjoffe")    return run_algo<T, BenJoffe>(o);
    std::cerr << "Unknown algorithm: " << algo << "\n";
    return 2;
}

int run(const Options& o, const std::string& algo) {
    switch (o.bits) {
        case 8:  return run_width<int8_t>(o, algo);
        case 16: return run_width<int16_t>(o, algo);
        case 32: return run_width<int32_t>(o, algo);
        case 64: return run_width<int64_t>(o, algo);
    }
    return 2;
}

const char* const ALGOS[] = {"traditional", "gm-v1", "gm-v2", "robison-v1", "robison-v2", "benjoffe"};

void print_usage() {
    std::cout << YELLOW
              << "\nRun again with -bits and -algo\n"
              << "Options:\n"
              << "  -bits:     8, 16, 32, 64\n"
              << "  -algo:     traditional, gm-v1, gm-v2, robison-v1, robison-v2, benjoffe, all\n"
              << "  -d=N:      test only the divisor N\n"
              << "  -boundary: only test numerators adjacent to multiples of the divisor\n"
              << "  -quick[=N]: N numerators each around zero, up to INT_MAX, down to INT_MIN, and random\n"
              << "\nExample:\n"
              << "  ./build/bin/div_test -bits=16 -algo=benjoffe\n"
              << "  ./build/bin/div_test -bits=64 -algo=robison-v2 -d=-7 -quick=1000000\n"
              << RESET << "\n";
}

}  // namespace

int main(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if      (a.rfind("-bits=", 0) == 0)  o.bits = std::atoi(a.c_str() + 6);
        else if (a.rfind("-algo=", 0) == 0)  o.algo = a.substr(6);
        else if (a.rfind("-d=", 0) == 0)     { o.has_d = true; o.d = std::stoll(a.substr(3)); }
        else if (a == "-boundary")           o.boundary = true;
        else if (a == "-quick")              o.quick = true;
        else if (a.rfind("-quick=", 0) == 0) { o.quick = true; o.quick_n = std::stoull(a.substr(7)); }
        else { std::cerr << "Unknown argument: " << a << "\n"; print_usage(); return 2; }
    }

    if (o.bits == 0 || o.algo.empty()) { print_usage(); return 0; }
    if (o.bits != 8 && o.bits != 16 && o.bits != 32 && o.bits != 64) {
        std::cerr << "Bit-width must be 8, 16, 32 or 64.\n";
        return 2;
    }

    if (o.has_d) {
        const long long lim = o.bits == 64 ? INT64_MAX : (1ll << (o.bits - 1)) - 1;
        if (o.d == 0 || o.d > lim || o.d < -lim - 1) {
            std::cerr << "-d=" << o.d << " is not a non-zero " << o.bits << "-bit signed divisor\n";
            return 2;
        }
    }
    if (o.bits == 64 && !o.quick && !o.boundary) {
        std::cerr << RED << "64-bit cannot be tested exhaustively; pass -quick[=N] or -boundary." << RESET << "\n";
        return 2;
    }
    if (o.bits == 32 && !o.quick && !o.boundary && !o.has_d)
        std::cout << YELLOW << "Note: a full 32-bit sweep checks every numerator for every divisor and will take a very long\n"
                  << "time. Use -boundary or -quick to cut it down; Ctrl+C stops at the next row." << RESET << "\n";

    std::signal(SIGINT, on_sigint);

    if (o.algo == "all") {
        for (const char* a : ALGOS) {
            int rc = run(o, a);
            if (rc) return rc;
        }
        return 0;
    }
    return run(o, o.algo);
}
