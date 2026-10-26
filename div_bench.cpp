// div_bench -- compile-time-divisor benchmarks for every algorithm (64-bit).
//
// Divisors: 100, 360, 86400, 1e5, 1e7, 86400000, 1e8, then 3..31 and -3..-31
// (powers of 2 skipped). All
// algorithms are registered back-to-back per divisor so results sit together.
//
//   ./div_bench                  throughput (independent operations)
//   ./div_bench -latency         dependent chain: each result feeds the next input
//   ./div_bench -batch=N         N operations per loop iteration (1..8). In
//                                throughput mode this widens the independent
//                                work per iteration; in latency mode the chain
//                                runs through all N.
//   ./div_bench -count=N         number of input values (default 65536)
//   ./div_bench -repeat=N        shorthand for --benchmark_repetitions=N
//   (any other flag is passed through to Google Benchmark, e.g.
//    --benchmark_filter=/d7$)
//
// Measurement guards:
//  * every input goes through a register-only asm barrier (OPAQUE_IN) and every
//    result into an asm sink (SINK), so nothing is hoisted, folded or dropped
//    and no spill code is added to the loop;
//  * the loop is not auto-vectorised (NOVECTOR), so each algorithm is timed as
//    scalar code with the same loop shape;
//  * the `base` benchmark is the identical loop with no division, to read off
//    loop overhead;
//  * CMake disables identical-code folding / function outlining for this
//    target so equal-looking functions are not merged.
#include <benchmark/benchmark.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "algorithms/benjoffe.hpp"
#include "algorithms/traditional.hpp"
#include "algorithms/granlund_montgomery_v1.hpp"
#include "algorithms/granlund_montgomery_v2.hpp"
#include "algorithms/robison.hpp"

using namespace fastdiv;

// Register-only optimisation barriers. Unlike
// benchmark::DoNotOptimize, which clang satisfies by spilling the value to the
// stack (a store, plus a reload for inputs, per operation), these cost zero
// instructions, so the loop contains only the arithmetic under test.
//   OPAQUE_IN(x): compiler must assume x changed (blocks hoisting/CSE).
//   SINK(x):      compiler must compute x (blocks dead-code elimination).
//   OPAQUE_READ(x): the value of x, read in a way the compiler can't see through.
// MSVC has no inline asm (and none at all on x64/ARM64), so fall back to
// volatile accesses: OPAQUE_READ is one load per op (x is never written, so it
// is stored once per outer iteration and the loads don't chain through
// memory), SINK is one store. Not free -- `base` measures it too.
#if defined(_MSC_VER) && !defined(__clang__)
#define OPAQUE_IN(x) ((void)0)
#define OPAQUE_READ(x) (*static_cast<const volatile decltype(x)*>(&(x)))
#define SINK(x) do { [[maybe_unused]] volatile auto sink_ = (x); } while (0)
// The per-op lambdas below must inline too, or each op becomes a call.
#define LAMBDA_INLINE [[msvc::forceinline]]
#else
#define OPAQUE_IN(x) __asm__ volatile("" : "+r"(x))
#define SINK(x) __asm__ volatile("" : : "r"(x))
#define OPAQUE_READ(x) (x)
#define LAMBDA_INLINE
#endif

#ifdef __clang__
#define NOVECTOR _Pragma("clang loop vectorize(disable)")
#elif defined(_MSC_VER)
// MSVC's loop(no_vector) rejects range-for; the volatile SINK above already
// stops it vectorising.
#define NOVECTOR
#else
#define NOVECTOR _Pragma("GCC novector")
#endif

namespace {

size_t COUNT = 1 << 16;
std::vector<int64_t> inputs;
bool g_latency = false;
int g_batch = 1;

template <typename F, size_t... I>
FASTDIV_INLINE inline void unroll_impl(F&& f, std::index_sequence<I...>) {
    (((void)I, f()), ...);
}
template <int N, typename F>
FASTDIV_INLINE inline void unroll(F&& f) {
    unroll_impl(static_cast<F&&>(f), std::make_index_sequence<N>{});
}

void init_inputs() {
    std::mt19937_64 rng(0xDEADBEEF);
    inputs.resize(COUNT);
    for (auto& v : inputs) v = int64_t(rng());
}

}  // namespace

// External linkage is load-bearing: when BM_run is in the anonymous namespace
// its only callers are the RegisterBenchmark calls from static init, so Clang's
// IPA infers `cold`, which silently implies `minsize`. Under `minsize` the
// AArch64 cost model prefers a single `sdiv` over a multi-instruction magic-
// multiply, turning off divide-by-constant strength reduction and tanking every
// benchmark's numbers (0.46 ns/op "traditional/d3" is a hardware sdiv, not a
// strength-reduced constant divide). Giving BM_run external (weak_odr) linkage
// hides the complete call set from IPA, so the inference doesn't fire.
template <class Op, int BATCH>
[[gnu::hot, gnu::used, gnu::noinline]] void BM_run(benchmark::State& state) {
    double total_ns = 0;
    for (auto _ : state) {
        const auto t0 = std::chrono::steady_clock::now();
        if (g_latency) {
            int64_t result = 0;
            NOVECTOR
            for (int64_t d : inputs) {
                unroll<BATCH>([&] LAMBDA_INLINE {
                    d ^= result;
                    result = Op::apply(d);
                    SINK(result);
                });
            }
        } else {
            NOVECTOR
            for (int64_t d : inputs) {
                unroll<BATCH>([&] LAMBDA_INLINE {
                    OPAQUE_IN(d);
                    SINK(Op::apply(OPAQUE_READ(d)));
                });
            }
        }
        total_ns += std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count();
    }
    state.counters["ns/d"] = total_ns / (double(state.iterations()) * double(COUNT) * double(BATCH));
}

namespace {

// `base` is pure loop overhead and amortises with -batch.
struct Base {
    static int64_t& apply(int64_t& d) { return d; }
#if defined(_MSC_VER) && !defined(__clang__)
    static int64_t apply(const volatile int64_t& d) { return d; }  // OPAQUE_READ
#endif
};

template <class Algo, int64_t D>
struct Op {
    FASTDIV_INLINE static inline int64_t apply(int64_t d) { return Algo::template floordiv_ct<D>(d); }
};

constexpr bool is_pow2(long long v) {
    if (v < 0) v = -v;
    return (v & (v - 1)) == 0;
}

template <int BATCH, class Algo, int64_t D>
void reg(const std::string& suffix) {
    benchmark::RegisterBenchmark((std::string(Algo::name) + suffix).c_str(), BM_run<Op<Algo, D>, BATCH>);
}

// One divisor, every algorithm back-to-back.
template <int BATCH, int64_t D>
void register_divisor() {
    if constexpr (D != 0 && !is_pow2(D)) {
        const std::string suffix = std::string("/i64/d") + (D < 0 ? "n" : "") + std::to_string(D < 0 ? -D : D);
        reg<BATCH, Traditional, D>(suffix);
        reg<BATCH, GranlundMontgomeryV1, D>(suffix);
        reg<BATCH, GranlundMontgomeryV2, D>(suffix);
        reg<BATCH, RobisonV1, D>(suffix);
        reg<BATCH, RobisonV2, D>(suffix);
        reg<BATCH, BenJoffe, D>(suffix);
    }
}

template <int BATCH, int... Is>
void register_all_impl(std::integer_sequence<int, Is...>) {
    benchmark::RegisterBenchmark("base", BM_run<Base, BATCH>);
    register_divisor<BATCH, 100>();
    register_divisor<BATCH, 360>();
    register_divisor<BATCH, 86400>();
    register_divisor<BATCH, 100000>();
    register_divisor<BATCH, 10000000>();
    register_divisor<BATCH, 86400000>();
    register_divisor<BATCH, 100000000>();

    (register_divisor<BATCH, int64_t(3 + Is)>(), ...);
    (register_divisor<BATCH, int64_t(-3 - Is)>(), ...);
}

template <int BATCH>
void register_all() {
    register_all_impl<BATCH>(std::make_integer_sequence<int, 29>{});
}

}  // namespace

int main(int argc, char** argv) {
    std::vector<char*> bm_argv{argv[0]};
    static char rep_buf[64];
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        // Matches "-name=V", "--name=V" or "-name V"; returns V or nullptr.
        auto value = [&](const std::string& name) -> const char* {
            for (const char* dashes : {"-", "--"}) {
                const std::string key = dashes + name;
                if (a.compare(0, key.size(), key) != 0) continue;
                if (a.size() > key.size() && a[key.size()] == '=') return argv[i] + key.size() + 1;
                if (a.size() == key.size() && i + 1 < argc) return argv[++i];
            }
            return nullptr;
        };
        if (a == "-latency" || a == "--latency") {
            g_latency = true;
        } else if (const char* v = value("batch")) {
            g_batch = std::atoi(v);
        } else if (const char* v = value("count")) {
            COUNT = std::stoul(v);
        } else if (const char* v = value("repeat")) {
            std::snprintf(rep_buf, sizeof(rep_buf), "--benchmark_repetitions=%s", v);
            bm_argv.push_back(rep_buf);
        } else {
            bm_argv.push_back(argv[i]);
        }
    }

    init_inputs();

    switch (g_batch) {
        case 1: register_all<1>(); break;
        case 2: register_all<2>(); break;
        case 3: register_all<3>(); break;
        case 4: register_all<4>(); break;
        case 5: register_all<5>(); break;
        case 6: register_all<6>(); break;
        case 7: register_all<7>(); break;
        case 8: register_all<8>(); break;
        default:
            std::fprintf(stderr, "unsupported -batch=%d (use 1..8)\n", g_batch);
            return 1;
    }

    int bm_argc = int(bm_argv.size());
    benchmark::Initialize(&bm_argc, bm_argv.data());
    if (benchmark::ReportUnrecognizedArguments(bm_argc, bm_argv.data())) return 1;
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
