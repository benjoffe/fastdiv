# fastdiv

Fast C++ implementations of signed floor division by a constant. Reference material for the blog post:

- **fast-floordiv** — https://www.benjoffe.com/fast-floordiv

## Algorithms

All live in `algorithms/` and compute `floor(x / d)` for signed `x`, for any non-zero signed divisor `d`, at 8, 16, 32 and 64 bits.

| Name | Header | Notes |
|---|---|---|
| `traditional` | `traditional.hpp` | divmod, then fix up the quotient when the remainder is non-zero and signs differ. The baseline. |
| `gm-v1` | `granlund_montgomery_v1.hpp` | Granlund & Montgomery (1994): `(x + (x < 0)) / D - (x < 0)` for `D > 0`. |
| `gm-v2` | `granlund_montgomery_v2.hpp` | Granlund & Montgomery XOR trick for positive divisors; negative divisors fall back to `gm-v1`. |
| `robison-v1` | `robison.hpp` | Robison (ARITH17, 2005) multiply-add, `(A*x + B) >> S`, in its "add-back" form: native-width multiply, with the correction for a too-wide multiplier added back to the high half. |
| `robison-v2` | `robison.hpp` | Robison multiply-add. Double-width multiply at 8–32 bits; at 64 bits, a native multiply with the carry from adding `B` computed explicitly. |
| `benjoffe` | `benjoffe.hpp` | Search over round-up / round-down multipliers with a one-compare bias, or a sign-selected multiplier, whichever is cheapest. |

Each algorithm has a runtime form (a params struct computed from the divisor, used by `div_test`) and a compile-time form `floordiv_ct<D>(x)` (used by `div_bench`). Powers of two (and ±1) are handled by `algorithms/trivial.hpp`.

## Build

Requires: CMake ≥ 3.18, a C++26-capable compiler (`__int128` support recommended; portable fallbacks are in `util/`), internet access on first build (fetches Google Benchmark).

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Using a specific compiler:

```sh
cmake -S . -B build-gcc -DCMAKE_CXX_COMPILER=g++-15
cmake --build build-gcc -j
```

Windows (clang or MSVC, Ninja):

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS="-DNOMINMAX"
cmake -S . -B build-msvc -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=cl -DCMAKE_CXX_FLAGS="/DNOMINMAX"
```

## Test

Verifies each algorithm against the naive reference, for every divisor of the chosen bit-width:

```sh
./build/bin/div_test -bits=16 -algo=benjoffe
./build/bin/div_test -bits=64 -algo=robison-v2 -d=-7 -quick=1000000
```

| Flag | Meaning |
|---|---|
| `-bits=8\|16\|32\|64` | Bit-width (required) |
| `-algo=NAME\|all` | `traditional`, `gm-v1`, `gm-v2`, `robison-v1`, `robison-v2`, `benjoffe` (required) |
| `-d=N` | Test only divisor `N` |
| `-boundary` | Only numerators adjacent to multiples of the divisor |
| `-quick[=N]` | `N` numerators each around zero, up to `INT_MAX`, down to `INT_MIN`, and random (default 1e9) |

Without `-d`, every divisor is tested. 8- and 16-bit are exhaustive over all numerators. 32-bit is allowed but takes a very long time, so use `-quick` or `-boundary`. 64-bit refuses to run without one of them, and testing every 64-bit divisor is not feasible, so pass `-d=N` there.

## Benchmark

Compile-time divisors at 64 bits: 100, 360, 86400, 100000, 10^7, 86400000 and 10^8, then 3..31 and -3..-31 (powers of two skipped).

```sh
./build/bin/div_bench
```

| Flag | Measures |
|---|---|
| *(none)* | Throughput — independent calls, CPU free to pipeline |
| `-latency` | Latency — serial dependency chain, each result feeds the next input |
| `-batch=N` | `N` operations per loop iteration (1..8) |

Additional flags: `-count=N` (input count, default 65536), `-repeat=N` (benchmark repetitions). Any other flag is passed to Google Benchmark, e.g. `--benchmark_filter=/d7$`. The `base` benchmark is the same loop with no division, to read off loop overhead.

Medians only:

```sh
./build/bin/div_bench -batch=8 -repeat=10 --benchmark_report_aggregates_only=true 2>&1 \
    | grep -E '^(Benchmark|---|Load|Run on|CPU|-+$)|_median'
```

## Benjoffe divisor statistics

`div_report` (`benjoffe` only) shows what fraction of divisors at each width fall into each `benjoffe` kind (unit, power-of-two, sign-selected, biased round-up/round-down):

```sh
./build/bin/div_report
./build/bin/div_report -min=-1000000 -max=1000000
```

## Attribution

The Granlund & Montgomery and Robison implementations are independent reimplementations of their published algorithms — see the respective header files for sources.

## License

MIT — Copyright 2026 Ben Joffe
