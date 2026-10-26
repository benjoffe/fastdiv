#pragma once
// Shared pieces for every algorithm in this folder.
//
// Every algorithm exposes the same shape:
//
//   struct Algo {
//     static constexpr const char* name;
//
//     // Runtime form -- lets the tester sweep every divisor without
//     // compiling one function per divisor.
//     template <class T> using Params = ...;
//     template <class T> static constexpr Params<T> params(T d);
//     template <class T> static T floordiv(T x, Params<T> p);
//     template <class T> static std::string describe(Params<T> p);
//
//     // Template form -- the divisor is a compile-time constant, so the
//     // params fold away and the compiler emits the ideal code. This is
//     // what the benchmarks use.
//     template <auto D> static decltype(D) floordiv_ct(decltype(D) x);
//   };
//
// Divisors of +-2^k (isTrivial) fall back to trivial.hpp in every algorithm
// except trivial itself.
//
// floordiv_ct<D> is literally `floordiv(x, params(D))` with the params
// constexpr, so the code that is benchmarked is the code that was tested.

#include <cstdint>
#include <string>

#include "util/int256.hpp"
#if !defined(__SIZEOF_INT128__)
#include <cstdlib>
#include "util/portable_int128.hpp"
#include "util/portable_uint128.hpp"
#endif

// Optimisation barrier: makes the compiler forget what it knows about `v`.
// Clang-only; GCC does not need it.
#if defined(__clang__)
#define FASTDIV_OPAQUE(v) __asm__("" : "+r"(v))
#else
#define FASTDIV_OPAQUE(v) ((void)0)
#endif

// Force inlining. MSVC ignores [[gnu::always_inline]], and in div_bench's huge
// TU it then leaves the hot path as an out-of-line call per operation.
#if defined(_MSC_VER) && !defined(__clang__)
#define FASTDIV_INLINE [[msvc::forceinline]]
#else
#define FASTDIV_INLINE [[gnu::always_inline]]
#endif

// Unreachable-in-practice failure; not constexpr, so it fails constant
// evaluation if ever reached there. MSVC has no __builtin_trap.
#if defined(__GNUC__) || defined(__clang__)
#define FASTDIV_TRAP() __builtin_trap()
#else
#define FASTDIV_TRAP() std::abort()
#endif

namespace fastdiv {

// Per-width types. T is always the signed type (int8_t .. int64_t).
template <class T> struct traits;
template <> struct traits<int8_t>  { using U = uint8_t;  using Wd = int16_t; using UWd = uint16_t; static constexpr int bits = 8;  };
template <> struct traits<int16_t> { using U = uint16_t; using Wd = int32_t; using UWd = uint32_t; static constexpr int bits = 16; };
template <> struct traits<int32_t> { using U = uint32_t; using Wd = int64_t; using UWd = uint64_t; static constexpr int bits = 32; };
template <> struct traits<int64_t> { using U = uint64_t; using Wd = i128;    using UWd = u128;     static constexpr int bits = 64; };

template <class T> using wide_t  = typename traits<T>::Wd;   // signed, 2x width
template <class T> using uwide_t = typename traits<T>::UWd;  // unsigned, 2x width
template <class T> using uns_t   = typename traits<T>::U;    // unsigned, same width

template <class T> constexpr T t_min() { return T(uns_t<T>(1) << (traits<T>::bits - 1)); }
template <class T> constexpr T t_max() { return T(~t_min<T>()); }

// Full double-width product: exactly wide_t<T>(a) * b where the wide type is
// native. For 64-bit on MSVC, i128 is a class type whose operator* is a full
// 128x128 multiply, so use the 64x64->128 intrinsics from
// util/portable_*int128.hpp instead.
template <class T>
FASTDIV_INLINE inline wide_t<T> mul_wide(T a, T b) {
#if !defined(__SIZEOF_INT128__)
    if constexpr (sizeof(T) == 8) {
        const int128_t r = s128_mul64(a, b);
        return i128(r.lo, uint64_t(r.hi));
    } else
#endif
    return wide_t<T>(a) * b;
}

// (a * b) >> s in the wide type, truncated to T: exactly T(mul_wide(a, b) >> s)
// where the wide type is native. On MSVC the class-type shift by a constant
// >= 64 compiles to a register shift (mov + asr), so shift the high word.
template <class T>
FASTDIV_INLINE inline T mul_shr(T a, T b, int s) {
#if !defined(__SIZEOF_INT128__)
    if constexpr (sizeof(T) == 8) {
        const int128_t r = s128_mul64(a, b);
        if (s >= 64) return T(r.hi >> (s - 64));
        return T((r.lo >> s) | (uint64_t(r.hi) << (64 - s)));
    } else
#endif
    return T(wide_t<T>(a) * b >> s);
}

template <class T>
FASTDIV_INLINE inline uns_t<T> umul_shr(uns_t<T> a, uns_t<T> b, int s) {
    using U = uns_t<T>;
#if !defined(__SIZEOF_INT128__)
    if constexpr (sizeof(T) == 8) {
        const uint128_t r = u128_mul64(a, b);
        if (s >= 64) return U(r.hi >> (s - 64));
        return U((r.lo >> s) | (r.hi << (64 - s)));
    } else
#endif
    return U(uwide_t<T>(a) * b >> s);
}

// |d| as an unsigned wide value (safe for d == INT_MIN).
template <class T> constexpr uwide_t<T> abs_wide(T d) {
    return uwide_t<T>(d < 0 ? -wide_t<T>(d) : wide_t<T>(d));
}

template <class T> constexpr int ctz_wide(uwide_t<T> v) {
    int c = 0;
    while (!(v & 1)) {
        v >>= 1;
        ++c;
    }
    return c;
}

// Reference floor division with a runtime divisor: exact, no tricks. Used by
// the tester as the oracle. Undefined for (INT_MIN, -1), which callers skip.
template <class T> inline T ref_floordiv(T x, T d) {
    int64_t q = int64_t(x) / int64_t(d);
    int64_t r = int64_t(x) % int64_t(d);
    if (r != 0 && ((r < 0) != (d < 0))) {
        --q;
    }
    return T(q);
}

}  // namespace fastdiv
