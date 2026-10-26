#pragma once
// Trivial floor division.
//
// Handles two categories:
// 1. Divisors of form +-2^k (truly trivial).
//
//    1:  x
//   -1: -x
//    2:  x >> 1
//   -2: -(x >> 1) - (x & 1)     (== floor(x / -2))
//   ... and so on for every power of two.
//
// 2. Also handled: |d| > 2^(W-2) (about half of all divisors). Then the quotient
// lies in [-2, 1] and is found by comparing against a few thresholds, e.g.
// for int8_t, d = 90:   (x >= 90) - (x < 0) - (x < -90)
// and for d = -90:      (x <= -90) - (x > 0) - (x > 90)
//
// Case 2. above has winning latency, but is not as high throughput as some
// alternatives.
//
// Any other divisor throws (a compile error when evaluated at compile time,
// e.g. floordiv_ct<D>).
#include <stdexcept>

#include "algorithms/base.hpp"

namespace fastdiv {

struct Trivial {
    static constexpr const char* name = "trivial";

    enum class Kind : uint8_t { One, NegOne, Shift, NegShift, Threshold, NegThreshold };

    template <class T> struct Params {
        T divisor;
        Kind kind;
        int shift;
        T m;  // Threshold kinds: |d|
    };

    // True for d == +-2^k (including +-1). The other algorithms fall back to
    // Trivial for exactly these divisors.
    template <class T> static constexpr bool is_pow2(T d) {
        const uwide_t<T> dabs = abs_wide<T>(d);
        return d != 0 && (dabs & (dabs - 1)) == 0;
    }

    template <class T> static constexpr Params<T> params(T d) {
        Params<T> p{};
        p.divisor = d;
        const uwide_t<T> dabs = abs_wide<T>(d);
        const bool neg = d < 0;
        const bool pow2 = is_pow2(d);
        const bool big = dabs > (uwide_t<T>(1) << (traits<T>::bits - 2));
        if (!pow2 && !big) {
            throw std::invalid_argument(
                "trivial: divisor must be +-2^k or have |d| > 2^(W-2)");
        }
        if (!pow2) {
            p.kind = neg ? Kind::NegThreshold : Kind::Threshold;
            p.m = T(dabs);
        }
        else if (dabs == 1) {
            p.kind = neg ? Kind::NegOne : Kind::One;
        }
        else {
            p.kind = neg ? Kind::NegShift : Kind::Shift;
            p.shift = ctz_wide<T>(dabs);
        }
        return p;
    }

    template <class T>
    FASTDIV_INLINE static inline T floordiv(T x, Params<T> p) {
        using U = uns_t<T>;
        switch (p.kind) {
            case Kind::One:
                return x;
            case Kind::NegOne:
                return T(U(0) - U(x));
            case Kind::Shift:
                return T(x >> p.shift);
            case Kind::NegShift: {
                const U mask = U((U(1) << p.shift) - 1);
                return T(U(0) - U(x >> p.shift) - U((U(x) & mask) != 0));
            }
            case Kind::Threshold:
                return T((x >= p.m) - (x < 0) - (x < T(-p.m)));
            case Kind::NegThreshold:
                return T((x <= T(-p.m)) - (x > 0) - (x > p.m));
        }
        return 0;
    }

    template <class T> static std::string describe(Params<T> p) {
        switch (p.kind) {
            case Kind::One:      return "x";
            case Kind::NegOne:   return "-x";
            case Kind::Shift:    return "x >> " + std::to_string(p.shift);
            case Kind::NegShift: return "-(x >> " + std::to_string(p.shift) + ") - (x & mask != 0)";
            case Kind::Threshold:    return "thresholds  |D|=" + std::to_string((long long)p.m);
            case Kind::NegThreshold: return "thresholds  |D|=" + std::to_string((long long)p.m) + "  (negative)";
        }
        return {};
    }

    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        constexpr Params<T> p = params<T>(D);
        using U = uns_t<T>;
        if constexpr (p.kind == Kind::One) {
            return x;
        }
        else if constexpr (p.kind == Kind::NegOne) {
            return T(U(0) - U(x));
        }
        else if constexpr (p.kind == Kind::Shift) {
            return T(x >> p.shift);
        }
        else if constexpr (p.kind == Kind::NegShift) {
            constexpr U mask = U((U(1) << p.shift) - 1);
            return T(U(0) - U(x >> p.shift) - U((U(x) & mask) != 0));
        }
        else if constexpr (p.kind == Kind::Threshold) {
            return T((x >= p.m) - (x < 0) - (x < T(-p.m)));
        }
        else {  // NegThreshold
            return T((x <= T(-p.m)) - (x > 0) - (x > p.m));
        }
    }
};

}  // namespace fastdiv
