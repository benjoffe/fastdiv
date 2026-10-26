#pragma once
// Granlund & Montgomery -- the positive-divisor XOR trick:
//
//   s = x >> (W-1)            // x < 0 ? -1 : 0
//   a = x ^ s                 // negatives reflect to ~x = -1 - x  (>= 0)
//   b = (a * M) >> S          // unsigned mul-shift, a in [0, 2^(W-1))
//   return b ^ s
//
// Because `a` only spans half the unsigned range, a plain mul-shift works for
// every positive divisor, cooperative or not. Only specified for D > 0;
// negative divisors fall back to v1.
#include "algorithms/base.hpp"
#include "algorithms/granlund_montgomery_v1.hpp"

namespace fastdiv {

struct GranlundMontgomeryV2 {
    static constexpr const char* name = "gm-v2";

    template <class T> struct Params {
        bool isFallback;  // negative or +-2^k: handled by v1
        GranlundMontgomeryV1::Params<T> fallback;
        int shift;        // S
        uns_t<T> mul;     // M, M < 2^W
    };

    template <class T> static constexpr Params<T> params(T d) {
        constexpr int W = traits<T>::bits;
        using UW = uwide_t<T>;
        Params<T> p{};

        const UW D = UW(d);

        if (d > 0 && !Trivial::is_pow2(d)) {
            // Only run for non-power-of-2 positive divisors

            // Smallest S such that M = ceil(2^S / D) fits W bits and the rounding
            // error e = M*D - 2^S satisfies (2^(W-1) - 1) * e < 2^S. The bound is
            // satisfiable from S = W-1 (a only spans W-1 bits), but starting the
            // search there instead of W lets some divisors (e.g. 9, 19 for W=64)
            // land on S < W, which forces a cross-register combine (shrd/EXTR) to
            // splice the low and high halves of the product instead of a plain
            // high-word extraction. Monotonicity (valid at S implies valid at
            // S+1) means starting at W costs nothing but that instruction.
            const UW N = (UW(1) << (W - 1)) - 1;
            const UW limit = UW(1) << W;
            for (int S = W; S < 2 * W; S++) {
                const UW two_s = UW(1) << S;
                const UW M = two_s / D + 1;
                const UW e = D - two_s % D;
                if (M < limit && N * e < two_s) {
                    p.shift = S;
                    p.mul = uns_t<T>(M);
                    return p;
                }
            }
        }

        // Negative, and trivial divisors fallback
        p.isFallback = true;
        p.fallback = GranlundMontgomeryV1::params<T>(d);
        return p;
    }

    template <class T>
    FASTDIV_INLINE static inline T floordiv(T x, Params<T> p) {
        constexpr int W = traits<T>::bits;
        using U = uns_t<T>;
        if (p.isFallback) {
            return GranlundMontgomeryV1::floordiv<T>(x, p.fallback);
        }
        const T s = T(x >> (W - 1));
        const U a = U(x ^ s);
        const U b = umul_shr<T>(a, p.mul, p.shift);
        return T(T(b) ^ s);
    }

    template <class T> static std::string describe(Params<T> p) {
        return p.isFallback ? GranlundMontgomeryV1::describe<T>(p.fallback) : std::string{};
    }

    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        static_assert(D != 0, "divisor must be non-zero");
        constexpr Params<T> p = params<T>(D);
        constexpr int W = traits<T>::bits;
        using U = uns_t<T>;
        if constexpr (p.isFallback) {
            // Negative, and trivial divisors
            return GranlundMontgomeryV1::floordiv_ct<D>(x);
        }
        else {
            const T s = T(x >> (W - 1));
            const U a = U(x ^ s);
            const U b = umul_shr<T>(a, p.mul, p.shift);
            return T(T(b) ^ s);
        }
    }
};

}  // namespace fastdiv
