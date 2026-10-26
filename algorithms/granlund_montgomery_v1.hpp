#pragma once
// Granlund & Montgomery, "Division by Invariant Integers using Multiplication"
// (1994), page 6 -- the flexible signed floor division:
//
//   D > 0:  (x + (x < 0)) / D - (x < 0)
//   D < 0:  (x - (x > 0)) / D - (x > 0)
//
// With a compile-time D, the compiler turns the single division into a
// mul-shift. There are no parameters to solve; the "params" are just D.
#include "algorithms/base.hpp"
#include "algorithms/trivial.hpp"

namespace fastdiv {

struct GranlundMontgomeryV1 {
    static constexpr const char* name = "gm-v1";

    template <class T> struct Params {
        T divisor;
        bool isTrivial;  // +-2^k: handled by Trivial
        Trivial::Params<T> trivial;
    };

    template <class T> static constexpr Params<T> params(T d) {
        Params<T> p{};
        p.divisor = d;
        p.isTrivial = Trivial::is_pow2(d);
        if (p.isTrivial) {
            p.trivial = Trivial::params<T>(d);
        }
        return p;
    }

    template <class T>
    FASTDIV_INLINE static inline T floordiv(T x, Params<T> p) {
        if (p.isTrivial) {
            return Trivial::floordiv<T>(x, p.trivial);
        }
        const T D = p.divisor;
        if (D > 0) {
            const bool adjust = x < 0;
            return T(T(x + adjust) / D - adjust);
        }
        else {
            const bool adjust = x > 0;
            return T(T(x - adjust) / D - adjust);
        }
    }

    template <class T> static std::string describe(Params<T> p) {
        return p.isTrivial ? Trivial::describe<T>(p.trivial) : std::string{};
    }

    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        static_assert(D != 0, "divisor must be non-zero");
        if constexpr (Trivial::is_pow2(D)) {
            return Trivial::floordiv_ct<D>(x);
        }
        else if constexpr (D > 0) {
            const bool adjust = x < 0;
            return T(T(x + adjust) / T(D) - adjust);
        }
        else {
            const bool adjust = x > 0;
            return T(T(x - adjust) / T(D) - adjust);
        }
    }
};

}  // namespace fastdiv
