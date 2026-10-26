#pragma once
// The "ubiquitous" divmod floor-division (the baseline everything else beats).
#include "algorithms/base.hpp"
#include "algorithms/trivial.hpp"

// ---------------------------------------------------------------------------
// The traditional solution: divmod, then fix up the quotient when the remainder is non-zero and the
// operands have opposite signs.
// ---------------------------------------------------------------------------
namespace fastdiv {

struct Traditional {
    static constexpr const char* name = "traditional";

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

        const T y = p.divisor;
        const T div = T(x / y);
        const T rem = T(x % y);
        return T(div - (rem != 0 && (x < 0) != (y < 0)));

        // The form below produces strictly better codegen on clang -- it
        // avoids the add-back fixup on divisors whose magic multiplier is
        // >= 2^(W-1) (e.g. 15), and keeps the quotient on a single dependency
        // chain for negative divisors. It is not used here because it is not
        // the shape production code is typically written in, and
        // `traditional` is meant to represent the realistic baseline.
        //
        // constexpr int W = traits<T>::bits;
        // const T rem = T(x - div * y);  // same value but often compiles better
        // if (y < 0) return T(div - (rem > 0));
        // else       return T(div + (rem >> (W - 1)));
    }

    template <class T> static std::string describe(Params<T> p) {
        return p.isTrivial ? Trivial::describe<T>(p.trivial) : std::string{};
    }

    // Body inlined (not `floordiv(x, params(D))`) because with -DNDEBUG some
    // Clang pass leaves the divisor as a register-held constant through the
    // Params-struct call, which blocks InstCombine's divide-by-constant
    // strength reduction and leaves a real sdiv in the hot loop.
    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        static_assert(D != 0, "divisor must be non-zero");
        if constexpr (Trivial::is_pow2(D)) {
            return Trivial::floordiv_ct<D>(x);
        }
        const T div = T(x / T(D));
        const T rem = T(x % T(D));
        return T(div - (rem != 0 && (x < 0) != (D < 0)));

        // Better codegen (see note on floordiv above), but not the shape
        // production code is written in -- kept out for a realistic baseline:
        //
        // const T rem = T(x - div * T(D));
        // if constexpr (D < 0) return T(div - (rem > 0));
        // else                 return T(div + (rem >> (traits<T>::bits - 1)));
    }
};

}  // namespace fastdiv
