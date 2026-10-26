#pragma once
// Ben Joffe's floor division by a constant.
//
// Two candidate (multiplier, shift) pairs are searched: round-up
// (M = ceil(2^S/|D|)) and round-down (M = floor(2^S/|D|)), each with a
// one-compare bias on x that makes it exact over the full signed range:
//
//   round-up    t = x + (x < -1)                        (D > 0)
//   round-down  t = x + (uint(x) < 2^(W-1) - 1)         (D > 0)
//   result      (t * M) >> S
//
// When both directions work and neither lands on a "clean" shift of W bits,
// the sign-selected variant is used instead (no compare at all):
//
//   m = M + (x >> (W-1));   result = (x * m) >> S
#include "algorithms/base.hpp"
#include "algorithms/trivial.hpp"

namespace fastdiv {

struct BenJoffe {
    static constexpr const char* name = "benjoffe";

    static constexpr int ModeTrivial      = 1;  // +-2^k: handled by trivial.hpp
    static constexpr int ModeSignSelected = 2;
    static constexpr int ModeRoundUp      = 3;
    static constexpr int ModeRoundDown    = 4;

    template <class T> struct Params {
        T divisor;
        int mode;            // Mode*
        Trivial::Params<T> trivial;  // ModeTrivial only
        int shift;           // S
        T mul;               // signed multiplier (carries the divisor's sign)
        int upShift;         // smallest working round-up S, 0 if none (div_report)
        int downShift;       // smallest working round-down S, 0 if none (div_report)
    };

    template <class T> static constexpr Params<T> params(T d) {
        constexpr int W = traits<T>::bits;
        using UW = uwide_t<T>;
        using SWd = wide_t<T>;
        Params<T> p{};
        p.divisor = d;

        const bool neg = d < 0;
        const UW dabs = abs_wide<T>(d);
        const UW maxint = UW(1) << (W - 1);

        if (Trivial::is_pow2(d)) {
            p.mode = ModeTrivial;
            p.trivial = Trivial::params<T>(d);
            return p;
        }

        // Search both rounding directions for the smallest working shift S >= W.
        // With err = the multiplier's rounding error at 2^S, a direction works
        // when 2^S > bound * err (and M fits), where `bound` is the largest
        // numerator magnitude, snapped down to a multiple of |D|, at which the
        // error could flip the quotient.
        const UW top = maxint + neg;  // largest magnitude to cover
        const UW boundUp = top / dabs * dabs - 1;
        const UW downMult = maxint / dabs * dabs + 1;
        const UW boundDown = downMult - (neg && downMult == maxint);

        struct Found {
            bool pass = false;
            int S = 0;
            UW M = 0;
            int mode = 0;
        };
        const auto search = [&](bool up) {
            for (int s = W; s < 2 * W; s++) {
                const UW two_s = UW(1) << s;
                const UW rem = two_s % dabs;
                const UW M = two_s / dabs + up;
                const UW err = up ? dabs - rem : rem;
                if (two_s > (up ? boundUp : boundDown) * err && M < maxint) {
                    return Found{true, s, M, up ? ModeRoundUp : ModeRoundDown};
                }
            }
            return Found{};
        };
        const Found up = search(true), down = search(false);
        if (!up.pass && !down.pass) {
            FASTDIV_TRAP();  // never happens; fails constexpr if it does
        }

        // slightly different decision tree for x86:
        bool isX86 = false;
#if (defined(__i386__) || defined(__x86_64__) || defined(_M_IX86) || defined(_M_X64))
        isX86 = true;
#endif

        p.upShift = up.pass ? up.S : 0;
        p.downShift = down.pass ? down.S : 0;

        Found chosen;
        if (down.pass && !up.pass) {
            // Forced solution
            chosen = down;
        }
        else if (up.pass && !down.pass) {
            // Forced solution
            chosen = up;
        }
        else if (down.pass && W > 8 && up.S != W && down.S == W) {
            // Favour clean 64-bit shift with round-down,
            // wins on all tested machines
            chosen = down;
        }
        else if (up.pass && (W > 8 && up.S == W && down.S != W) && !isX86) {
            // Favour clean 64-bit shift with round-up,
            // except on x86, as testing indicates it either
            // ties with sign-selected, or sometimes loses to it
            chosen = up;
        }
        else {
            chosen.S = up.S > down.S ? up.S : down.S;
            chosen.M = (UW(1) << chosen.S) / dabs + !neg;
            chosen.mode = ModeSignSelected;
        }

        p.mode = chosen.mode;
        p.shift = chosen.S;
        p.mul = neg ? T(-SWd(chosen.M)) : T(SWd(chosen.M));
        return p;
    }

private:
    // Hide M's sign so clang can't prove m = M + (x >> W-1) >= 0
    // (it then swaps the signed multiply-high for umulh + madd).
    // The barrier is on the loop-invariant M, not on m, so it
    // hoists out of loops and doesn't block unrolling.
    template <class T>
    FASTDIV_INLINE static inline T sign_selected_div(T x, const Params<T>& p) {
        constexpr int W = traits<T>::bits;
        T M = p.mul;
        FASTDIV_OPAQUE(M);
        const T m = T(M + (x >> (W - 1)));
        return mul_shr<T>(x, m, p.shift);
    }

    // Round-up:   t = x + (x < -1)  (D > 0),  x - (x > 1)  (D < 0)
    // Round-down: t = x + (uint(x) < 2^(W-1) - 1)  (D > 0),  x - (uint(x) > 2^(W-1))  (D < 0)
    template <class T>
    FASTDIV_INLINE static inline T biased_div(T x, const Params<T>& p) {
        constexpr int W = traits<T>::bits;
        using U = uns_t<T>;
        const bool neg = p.divisor < 0;
        constexpr U maxint = U(1) << (W - 1);
        bool bit;
        if (p.mode == ModeRoundUp) {
            bit = neg ? (x > T(1)) : (x < T(-1));
        }
        else {
            bit = neg ? (U(x) > maxint) : (U(x) < maxint - 1);
        }
        // x - bit == (x - 1) + !bit. AArch64 has "increment if" (cinc)
        // but no "decrement if", so x - bit compiles to cmp+cset+sub in
        // series. Computing x-1 separately gives cmp+cinc on the
        // critical path (x-1 runs in parallel), same as x + bit. The
        // barrier stops the compiler folding the -1 back after the cinc.
        T t;
        if (!neg) {
            t = T(U(x) + U(bit));
        }
        else {
            U xm1 = U(U(x) - U(1));
            FASTDIV_OPAQUE(xm1);
            t = T(U(xm1 + U(!bit)));
        }
        // Hide the (positive) multiplier's sign too, otherwise clang
        // rewrites the signed multiply-high as umulh + correction.
        T M = p.mul;
        FASTDIV_OPAQUE(M);
        return mul_shr<T>(t, M, p.shift);
    }

public:
    template <class T>
    FASTDIV_INLINE static inline T floordiv(T x, Params<T> p) {
        if (p.mode == ModeTrivial) {
            return Trivial::floordiv<T>(x, p.trivial);
        }
        if (p.mode == ModeSignSelected) {
            return sign_selected_div<T>(x, p);
        }
        return biased_div<T>(x, p);
    }

    template <class T> static std::string describe(Params<T> p) {
        if (p.mode == ModeTrivial) {
            return Trivial::describe<T>(p.trivial);
        }
        const std::string ms = "M=" + std::to_string((long long)p.mul) + "  S=" + std::to_string(p.shift);
        if (p.mode == ModeSignSelected) {
            return "sign-selected  " + ms;
        }
        return std::string(p.mode == ModeRoundUp ? "round-up  " : "round-down  ") + ms;
    }

    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        static_assert(D != 0, "divisor must be non-zero");
        constexpr Params<T> p = params<T>(D);
        if constexpr (p.mode == ModeTrivial) {
            return Trivial::floordiv_ct<D>(x);
        }
        else if constexpr (p.mode == ModeSignSelected) {
            return sign_selected_div<T>(x, p);
        }
        else {
            return biased_div<T>(x, p);
        }
    }
};

}  // namespace fastdiv
