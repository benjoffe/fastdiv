#pragma once
// Arch D. Robison, "N-Bit Unsigned Division Via N-Bit Multiply-Add"
// (ARITH17, 2005) -- Algorithm 3, signed-numerator floor variant:
//
//   floor(x / d) == (A*x + B) >> S      for every signed W-bit x
//
// The smallest working shift is searched directly (starting at S = W), using
// the exact valid-B interval for each candidate A in {floor(2^S/d), +1}.
// A can need up to W+1 bits. Two forms handle that:
//
//   RobisonV1  the "Add-Back" form, defined below: keeps the multiply at native
//              width and adds k*x back to the high half. The 128-bit add of B
//              is left to the compiler.
//   RobisonV2  8/16/32-bit: the multiply-add is done in the double-width type
//              (the "Wider Multiply" form). 64-bit: A is split into a native
//              signed An plus k*2^64, and the carry out of the low word when
//              adding B is computed explicitly (see mul_add_64), so no 128-bit
//              addition is needed.
//
// RobisonV1 reuses RobisonV2's parameter solver, so RobisonV2 is defined first.
#include "algorithms/base.hpp"
#include "algorithms/trivial.hpp"

namespace fastdiv {

// Solver precision: must hold 2^S * (small factor) and products like A*r.
template <int W> struct RobisonSolverWide { using type = i128; };
template <> struct RobisonSolverWide<8>  { using type = int64_t; };
template <> struct RobisonSolverWide<16> { using type = int64_t; };
template <> struct RobisonSolverWide<64> { using type = Int256; };

struct RobisonV2 {
    static constexpr const char* name = "robison-v2";

    template <class T> struct Params {
        T divisor;
        bool isTrivial;  // +-2^k: handled by Trivial
        Trivial::Params<T> trivial;
        int shift;
        wide_t<T> A;  // carries the divisor's sign
        wide_t<T> B;
        T An;         // A wrapped to a signed W-bit value (used by the 64-bit path)
        int k;        // (A - An) / 2^W, in {-1, 0, 1}; 0 when A already fits
    };

private:
    template <class SW> struct Window { SW q, r; };

    template <class SW> static constexpr Window<SW> make_window(SW xEnd, SW d) {
        SW q = xEnd / d, r = xEnd % d;
        // floor division
        if (r < SW(0)) {
            r = r + d;
            q = q - SW(1);
        }
        return {q, r};
    }

    // Min/max of f = q*k - A*r over a window of d consecutive values ending at xEnd.
    template <class SW>
    static constexpr void window_min_max(const Window<SW>& w, SW d, SW A, SW k, SW& lo, SW& hi) {
        lo = w.q * k - A * w.r;
        hi = w.q * k;
        if (w.r < d - SW(1)) {
            const SW qPrev = w.q - SW(1);
            const SW loPrev = qPrev * k - A * (d - SW(1));
            const SW hiPrev = qPrev * k - A * (w.r + SW(1));
            if (loPrev < lo) {
                lo = loPrev;
            }
            if (hiPrev > hi) {
                hi = hiPrev;
            }
        }
    }

    // Valid B in [bMin, bMaxExcl) for (A, shift); false if empty.
    template <class SW>
    static constexpr bool valid_b_range(const Window<SW>& wR, const Window<SW>& wL, SW d, SW A,
                                        int shift, SW& bMin, SW& bMaxExcl) {
        const SW two_s = SW(1) << shift;
        const SW k = two_s - A * d;
        SW rMin, rMax, lMin, lMax;
        window_min_max(wR, d, A, k, rMin, rMax);
        window_min_max(wL, d, A, k, lMin, lMax);
        bMin = rMax > lMax ? rMax : lMax;
        bMaxExcl = (rMin < lMin ? rMin : lMin) + two_s;
        return bMin < bMaxExcl;
    }

public:
    template <class T> static constexpr Params<T> params(T d) {
        constexpr int W = traits<T>::bits;
        using SW = typename RobisonSolverWide<W>::type;
        using Wd = wide_t<T>;
        Params<T> p{};
        p.divisor = d;

        const bool neg = d < 0;
        const uwide_t<T> dabsU = abs_wide<T>(d);

        if (Trivial::is_pow2(d)) {
            p.isTrivial = true;
            p.trivial = Trivial::params<T>(d);
            return p;
        }

        const SW dabs = SW(dabsU);
        const SW xlo = -(SW(1) << (W - 1));
        const SW xhi = (SW(1) << (W - 1)) - SW(1);
        const Window<SW> wR = make_window<SW>(xhi, dabs);
        const Window<SW> wL = make_window<SW>(xlo + dabs - SW(1), dabs);

        SW A = 0, B = 0, bLo = 0, bHi = 0;
        int shift = W;
        for (;; shift++) {
            if (shift > 2 * W + 2) {
                FASTDIV_TRAP();  // never happens
            }
            const SW floorA = (SW(1) << shift) / dabs;
            bool found = false;
            for (int ceil = 0; ceil < 2 && !found; ceil++) {
                const SW a = floorA + SW(ceil);
                const SW b = a / SW(2);
                SW lo, hi;
                if (valid_b_range<SW>(wR, wL, dabs, a, shift, lo, hi) && b >= lo && b < hi) {
                    A = a; B = b; bLo = lo; bHi = hi; found = true;
                }
            }
            if (found) {
                break;
            }
        }

        // 32/64-bit: prefer a B with fewer non-zero 16-bit chunks (cheaper to
        // materialise on ARM): clear low chunks while B stays in range.
        if (W >= 32 && B != SW(0)) {
            int bitLen = 0;
            for (SW v = B; v != SW(0); v = v >> 1) {
                bitLen++;
            }
            for (int cleared = (bitLen - 1) / 16 * 16; cleared > 0; cleared -= 16) {
                const SW cand = (B >> cleared) << cleared;
                if (cand >= bLo && cand < bHi) {
                    B = cand;
                    break;
                }
            }
        }

        p.shift = shift;
        p.A = neg ? Wd(-A) : Wd(A);
        p.B = Wd(B);
        p.An = T(p.A);
        p.k = int((p.A - Wd(p.An)) >> W);
        return p;
    }

    template <class T>
    FASTDIV_INLINE static inline T floordiv(T x, Params<T> p) {
        using Wd = wide_t<T>;
        if (p.isTrivial) {
            return Trivial::floordiv<T>(x, p.trivial);
        }
        if constexpr (sizeof(T) == 8) {
            return mul_add_64<T>(x, p);
        }
        else {
            return T((Wd(x) * p.A + p.B) >> p.shift);
        }
    }

    // (x*A + B) >> S for W = 64, where the compiler's generic i128 signed*
    // unsigned multiply is poor when A >= 2^63 (three multiplies). Splitting
    // A = An + k*2^64 with An a native signed value gives one smulh + one mul
    // for any A, and the k*x term is a single add/sub:
    //   (x*A + B) >> 64  ==  smulh(x, An) + k*x + carry(lo(x*An) + B)
    // then the remaining S-64 bits are a plain arithmetic shift. Exact (floor
    // of floor), and for cooperative divisors k == 0 so it is just
    // smulh + carry. B < 2^63 so the carry only comes from the low word.
    template <class T>
    FASTDIV_INLINE static inline T mul_add_64(T x, Params<T> p) {
        // Hide An's sign: with a known-negative constant, clang rewrites
        // smulh(x, An) as umulh + madd corrections inside loops. The barrier is
        // on the loop-invariant An, so it hoists out of loops.
        T An = p.An;
        FASTDIV_OPAQUE(An);
        const i128 m = mul_wide<T>(x, An);
        const uint64_t lo = uint64_t(m);
        int64_t t = int64_t(m >> 64);
        if (p.k > 0) {
            t += x;
        }
        else if (p.k < 0) {
            t = int64_t(uint64_t(t) - uint64_t(x));
        }
        t += (lo + uint64_t(p.B) < lo);
        return T(t >> (p.shift - 64));
    }

    template <class T> static std::string describe(Params<T> p) {
        return p.isTrivial ? Trivial::describe<T>(p.trivial) : std::string{};
    }

    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        static_assert(D != 0, "divisor must be non-zero");
        constexpr Params<T> p = params<T>(D);
        using Wd = wide_t<T>;
        if constexpr (p.isTrivial) {
            return Trivial::floordiv_ct<D>(x);
        }
        else if constexpr (sizeof(T) == 8) {
            return mul_add_64<T>(x, p);
        }
        else {
            return T((Wd(x) * p.A + p.B) >> p.shift);
        }
    }
};

// ---------------------------------------------------------------------------
// RobisonV1: the "Add-Back" form.
//
// When A needs more than a native signed W-bit register (|A| > INT_MAX), keep
// the multiply native instead of widening to 2W. With T the true A and
// A_n = T wrapped to a signed W-bit value, T = A_n + k * 2^W (k is -1, 0 or +1):
//
//   (T*x + B) >> W  ==  ((A_n*x + B) >> W) + k*x
//
// because k*2^W*x is an exact multiple of 2^W. The first term is the "multiply,
// keep the high half" a native-width ISA gives you directly; the correction is
// one add or subtract of x. The result is then shifted the remaining S - W bits.
// For cooperative divisors (A already native, k == 0) this is the same as
// RobisonV2 and uses a plain multiply-add.
// ---------------------------------------------------------------------------
struct RobisonV1 {
    static constexpr const char* name = "robison-v1";

    template <class T> struct Params {
        T divisor;
        bool isTrivial;  // +-2^k: handled by Trivial
        Trivial::Params<T> trivial;
        int shift;
        wide_t<T> A;  // k == 0: the full multiplier
        wide_t<T> B;
        T An;         // k != 0: A wrapped to a native signed W-bit value
        int k;        // (A - An) / 2^W, in {-1, 0, 1}; non-zero means add-back
    };

private:
    // k == 0: (x*A + B) >> S with A known to fit T. Identical to the plain
    // wide expression; on MSVC avoid i128's full 128x128 operator*. The
    // 64-bit path applies the same barrier as the k != 0 case: when A is a
    // visible compile-time constant, clang rewrites smulh(x, A) into
    // umulh + sign-correction, adding instructions to every iteration.
    template <class T>
    FASTDIV_INLINE static inline T mul_add_k0(T x, const Params<T>& p) {
        if constexpr (sizeof(T) == 8) {
            T An = T(p.A);
            FASTDIV_OPAQUE(An);
            return T((mul_wide<T>(x, An) + p.B) >> p.shift);
        } else {
            return T((wide_t<T>(x) * p.A + p.B) >> p.shift);
        }
    }

public:
    template <class T> static constexpr Params<T> params(T d) {
        constexpr int W = traits<T>::bits;
        using Wd = wide_t<T>;
        const auto r = RobisonV2::params<T>(d);
        Params<T> p{};
        p.divisor = r.divisor;
        p.isTrivial = r.isTrivial;
        p.trivial = r.trivial;
        p.shift = r.shift;
        p.A = r.A;
        p.B = r.B;
        if (!r.isTrivial) {
            p.An = T(r.A);  // two's-complement wrap to W bits
            p.k = int((r.A - Wd(p.An)) >> W);
        }
        return p;
    }

    template <class T>
    FASTDIV_INLINE static inline T floordiv(T x, Params<T> p) {
        constexpr int W = traits<T>::bits;
        using U = uns_t<T>;
        if (p.isTrivial) {
            return Trivial::floordiv<T>(x, p.trivial);
        }
        if (p.k == 0) {
            return mul_add_k0<T>(x, p);
        }
        T An = p.An;
        // see RobisonV2::mul_add_64
        if constexpr (sizeof(T) == 8) {
            FASTDIV_OPAQUE(An);
        }
        const T h = T((mul_wide<T>(x, An) + p.B) >> W);
        const T r = T(U(h) + U(U(p.k) * U(x)));  // h + k*x, wrap-safe
        return T(r >> (p.shift - W));
    }

    template <class T> static std::string describe(Params<T> p) {
        return p.isTrivial ? Trivial::describe<T>(p.trivial) : std::string{};
    }

    template <auto D> FASTDIV_INLINE static inline decltype(D) floordiv_ct(decltype(D) x) {
        using T = decltype(D);
        static_assert(D != 0, "divisor must be non-zero");
        constexpr Params<T> p = params<T>(D);
        constexpr int W = traits<T>::bits;
        using U = uns_t<T>;
        if constexpr (p.isTrivial) {
            return Trivial::floordiv_ct<D>(x);
        }
        else if constexpr (p.k == 0) {
            return mul_add_k0<T>(x, p);
        }
        else {  // add-back
            T An = p.An;
            if constexpr (sizeof(T) == 8) {
                FASTDIV_OPAQUE(An);
            }
            const T h = T((mul_wide<T>(x, An) + p.B) >> W);
            const T r = T(U(h) + U(U(p.k) * U(x)));
            return T(r >> (p.shift - W));
        }
    }
};

}  // namespace fastdiv
