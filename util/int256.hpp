#pragma once
#include <compare>
#include <cstdint>
#include <type_traits>

// Minimal signed 256-bit integer (two's complement, little-endian limbs).
//
// Only exists so the 64-bit parameter solvers (Robison in particular, whose
// intermediate values reach ~2^(2W) * 2^W) can run in constexpr without
// overflowing __int128. Supports exactly the operations those solvers use.
// Not fast: division is bitwise long division. It is only ever used while
// computing parameters, never on the hot path.

#if defined(__SIZEOF_INT128__)
__extension__ typedef __int128          i128;
__extension__ typedef unsigned __int128 u128;
#else
// MSVC has no native 128-bit integer. Its STL ships constexpr class types with
// the full operator set (used internally by <ranges>), which is what the
// solvers need. They are slow at runtime (software multiply on ARM64), so the
// hot paths go through mul_wide()/umul_shr() in algorithms/base.hpp instead.
#include <__msvc_int128.hpp>
using i128 = std::_Signed128;
using u128 = std::_Unsigned128;
#endif

struct Int256 {
    uint64_t w[4] = {0, 0, 0, 0};

    constexpr Int256() = default;

    template <class T>
        requires(sizeof(T) <= 16 && !std::is_floating_point_v<T>)
    constexpr Int256(T v) {
        u128 u = static_cast<u128>(v);
        w[0] = static_cast<uint64_t>(u);
        w[1] = static_cast<uint64_t>(u >> 64);
        const bool neg = (T(-1) < T(0)) && (v < T(0));
        w[2] = w[3] = neg ? ~uint64_t(0) : 0;
    }

    // Truncating conversion to the low 128 bits.
    explicit constexpr operator i128() const {
        return static_cast<i128>((static_cast<u128>(w[1]) << 64) | w[0]);
    }

    constexpr bool is_neg() const { return (w[3] >> 63) != 0; }

    friend constexpr Int256 operator~(Int256 a) {
        for (auto& l : a.w) l = ~l;
        return a;
    }
    friend constexpr Int256 operator+(Int256 a, Int256 b) {
        Int256 r;
        u128 carry = 0;
        for (int i = 0; i < 4; i++) {
            u128 s = static_cast<u128>(a.w[i]) + b.w[i] + carry;
            r.w[i] = static_cast<uint64_t>(s);
            carry = s >> 64;
        }
        return r;
    }
    friend constexpr Int256 operator-(Int256 a) { return ~a + Int256(1); }
    friend constexpr Int256 operator-(Int256 a, Int256 b) { return a + (-b); }
    friend constexpr Int256 operator*(Int256 a, Int256 b) {
        Int256 r;
        for (int i = 0; i < 4; i++) {
            u128 carry = 0;
            for (int j = 0; i + j < 4; j++) {
                u128 cur = static_cast<u128>(a.w[i]) * b.w[j] + r.w[i + j] + carry;
                r.w[i + j] = static_cast<uint64_t>(cur);
                carry = cur >> 64;
            }
        }
        return r;
    }
    friend constexpr Int256 operator<<(Int256 a, int n) {
        Int256 r;
        const int q = n / 64, b = n % 64;
        for (int i = 3; i >= q; i--) {
            uint64_t v = a.w[i - q] << b;
            if (b && i - q - 1 >= 0) v |= a.w[i - q - 1] >> (64 - b);
            r.w[i] = v;
        }
        return r;
    }
    // Arithmetic right shift.
    friend constexpr Int256 operator>>(Int256 a, int n) {
        const uint64_t fill = a.is_neg() ? ~uint64_t(0) : 0;
        Int256 r;
        const int q = n / 64, b = n % 64;
        for (int i = 0; i < 4; i++) {
            uint64_t lo = (i + q < 4) ? a.w[i + q] : fill;
            uint64_t hi = (i + q + 1 < 4) ? a.w[i + q + 1] : fill;
            r.w[i] = b ? (lo >> b) | (hi << (64 - b)) : lo;
        }
        return r;
    }
    friend constexpr Int256 operator&(Int256 a, Int256 b) {
        for (int i = 0; i < 4; i++) a.w[i] &= b.w[i];
        return a;
    }

    friend constexpr bool operator==(const Int256& a, const Int256& b) {
        return a.w[0] == b.w[0] && a.w[1] == b.w[1] && a.w[2] == b.w[2] && a.w[3] == b.w[3];
    }
    friend constexpr std::strong_ordering operator<=>(const Int256& a, const Int256& b) {
        const int64_t ah = static_cast<int64_t>(a.w[3]), bh = static_cast<int64_t>(b.w[3]);
        if (ah != bh) return ah <=> bh;
        for (int i = 2; i >= 0; i--)
            if (a.w[i] != b.w[i]) return a.w[i] <=> b.w[i];
        return std::strong_ordering::equal;
    }

    // Truncating division (C semantics): quotient rounds toward zero, the
    // remainder takes the dividend's sign.
    static constexpr void divmod(Int256 n, Int256 d, Int256& q, Int256& r) {
        const bool nn = n.is_neg(), dn = d.is_neg();
        if (nn) n = -n;
        if (dn) d = -d;
        Int256 qq, rr;
        for (int i = 255; i >= 0; i--) {
            rr = rr << 1;
            rr.w[0] |= (n.w[i / 64] >> (i % 64)) & 1;
            if (!ult(rr, d)) {
                rr = rr - d;
                qq.w[i / 64] |= uint64_t(1) << (i % 64);
            }
        }
        q = (nn != dn) ? -qq : qq;
        r = nn ? -rr : rr;
    }
    friend constexpr Int256 operator/(Int256 a, Int256 b) {
        Int256 q, r;
        divmod(a, b, q, r);
        return q;
    }
    friend constexpr Int256 operator%(Int256 a, Int256 b) {
        Int256 q, r;
        divmod(a, b, q, r);
        return r;
    }

private:
    static constexpr bool ult(const Int256& a, const Int256& b) {
        for (int i = 3; i >= 0; i--)
            if (a.w[i] != b.w[i]) return a.w[i] < b.w[i];
        return false;
    }
};
