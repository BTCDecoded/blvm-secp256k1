/* BLVM secp256k1 scalar mod n — 4×u64 LE limbs. */
#pragma once
#include "blvm_hd.h"
#include <stdint.h>

struct Sc {
    uint64_t d[4];
};

/* n = FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE BAAEDCE6 AF48A03B BFD25E8C D0364141 */
BLVM_HD_INLINE void sc_n(uint64_t n[4]) {
    n[0] = 0xBFD25E8CD0364141ULL;
    n[1] = 0xBAAEDCE6AF48A03BULL;
    n[2] = 0xFFFFFFFFFFFFFFFEULL;
    n[3] = 0xFFFFFFFFFFFFFFFFULL;
}

/* 2^256 - n */
BLVM_HD_INLINE void sc_k0(uint64_t k[4]) {
    k[0] = 0x402DA1732FC9BEBFULL;
    k[1] = 0x4551231950B75FC4ULL;
    k[2] = 0x0000000000000001ULL;
    k[3] = 0x0000000000000000ULL;
}

BLVM_HD_INLINE void sc_set0(Sc* r) { r->d[0] = r->d[1] = r->d[2] = r->d[3] = 0; }

BLVM_HD_INLINE int sc_is0(const Sc* a) {
    return (a->d[0] | a->d[1] | a->d[2] | a->d[3]) == 0;
}

BLVM_HD_INLINE int sc_eq(const Sc* a, const Sc* b) {
    return ((a->d[0] ^ b->d[0]) | (a->d[1] ^ b->d[1]) | (a->d[2] ^ b->d[2]) |
            (a->d[3] ^ b->d[3])) == 0;
}

BLVM_HD_INLINE int sc_cmp(const uint64_t* a, const uint64_t* b) {
    for (int i = 3; i >= 0; --i) {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

BLVM_HD_INLINE void sc_from_be32(Sc* r, const uint8_t b[32]) {
    r->d[3] = ((uint64_t)b[0] << 56) | ((uint64_t)b[1] << 48) | ((uint64_t)b[2] << 40) |
              ((uint64_t)b[3] << 32) | ((uint64_t)b[4] << 24) | ((uint64_t)b[5] << 16) |
              ((uint64_t)b[6] << 8) | (uint64_t)b[7];
    r->d[2] = ((uint64_t)b[8] << 56) | ((uint64_t)b[9] << 48) | ((uint64_t)b[10] << 40) |
              ((uint64_t)b[11] << 32) | ((uint64_t)b[12] << 24) | ((uint64_t)b[13] << 16) |
              ((uint64_t)b[14] << 8) | (uint64_t)b[15];
    r->d[1] = ((uint64_t)b[16] << 56) | ((uint64_t)b[17] << 48) | ((uint64_t)b[18] << 40) |
              ((uint64_t)b[19] << 32) | ((uint64_t)b[20] << 24) | ((uint64_t)b[21] << 16) |
              ((uint64_t)b[22] << 8) | (uint64_t)b[23];
    r->d[0] = ((uint64_t)b[24] << 56) | ((uint64_t)b[25] << 48) | ((uint64_t)b[26] << 40) |
              ((uint64_t)b[27] << 32) | ((uint64_t)b[28] << 24) | ((uint64_t)b[29] << 16) |
              ((uint64_t)b[30] << 8) | (uint64_t)b[31];
}

BLVM_HD_INLINE void sc_reduce(Sc* r) {
    uint64_t n[4];
    sc_n(n);
    while (sc_cmp(r->d, n) >= 0) {
        unsigned __int128 borrow = 0;
        for (int i = 0; i < 4; ++i) {
            unsigned __int128 left = r->d[i];
            unsigned __int128 right = (unsigned __int128)n[i] + borrow;
            r->d[i] = (uint64_t)(left - right);
            borrow = left < right ? 1 : 0;
        }
    }
}

/* Returns 1 if input was >= n (overflow), 0 otherwise. Always reduces into range. */
BLVM_HD_INLINE int sc_set_b32(Sc* r, const uint8_t b[32]) {
    sc_from_be32(r, b);
    uint64_t n[4];
    sc_n(n);
    int overflow = sc_cmp(r->d, n) >= 0;
    if (overflow) sc_reduce(r);
    return overflow;
}

BLVM_HD_INLINE int sc_is_high(const Sc* a) {
    /* n/2 */
    const uint64_t half[4] = {0xDFE92F46681B20A0ULL, 0x5D576E7357A4501DULL, 0xFFFFFFFFFFFFFFFFULL,
                              0x7FFFFFFFFFFFFFFFULL};
    return sc_cmp(a->d, half) > 0;
}

BLVM_HD_INLINE void sc_add(Sc* r, const Sc* a, const Sc* b) {
    unsigned __int128 c = 0;
    for (int i = 0; i < 4; ++i) {
        c += (unsigned __int128)a->d[i] + b->d[i];
        r->d[i] = (uint64_t)c;
        c >>= 64;
    }
    sc_reduce(r);
}

BLVM_HD_INLINE void sc_neg(Sc* r, const Sc* a) {
    if (sc_is0(a)) {
        sc_set0(r);
        return;
    }
    uint64_t n[4];
    sc_n(n);
    unsigned __int128 borrow = 0;
    for (int i = 0; i < 4; ++i) {
        unsigned __int128 left = n[i];
        unsigned __int128 right = (unsigned __int128)a->d[i] + borrow;
        r->d[i] = (uint64_t)(left - right);
        borrow = left < right ? 1 : 0;
    }
}

BLVM_HD_INLINE void sc_mul_512(uint64_t limbs[8], const Sc* a, const Sc* b) {
    for (int i = 0; i < 8; ++i) limbs[i] = 0;
    for (int i = 0; i < 4; ++i) {
        unsigned __int128 carry = 0;
        for (int j = 0; j < 4; ++j) {
            unsigned __int128 s =
                (unsigned __int128)limbs[i + j] + (unsigned __int128)a->d[i] * b->d[j] + carry;
            limbs[i + j] = (uint64_t)s;
            carry = s >> 64;
        }
        int k = i + 4;
        while (carry) {
            unsigned __int128 s = (unsigned __int128)limbs[k] + carry;
            limbs[k] = (uint64_t)s;
            carry = s >> 64;
            ++k;
        }
    }
}

BLVM_HD_INLINE void sc_mul(Sc* r, const Sc* a, const Sc* b) {
    uint64_t limbs[8];
    sc_mul_512(limbs, a, b);
    uint64_t k0[4];
    sc_k0(k0);
    /* Fold: r = lo + hi * (2^256 - n) mod n */
    for (int round = 0; round < 3; ++round) {
        if ((limbs[4] | limbs[5] | limbs[6] | limbs[7]) == 0) break;
        uint64_t mid[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        for (int i = 0; i < 4; ++i) {
            unsigned __int128 carry = 0;
            for (int j = 0; j < 4; ++j) {
                unsigned __int128 s = (unsigned __int128)mid[i + j] +
                                     (unsigned __int128)limbs[i + 4] * k0[j] + carry;
                mid[i + j] = (uint64_t)s;
                carry = s >> 64;
            }
            int kk = i + 4;
            while (carry) {
                unsigned __int128 s = (unsigned __int128)mid[kk] + carry;
                mid[kk] = (uint64_t)s;
                carry = s >> 64;
                ++kk;
            }
        }
        unsigned __int128 carry = 0;
        for (int i = 0; i < 4; ++i) {
            unsigned __int128 s = (unsigned __int128)limbs[i] + mid[i] + carry;
            limbs[i] = (uint64_t)s;
            carry = s >> 64;
        }
        for (int i = 0; i < 4; ++i) {
            unsigned __int128 s = (unsigned __int128)mid[i + 4] + carry;
            limbs[i + 4] = (uint64_t)s;
            carry = s >> 64;
        }
        if (carry) limbs[4] += (uint64_t)carry; /* next round folds */
    }
    r->d[0] = limbs[0];
    r->d[1] = limbs[1];
    r->d[2] = limbs[2];
    r->d[3] = limbs[3];
    sc_reduce(r);
}

/* Add 2^bit to r (no reduce). Used by sc_mul_shift_var. */
BLVM_HD_INLINE void sc_cadd_bit(Sc* r, unsigned bit, int flag) {
    if (!flag || bit >= 256) return;
    unsigned limb = bit >> 6;
    unsigned __int128 t = (unsigned __int128)r->d[0] + (limb == 0 ? (1ull << (bit & 63)) : 0);
    r->d[0] = (uint64_t)t;
    t >>= 64;
    t += (unsigned __int128)r->d[1] + (limb == 1 ? (1ull << (bit & 63)) : 0);
    r->d[1] = (uint64_t)t;
    t >>= 64;
    t += (unsigned __int128)r->d[2] + (limb == 2 ? (1ull << (bit & 63)) : 0);
    r->d[2] = (uint64_t)t;
    t >>= 64;
    t += (unsigned __int128)r->d[3] + (limb == 3 ? (1ull << (bit & 63)) : 0);
    r->d[3] = (uint64_t)t;
}

/* r = (a*b) >> shift, rounded (add bit shift-1). shift >= 256. Mirrors Rust scalar_mul_shift_var. */
BLVM_HD_INLINE void sc_mul_shift_var(Sc* r, const Sc* a, const Sc* b, unsigned shift) {
    uint64_t l[8];
    sc_mul_512(l, a, b);
    unsigned shiftlimbs = shift >> 6;
    unsigned shiftlow = shift & 63u;
    unsigned shifthigh = 64u - shiftlow;
    r->d[0] = shift < 512
                  ? ((l[shiftlimbs] >> shiftlow) |
                     ((shift < 448 && shiftlow != 0) ? (l[1 + shiftlimbs] << shifthigh) : 0))
                  : 0;
    r->d[1] = shift < 448
                  ? ((l[1 + shiftlimbs] >> shiftlow) |
                     ((shift < 384 && shiftlow != 0) ? (l[2 + shiftlimbs] << shifthigh) : 0))
                  : 0;
    r->d[2] = shift < 384
                  ? ((l[2 + shiftlimbs] >> shiftlow) |
                     ((shift < 320 && shiftlow != 0) ? (l[3 + shiftlimbs] << shifthigh) : 0))
                  : 0;
    r->d[3] = shift < 320 ? (l[3 + shiftlimbs] >> shiftlow) : 0;
    unsigned bitpos = shift - 1;
    int bit = (int)((l[bitpos >> 6] >> (bitpos & 63u)) & 1ull);
    sc_cadd_bit(r, 0, bit);
}

BLVM_HD_INLINE uint32_t sc_get_bits(const Sc* a, unsigned offset, unsigned count) {
    unsigned limb = offset >> 6;
    unsigned shift = offset & 63u;
    uint32_t mask = (count == 32) ? 0xffffffffu : ((1u << count) - 1u);
    if (((offset + count - 1) >> 6) == limb) {
        return ((uint32_t)(a->d[limb] >> shift)) & mask;
    }
    uint64_t lo = a->d[limb] >> shift;
    uint64_t hi = a->d[limb + 1] << (64u - shift);
    return (uint32_t)(lo | hi) & mask;
}

/* k = r1 + r2*lambda (mod n). Mirrors Rust Scalar::split_lambda. */
BLVM_HD_INLINE void sc_split_lambda(Sc* r1, Sc* r2, const Sc* k) {
    const Sc minus_b1 = {{0x6F547FA90ABFE4C3ULL, 0xE4437ED6010E8828ULL, 0, 0}};
    /* libsecp SECP256K1_SCALAR_CONST(FFF…, FFF…, FFF…, FFFFFFFE, 8A280AC5, …) */
    const Sc minus_b2 = {{0xD765CDA83DB1562CULL, 0x8A280AC50774346DULL, 0xFFFFFFFFFFFFFFFEULL,
                          0xFFFFFFFFFFFFFFFFULL}};
    const Sc g1 = {{0xE893209A45DBB031ULL, 0x3DAA8A1471E8CA7FULL, 0xE86C90E49284EB15ULL,
                    0x3086D221A7D46BCDULL}};
    const Sc g2 = {{0x1571B4AE8AC47F71ULL, 0x221208AC9DF506C6ULL, 0x6F547FA90ABFE4C4ULL,
                    0xE4437ED6010E8828ULL}};
    const Sc lambda = {{0xDF02967C1B23BD72ULL, 0x122E22EA20816678ULL, 0xA5261C028812645AULL,
                        0x5363AD4CC05C30E0ULL}};
    Sc c1, c2, t, neg;
    sc_mul_shift_var(&c1, k, &g1, 384);
    sc_mul_shift_var(&c2, k, &g2, 384);
    sc_mul(&t, &c1, &minus_b1);
    c1 = t;
    sc_mul(&t, &c2, &minus_b2);
    c2 = t;
    sc_add(r2, &c1, &c2);
    sc_mul(r1, r2, &lambda);
    sc_neg(&neg, r1);
    sc_add(r1, &neg, k);
}

BLVM_HD_INLINE void sc_split_128(Sc* r1, Sc* r2, const Sc* k) {
    r1->d[0] = k->d[0];
    r1->d[1] = k->d[1];
    r1->d[2] = 0;
    r1->d[3] = 0;
    r2->d[0] = k->d[2];
    r2->d[1] = k->d[3];
    r2->d[2] = 0;
    r2->d[3] = 0;
}

BLVM_HD_INLINE int sc_bit(const Sc* a, int i) {
    return (int)((a->d[i >> 6] >> (i & 63)) & 1ull);
}

BLVM_HD_INLINE int sc_bitlen(const Sc* a) {
    for (int limb = 3; limb >= 0; --limb) {
        if (a->d[limb] == 0) continue;
        uint64_t v = a->d[limb];
        for (int b = 63; b >= 0; --b) {
            if ((v >> b) & 1ull) return limb * 64 + b + 1;
        }
    }
    return 0;
}

/* k = k1 + k2*lambda (mod n); k1/k2 shortened; neg flags flip base Y. */
struct ScGlv {
    Sc k1, k2;
    int k1_neg, k2_neg;
};

#if defined(__CUDACC__)
__attribute__((noinline)) __attribute__((optnone))
#endif
BLVM_HD void sc_glv_decompose(ScGlv* r, const Sc* k) {
    const Sc lambda = {{0xDF02967C1B23BD72ULL, 0x122E22EA20816678ULL, 0xA5261C028812645AULL,
                        0x5363AD4CC05C30E0ULL}};
    Sc k2_mod, k2_neg_val, k1_mod, k1_neg_val, lk2, neg;
    sc_split_lambda(&k1_mod, &k2_mod, k);

    sc_neg(&k2_neg_val, &k2_mod);
    if (sc_bitlen(&k2_neg_val) < sc_bitlen(&k2_mod)) {
        r->k2 = k2_neg_val;
        r->k2_neg = 1;
    } else {
        r->k2 = k2_mod;
        r->k2_neg = 0;
    }

    sc_mul(&lk2, &k2_mod, &lambda);
    sc_neg(&neg, &lk2);
    sc_add(&k1_mod, &neg, k);

    sc_neg(&k1_neg_val, &k1_mod);
    if (sc_bitlen(&k1_neg_val) < sc_bitlen(&k1_mod)) {
        r->k1 = k1_neg_val;
        r->k1_neg = 1;
    } else {
        r->k1 = k1_mod;
        r->k1_neg = 0;
    }
}

BLVM_HD void sc_inv(Sc* r, const Sc* a) {
    /* a^(n-2) square-and-multiply */
    uint64_t e[4];
    sc_n(e);
    /* e = n - 2 */
    {
        unsigned __int128 borrow = 2;
        for (int i = 0; i < 4; ++i) {
            unsigned __int128 left = e[i];
            e[i] = (uint64_t)(left - borrow);
            borrow = left < borrow ? 1 : 0;
        }
    }
    Sc base = *a;
    sc_set0(r);
    r->d[0] = 1;
    for (int bit = 255; bit >= 0; --bit) {
        sc_mul(r, r, r);
        int limb = bit / 64;
        int b = bit % 64;
        if ((e[limb] >> b) & 1ull) sc_mul(r, r, &base);
    }
}
