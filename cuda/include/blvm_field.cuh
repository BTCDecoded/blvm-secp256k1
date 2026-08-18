/* BLVM secp256k1 prime field — 4×u64 LE limbs. p = 2^256 - 2^32 - 977. */
#pragma once
#include "blvm_hd.h"
#include <stdint.h>

struct Fe {
    uint64_t d[4];
};

BLVM_HD_INLINE void fe_set0(Fe* r) { r->d[0] = r->d[1] = r->d[2] = r->d[3] = 0; }

BLVM_HD_INLINE int fe_is0(const Fe* a) {
    return (a->d[0] | a->d[1] | a->d[2] | a->d[3]) == 0;
}

BLVM_HD_INLINE int fe_eq(const Fe* a, const Fe* b) {
    return ((a->d[0] ^ b->d[0]) | (a->d[1] ^ b->d[1]) | (a->d[2] ^ b->d[2]) |
            (a->d[3] ^ b->d[3])) == 0;
}

BLVM_HD_INLINE void fe_from_be32(Fe* r, const uint8_t b[32]) {
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

BLVM_HD_INLINE void fe_to_be32(uint8_t b[32], const Fe* a) {
    for (int i = 0; i < 4; ++i) {
        uint64_t w = a->d[3 - i];
        int o = i * 8;
        b[o] = (uint8_t)(w >> 56);
        b[o + 1] = (uint8_t)(w >> 48);
        b[o + 2] = (uint8_t)(w >> 40);
        b[o + 3] = (uint8_t)(w >> 32);
        b[o + 4] = (uint8_t)(w >> 24);
        b[o + 5] = (uint8_t)(w >> 16);
        b[o + 6] = (uint8_t)(w >> 8);
        b[o + 7] = (uint8_t)w;
    }
}

BLVM_HD_INLINE void fe_p(uint64_t p[4]) {
    p[0] = 0xFFFFFFFEFFFFFC2FULL;
    p[1] = 0xFFFFFFFFFFFFFFFFULL;
    p[2] = 0xFFFFFFFFFFFFFFFFULL;
    p[3] = 0xFFFFFFFFFFFFFFFFULL;
}

BLVM_HD_INLINE int fe_cmp_raw(const uint64_t* a, const uint64_t* b) {
    for (int i = 3; i >= 0; --i) {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

/** Load wire encoding; return 0 if value >= p (BIP340 / compressed-pk reject). */
BLVM_HD_INLINE int fe_set_b32_limit(Fe* r, const uint8_t b[32]) {
    fe_from_be32(r, b);
    uint64_t p[4];
    fe_p(p);
    if (fe_cmp_raw(r->d, p) >= 0) return 0;
    return 1;
}

BLVM_HD_INLINE void fe_normalize(Fe* r) {
    uint64_t p[4];
    fe_p(p);
    if (fe_cmp_raw(r->d, p) < 0) return;
    unsigned __int128 borrow = 0;
    for (int i = 0; i < 4; ++i) {
        unsigned __int128 left = r->d[i];
        unsigned __int128 right = (unsigned __int128)p[i] + borrow;
        r->d[i] = (uint64_t)(left - right);
        borrow = left < right ? 1 : 0;
    }
}

#if defined(__CUDA_ARCH__)

__device__ __forceinline__ void fe_mul_comba32(const Fe* a, const Fe* b, uint32_t t32[16]) {
    uint32_t a32[8], b32[8];
#pragma unroll
    for (int i = 0; i < 4; i++) {
        a32[2 * i] = (uint32_t)(a->d[i]);
        a32[2 * i + 1] = (uint32_t)(a->d[i] >> 32);
        b32[2 * i] = (uint32_t)(b->d[i]);
        b32[2 * i + 1] = (uint32_t)(b->d[i] >> 32);
    }

    uint32_t r0 = 0, r1 = 0, r2 = 0;

#define MUL32_ACC(ai, bj)                       \
    asm volatile(                               \
        "mad.lo.cc.u32 %0, %3, %4, %0; \n\t"    \
        "madc.hi.cc.u32 %1, %3, %4, %1; \n\t"   \
        "addc.u32 %2, %2, 0; \n\t"              \
        : "+r"(r0), "+r"(r1), "+r"(r2)          \
        : "r"(a32[ai]), "r"(b32[bj]));

    MUL32_ACC(0, 0);
    t32[0] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 1);
    MUL32_ACC(1, 0);
    t32[1] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 2);
    MUL32_ACC(1, 1);
    MUL32_ACC(2, 0);
    t32[2] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 3);
    MUL32_ACC(1, 2);
    MUL32_ACC(2, 1);
    MUL32_ACC(3, 0);
    t32[3] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 4);
    MUL32_ACC(1, 3);
    MUL32_ACC(2, 2);
    MUL32_ACC(3, 1);
    MUL32_ACC(4, 0);
    t32[4] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 5);
    MUL32_ACC(1, 4);
    MUL32_ACC(2, 3);
    MUL32_ACC(3, 2);
    MUL32_ACC(4, 1);
    MUL32_ACC(5, 0);
    t32[5] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 6);
    MUL32_ACC(1, 5);
    MUL32_ACC(2, 4);
    MUL32_ACC(3, 3);
    MUL32_ACC(4, 2);
    MUL32_ACC(5, 1);
    MUL32_ACC(6, 0);
    t32[6] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(0, 7);
    MUL32_ACC(1, 6);
    MUL32_ACC(2, 5);
    MUL32_ACC(3, 4);
    MUL32_ACC(4, 3);
    MUL32_ACC(5, 2);
    MUL32_ACC(6, 1);
    MUL32_ACC(7, 0);
    t32[7] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(1, 7);
    MUL32_ACC(2, 6);
    MUL32_ACC(3, 5);
    MUL32_ACC(4, 4);
    MUL32_ACC(5, 3);
    MUL32_ACC(6, 2);
    MUL32_ACC(7, 1);
    t32[8] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(2, 7);
    MUL32_ACC(3, 6);
    MUL32_ACC(4, 5);
    MUL32_ACC(5, 4);
    MUL32_ACC(6, 3);
    MUL32_ACC(7, 2);
    t32[9] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(3, 7);
    MUL32_ACC(4, 6);
    MUL32_ACC(5, 5);
    MUL32_ACC(6, 4);
    MUL32_ACC(7, 3);
    t32[10] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(4, 7);
    MUL32_ACC(5, 6);
    MUL32_ACC(6, 5);
    MUL32_ACC(7, 4);
    t32[11] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(5, 7);
    MUL32_ACC(6, 6);
    MUL32_ACC(7, 5);
    t32[12] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(6, 7);
    MUL32_ACC(7, 6);
    t32[13] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    MUL32_ACC(7, 7);
    t32[14] = r0;
    t32[15] = r1;

#undef MUL32_ACC
}

__device__ __forceinline__ void fe_sqr_comba32(const Fe* a, uint32_t t32[16]) {
    uint32_t a32[8];
#pragma unroll
    for (int i = 0; i < 4; i++) {
        a32[2 * i] = (uint32_t)(a->d[i]);
        a32[2 * i + 1] = (uint32_t)(a->d[i] >> 32);
    }

    uint32_t r0 = 0, r1 = 0, r2 = 0;

#define SQR32_DIAG(ai)                          \
    asm volatile(                               \
        "mad.lo.cc.u32 %0, %3, %3, %0; \n\t"    \
        "madc.hi.cc.u32 %1, %3, %3, %1; \n\t"   \
        "addc.u32 %2, %2, 0; \n\t"              \
        : "+r"(r0), "+r"(r1), "+r"(r2)          \
        : "r"(a32[ai]));

#define SQR32_MUL2(ai, aj)                      \
    {                                           \
        uint32_t lo, hi;                        \
        asm volatile("mul.lo.u32 %0, %2, %3; \n\t" "mul.hi.u32 %1, %2, %3; \n\t" : "=r"(lo), "=r"(hi) \
                     : "r"(a32[ai]), "r"(a32[aj])); \
        asm volatile(                           \
            "add.cc.u32 %0, %0, %3; \n\t"       \
            "addc.cc.u32 %1, %1, %4; \n\t"      \
            "addc.u32 %2, %2, 0; \n\t"          \
            "add.cc.u32 %0, %0, %3; \n\t"       \
            "addc.cc.u32 %1, %1, %4; \n\t"      \
            "addc.u32 %2, %2, 0; \n\t"          \
            : "+r"(r0), "+r"(r1), "+r"(r2)      \
            : "r"(lo), "r"(hi));                \
    }

    SQR32_DIAG(0);
    t32[0] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 1);
    t32[1] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 2);
    SQR32_DIAG(1);
    t32[2] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 3);
    SQR32_MUL2(1, 2);
    t32[3] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 4);
    SQR32_MUL2(1, 3);
    SQR32_DIAG(2);
    t32[4] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 5);
    SQR32_MUL2(1, 4);
    SQR32_MUL2(2, 3);
    t32[5] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 6);
    SQR32_MUL2(1, 5);
    SQR32_MUL2(2, 4);
    SQR32_DIAG(3);
    t32[6] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(0, 7);
    SQR32_MUL2(1, 6);
    SQR32_MUL2(2, 5);
    SQR32_MUL2(3, 4);
    t32[7] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(1, 7);
    SQR32_MUL2(2, 6);
    SQR32_MUL2(3, 5);
    SQR32_DIAG(4);
    t32[8] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(2, 7);
    SQR32_MUL2(3, 6);
    SQR32_MUL2(4, 5);
    t32[9] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(3, 7);
    SQR32_MUL2(4, 6);
    SQR32_DIAG(5);
    t32[10] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(4, 7);
    SQR32_MUL2(5, 6);
    t32[11] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(5, 7);
    SQR32_DIAG(6);
    t32[12] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_MUL2(6, 7);
    t32[13] = r0;
    r0 = r1;
    r1 = r2;
    r2 = 0;
    SQR32_DIAG(7);
    t32[14] = r0;
    t32[15] = r1;

#undef SQR32_DIAG
#undef SQR32_MUL2
}

__device__ __forceinline__ void fe_reduce512_hybrid(uint32_t t32[16], Fe* r) {
    uint32_t t0 = t32[0], t1 = t32[1], t2 = t32[2], t3 = t32[3];
    uint32_t t4 = t32[4], t5 = t32[5], t6 = t32[6], t7 = t32[7];
    const uint32_t t8 = t32[8], t9 = t32[9], t10 = t32[10], t11 = t32[11];
    const uint32_t t12 = t32[12], t13 = t32[13], t14 = t32[14], t15 = t32[15];

    uint32_t a0, a1, a2, a3, a4, a5, a6, a7, a8;
    asm volatile(
        "mul.lo.u32 %0, %9, 977;\n\t"
        "mul.hi.u32 %1, %9, 977;\n\t"
        "mad.lo.cc.u32 %1, %10, 977, %1;\n\t"
        "madc.hi.u32 %2, %10, 977, 0;\n\t"
        "mad.lo.cc.u32 %2, %11, 977, %2;\n\t"
        "madc.hi.u32 %3, %11, 977, 0;\n\t"
        "mad.lo.cc.u32 %3, %12, 977, %3;\n\t"
        "madc.hi.u32 %4, %12, 977, 0;\n\t"
        "mad.lo.cc.u32 %4, %13, 977, %4;\n\t"
        "madc.hi.u32 %5, %13, 977, 0;\n\t"
        "mad.lo.cc.u32 %5, %14, 977, %5;\n\t"
        "madc.hi.u32 %6, %14, 977, 0;\n\t"
        "mad.lo.cc.u32 %6, %15, 977, %6;\n\t"
        "madc.hi.u32 %7, %15, 977, 0;\n\t"
        "mad.lo.cc.u32 %7, %16, 977, %7;\n\t"
        "madc.hi.u32 %8, %16, 977, 0;\n\t"
        : "=r"(a0), "=r"(a1), "=r"(a2), "=r"(a3), "=r"(a4), "=r"(a5), "=r"(a6), "=r"(a7), "=r"(a8)
        : "r"(t8), "r"(t9), "r"(t10), "r"(t11), "r"(t12), "r"(t13), "r"(t14), "r"(t15));

    uint32_t a9;
    asm volatile(
        "add.cc.u32 %0, %0, %9;\n\t"
        "addc.cc.u32 %1, %1, %10;\n\t"
        "addc.cc.u32 %2, %2, %11;\n\t"
        "addc.cc.u32 %3, %3, %12;\n\t"
        "addc.cc.u32 %4, %4, %13;\n\t"
        "addc.cc.u32 %5, %5, %14;\n\t"
        "addc.cc.u32 %6, %6, %15;\n\t"
        "addc.cc.u32 %7, %7, %16;\n\t"
        "addc.u32 %8, 0, 0;\n\t"
        : "+r"(a1), "+r"(a2), "+r"(a3), "+r"(a4), "+r"(a5), "+r"(a6), "+r"(a7), "+r"(a8), "=r"(a9)
        : "r"(t8), "r"(t9), "r"(t10), "r"(t11), "r"(t12), "r"(t13), "r"(t14), "r"(t15));

    uint32_t carry;
    asm volatile(
        "add.cc.u32 %0, %0, %9;\n\t"
        "addc.cc.u32 %1, %1, %10;\n\t"
        "addc.cc.u32 %2, %2, %11;\n\t"
        "addc.cc.u32 %3, %3, %12;\n\t"
        "addc.cc.u32 %4, %4, %13;\n\t"
        "addc.cc.u32 %5, %5, %14;\n\t"
        "addc.cc.u32 %6, %6, %15;\n\t"
        "addc.cc.u32 %7, %7, %16;\n\t"
        "addc.u32 %8, 0, 0;\n\t"
        : "+r"(t0), "+r"(t1), "+r"(t2), "+r"(t3), "+r"(t4), "+r"(t5), "+r"(t6), "+r"(t7), "=r"(carry)
        : "r"(a0), "r"(a1), "r"(a2), "r"(a3), "r"(a4), "r"(a5), "r"(a6), "r"(a7));

    uint32_t e_lo, e_carry;
    asm volatile("add.cc.u32 %0, %2, %3;\n\t" "addc.u32 %1, 0, 0;\n\t" : "=r"(e_lo), "=r"(e_carry)
                 : "r"(a8), "r"(carry));
    uint32_t e_hi = a9 + e_carry;

    uint32_t p_lo, p_hi;
    asm volatile("mul.lo.u32 %0, %2, 977;\n\t" "mul.hi.u32 %1, %2, 977;\n\t" : "=r"(p_lo), "=r"(p_hi)
                 : "r"(e_lo));
    uint32_t m1 = p_hi + e_hi * 977u;

    uint32_t ek0 = p_lo;
    uint32_t ek1, ek1_carry;
    asm volatile("add.cc.u32 %0, %2, %3;\n\t" "addc.u32 %1, 0, 0;\n\t" : "=r"(ek1), "=r"(ek1_carry)
                 : "r"(m1), "r"(e_lo));
    uint32_t ek2 = e_hi + ek1_carry;

    uint64_t rd0 = ((uint64_t)t1 << 32) | t0;
    uint64_t rd1 = ((uint64_t)t3 << 32) | t2;
    uint64_t rd2 = ((uint64_t)t5 << 32) | t4;
    uint64_t rd3 = ((uint64_t)t7 << 32) | t6;
    uint64_t ek_lo = ((uint64_t)ek1 << 32) | ek0;
    uint64_t ek_hi = (uint64_t)ek2;

    uint64_t c;
    asm volatile(
        "add.cc.u64 %0, %0, %5;\n\t"
        "addc.cc.u64 %1, %1, %6;\n\t"
        "addc.cc.u64 %2, %2, 0;\n\t"
        "addc.cc.u64 %3, %3, 0;\n\t"
        "addc.u64 %4, 0, 0;\n\t"
        : "+l"(rd0), "+l"(rd1), "+l"(rd2), "+l"(rd3), "=l"(c)
        : "l"(ek_lo), "l"(ek_hi));

    if (c) {
        asm volatile(
            "add.cc.u64 %0, %0, %4;\n\t"
            "addc.cc.u64 %1, %1, 0;\n\t"
            "addc.cc.u64 %2, %2, 0;\n\t"
            "addc.u64 %3, %3, 0;\n\t"
            : "+l"(rd0), "+l"(rd1), "+l"(rd2), "+l"(rd3)
            : "l"(0x1000003D1ULL));
    }

    uint64_t s0, s1, s2, s3, borrow;
    asm volatile(
        "sub.cc.u64 %0, %5, %9;\n\t"
        "subc.cc.u64 %1, %6, %10;\n\t"
        "subc.cc.u64 %2, %7, %11;\n\t"
        "subc.cc.u64 %3, %8, %12;\n\t"
        "subc.u64 %4, 0, 0;\n\t"
        : "=l"(s0), "=l"(s1), "=l"(s2), "=l"(s3), "=l"(borrow)
        : "l"(rd0), "l"(rd1), "l"(rd2), "l"(rd3), "l"(0xFFFFFFFEFFFFFC2FULL), "l"(0xFFFFFFFFFFFFFFFFULL),
          "l"(0xFFFFFFFFFFFFFFFFULL), "l"(0xFFFFFFFFFFFFFFFFULL));

    if (borrow == 0) {
        r->d[0] = s0;
        r->d[1] = s1;
        r->d[2] = s2;
        r->d[3] = s3;
    } else {
        r->d[0] = rd0;
        r->d[1] = rd1;
        r->d[2] = rd2;
        r->d[3] = rd3;
    }
}

#endif /* __CUDA_ARCH__ */

BLVM_HD_INLINE void fe_add(Fe* r, const Fe* a, const Fe* b) {
#if defined(__CUDA_ARCH__)
    uint64_t r0, r1, r2, r3, carry;
    asm volatile(
        "add.cc.u64 %0, %5, %9; \n\t"
        "addc.cc.u64 %1, %6, %10; \n\t"
        "addc.cc.u64 %2, %7, %11; \n\t"
        "addc.cc.u64 %3, %8, %12; \n\t"
        "addc.u64 %4, 0, 0; \n\t"
        : "=l"(r0), "=l"(r1), "=l"(r2), "=l"(r3), "=l"(carry)
        : "l"(a->d[0]), "l"(a->d[1]), "l"(a->d[2]), "l"(a->d[3]), "l"(b->d[0]), "l"(b->d[1]),
          "l"(b->d[2]), "l"(b->d[3]));

    uint64_t t0, t1, t2, t3, borrow;
    asm volatile(
        "sub.cc.u64 %0, %5, %9; \n\t"
        "subc.cc.u64 %1, %6, %10; \n\t"
        "subc.cc.u64 %2, %7, %11; \n\t"
        "subc.cc.u64 %3, %8, %12; \n\t"
        "subc.u64 %4, 0, 0; \n\t"
        : "=l"(t0), "=l"(t1), "=l"(t2), "=l"(t3), "=l"(borrow)
        : "l"(r0), "l"(r1), "l"(r2), "l"(r3), "l"(0xFFFFFFFEFFFFFC2FULL), "l"(0xFFFFFFFFFFFFFFFFULL),
          "l"(0xFFFFFFFFFFFFFFFFULL), "l"(0xFFFFFFFFFFFFFFFFULL));

    if (carry || borrow == 0) {
        r->d[0] = t0;
        r->d[1] = t1;
        r->d[2] = t2;
        r->d[3] = t3;
    } else {
        r->d[0] = r0;
        r->d[1] = r1;
        r->d[2] = r2;
        r->d[3] = r3;
    }
#else
    unsigned __int128 c = 0;
    for (int i = 0; i < 4; ++i) {
        c += (unsigned __int128)a->d[i] + b->d[i];
        r->d[i] = (uint64_t)c;
        c >>= 64;
    }
    if (c) {
        c = (unsigned __int128)r->d[0] + 0x1000003D1ULL;
        r->d[0] = (uint64_t)c;
        c >>= 64;
        for (int i = 1; i < 4; ++i) {
            c += r->d[i];
            r->d[i] = (uint64_t)c;
            c >>= 64;
        }
    }
    fe_normalize(r);
#endif
}

BLVM_HD_INLINE void fe_neg(Fe* r, const Fe* a) {
    if (fe_is0(a)) {
        fe_set0(r);
        return;
    }
    uint64_t p[4];
    fe_p(p);
    unsigned __int128 borrow = 0;
    for (int i = 0; i < 4; ++i) {
        unsigned __int128 left = p[i];
        unsigned __int128 right = (unsigned __int128)a->d[i] + borrow;
        r->d[i] = (uint64_t)(left - right);
        borrow = left < right ? 1 : 0;
    }
}

BLVM_HD_INLINE void fe_sub(Fe* r, const Fe* a, const Fe* b) {
#if defined(__CUDA_ARCH__)
    uint64_t r0, r1, r2, r3, borrow;
    asm volatile(
        "sub.cc.u64 %0, %5, %9; \n\t"
        "subc.cc.u64 %1, %6, %10; \n\t"
        "subc.cc.u64 %2, %7, %11; \n\t"
        "subc.cc.u64 %3, %8, %12; \n\t"
        "subc.u64 %4, 0, 0; \n\t"
        : "=l"(r0), "=l"(r1), "=l"(r2), "=l"(r3), "=l"(borrow)
        : "l"(a->d[0]), "l"(a->d[1]), "l"(a->d[2]), "l"(a->d[3]), "l"(b->d[0]), "l"(b->d[1]),
          "l"(b->d[2]), "l"(b->d[3]));

    if (borrow) {
        uint64_t a0, a1, a2, a3;
        asm volatile(
            "add.cc.u64 %0, %4, %8; \n\t"
            "addc.cc.u64 %1, %5, %9; \n\t"
            "addc.cc.u64 %2, %6, %10; \n\t"
            "addc.u64 %3, %7, %11; \n\t"
            : "=l"(a0), "=l"(a1), "=l"(a2), "=l"(a3)
            : "l"(r0), "l"(r1), "l"(r2), "l"(r3), "l"(0xFFFFFFFEFFFFFC2FULL), "l"(0xFFFFFFFFFFFFFFFFULL),
              "l"(0xFFFFFFFFFFFFFFFFULL), "l"(0xFFFFFFFFFFFFFFFFULL));
        r->d[0] = a0;
        r->d[1] = a1;
        r->d[2] = a2;
        r->d[3] = a3;
    } else {
        r->d[0] = r0;
        r->d[1] = r1;
        r->d[2] = r2;
        r->d[3] = r3;
    }
#else
    Fe t;
    fe_neg(&t, b);
    fe_add(r, a, &t);
#endif
}

BLVM_HD_INLINE void fe_reduce512(Fe* r, uint64_t limbs[8]) {
    const uint64_t R = 0x1000003D1ULL;
    for (int round = 0; round < 4; ++round) {
        unsigned __int128 c0 = (unsigned __int128)limbs[0] + (unsigned __int128)limbs[4] * R;
        unsigned __int128 c1 =
            (unsigned __int128)limbs[1] + (unsigned __int128)limbs[5] * R + (c0 >> 64);
        unsigned __int128 c2 =
            (unsigned __int128)limbs[2] + (unsigned __int128)limbs[6] * R + (c1 >> 64);
        unsigned __int128 c3 =
            (unsigned __int128)limbs[3] + (unsigned __int128)limbs[7] * R + (c2 >> 64);
        limbs[0] = (uint64_t)c0;
        limbs[1] = (uint64_t)c1;
        limbs[2] = (uint64_t)c2;
        limbs[3] = (uint64_t)c3;
        unsigned __int128 hi = c3 >> 64;
        limbs[4] = (uint64_t)hi;
        limbs[5] = 0;
        limbs[6] = 0;
        limbs[7] = 0;
        if (limbs[4] == 0) break;
    }
    r->d[0] = limbs[0];
    r->d[1] = limbs[1];
    r->d[2] = limbs[2];
    r->d[3] = limbs[3];
    // Phase 1c: one normalize suffices after reduce512 (was double).
    fe_normalize(r);
}

BLVM_HD_INLINE void fe_mul(Fe* r, const Fe* a, const Fe* b) {
#if defined(__CUDA_ARCH__)
    uint32_t t32[16];
    fe_mul_comba32(a, b, t32);
    fe_reduce512_hybrid(t32, r);
#else
    uint64_t limbs[8] = {0, 0, 0, 0, 0, 0, 0, 0};
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
    fe_reduce512(r, limbs);
#endif
}

BLVM_HD_INLINE void fe_sqr(Fe* r, const Fe* a) {
#if defined(__CUDA_ARCH__)
    uint32_t t32[16];
    fe_sqr_comba32(a, t32);
    fe_reduce512_hybrid(t32, r);
#else
    fe_mul(r, a, a);
#endif
}

BLVM_HD_INLINE void fe_sqr_n(Fe* a, int n) {
#pragma unroll 1
    for (int i = 0; i < n; ++i) {
        fe_sqr(a, a);
    }
}

BLVM_HD void fe_inv(Fe* r, const Fe* a) {
    /* Fermat a^(p-2): 255 sqr + 16 mul (Ultrafast-style addition chain). */
    Fe x0, x1, x2, x3, x4, x5, t;

    fe_sqr(&x0, a);
    fe_mul(&x0, &x0, a);

    fe_sqr(&x1, &x0);
    fe_mul(&x1, a, &x1);

    fe_sqr(&x2, &x1);
    fe_sqr(&x2, &x2);
    fe_sqr(&x2, &x2);
    fe_mul(&x2, &x1, &x2);

    fe_sqr(&x3, &x2);
    fe_sqr(&x3, &x3);
    fe_sqr(&x3, &x3);
    fe_mul(&x3, &x1, &x3);

    fe_sqr(&x4, &x3);
    fe_sqr(&x4, &x4);
    fe_mul(&x4, &x0, &x4);

    t = x4;
    fe_sqr_n(&t, 11);
    fe_mul(&x3, &t, &x4);

    t = x3;
    fe_sqr_n(&t, 22);
    fe_mul(&x4, &t, &x3);

    t = x4;
    fe_sqr_n(&t, 44);
    fe_mul(&x5, &t, &x4);

    t = x5;
    fe_sqr_n(&t, 88);
    fe_mul(&x5, &t, &x5);

    fe_sqr_n(&x5, 44);
    fe_mul(&x5, &x5, &x4);

    fe_sqr(&x5, &x5);
    fe_sqr(&x5, &x5);
    fe_sqr(&x5, &x5);
    fe_mul(&x5, &x1, &x5);

    fe_sqr(&t, &x5);
    fe_sqr_n(&t, 22);
    fe_mul(&t, &t, &x3);

    fe_sqr(&t, &t);
    fe_sqr(&t, &t);
    fe_sqr(&t, &t);
    fe_sqr(&t, &t);

    fe_sqr(&t, &t);
    fe_mul(&t, a, &t);
    fe_sqr(&t, &t);
    fe_sqr(&t, &t);
    fe_mul(&t, a, &t);
    fe_sqr(&t, &t);
    fe_mul(&t, a, &t);
    fe_sqr(&t, &t);
    fe_sqr(&t, &t);
    fe_mul(r, &t, a);
}

BLVM_HD_INLINE int fe_is_odd(const Fe* a) {
    Fe t = *a;
    fe_normalize(&t);
    return (int)(t.d[0] & 1ull);
}

BLVM_HD int fe_sqrt(Fe* r, const Fe* a) {
    /* a^((p+1)/4) — libsecp256k1 sliding-window chain */
    Fe x2, x3, x6, x9, x11, x22, x44, x88, x176, x220, x223, t;
    fe_sqr(&x2, a);
    fe_mul(&x2, &x2, a);
    fe_sqr(&x3, &x2);
    fe_mul(&x3, &x3, a);
    x6 = x3;
    for (int i = 0; i < 3; ++i) fe_sqr(&x6, &x6);
    fe_mul(&x6, &x6, &x3);
    x9 = x6;
    for (int i = 0; i < 3; ++i) fe_sqr(&x9, &x9);
    fe_mul(&x9, &x9, &x3);
    x11 = x9;
    for (int i = 0; i < 2; ++i) fe_sqr(&x11, &x11);
    fe_mul(&x11, &x11, &x2);
    x22 = x11;
    for (int i = 0; i < 11; ++i) fe_sqr(&x22, &x22);
    fe_mul(&x22, &x22, &x11);
    x44 = x22;
    for (int i = 0; i < 22; ++i) fe_sqr(&x44, &x44);
    fe_mul(&x44, &x44, &x22);
    x88 = x44;
    for (int i = 0; i < 44; ++i) fe_sqr(&x88, &x88);
    fe_mul(&x88, &x88, &x44);
    x176 = x88;
    for (int i = 0; i < 88; ++i) fe_sqr(&x176, &x176);
    fe_mul(&x176, &x176, &x88);
    x220 = x176;
    for (int i = 0; i < 44; ++i) fe_sqr(&x220, &x220);
    fe_mul(&x220, &x220, &x44);
    x223 = x220;
    for (int i = 0; i < 3; ++i) fe_sqr(&x223, &x223);
    fe_mul(&x223, &x223, &x3);
    t = x223;
    for (int i = 0; i < 23; ++i) fe_sqr(&t, &t);
    fe_mul(&t, &t, &x22);
    for (int i = 0; i < 6; ++i) fe_sqr(&t, &t);
    fe_mul(&t, &t, &x2);
    fe_sqr(&t, &t);
    fe_sqr(r, &t);
    Fe chk;
    fe_sqr(&chk, r);
    return fe_eq(&chk, a);
}
