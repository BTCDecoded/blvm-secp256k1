/* BLVM ECDSA + BIP-340 Schnorr verify (one item) + batch kernels. */
#pragma once
#include "blvm_ecmult.cuh"
#include "blvm_sha256.cuh"

BLVM_HD int blvm_ecdsa_verify_one(const uint8_t msg32[32], const uint8_t pk33[33],
                                  const uint8_t sig64[64], const GeStorage* pre_g,
                                  const GeStorage* pre_g_128) {
    Sc r, s, z, sinv, u1, u2;
    if (sc_set_b32(&r, sig64)) return 0;
    if (sc_set_b32(&s, sig64 + 32)) return 0;
    if (sc_is0(&r) || sc_is0(&s) || sc_is_high(&s)) return 0;
    if (sc_set_b32(&z, msg32)) return 0;
    Gej q;
    if (!ge_from_compressed(&q, pk33)) return 0;
    sc_inv(&sinv, &s);
    sc_mul(&u1, &z, &sinv);
    sc_mul(&u2, &r, &sinv);
    Gej R;
    gej_lincomb(&R, &u1, &u2, &q, pre_g, pre_g_128);
    if (R.infinity) return 0;
    Fe rx;
    gej_get_x(&rx, &R);
    uint8_t rxb[32];
    fe_to_be32(rxb, &rx);
    Sc rxn;
    sc_set_b32(&rxn, rxb);
    return sc_eq(&rxn, &r);
}

BLVM_HD int ge_lift_x(Gej* r, const uint8_t x32[32]) {
    Fe x, x2, x3, y2, y, seven;
    if (!fe_set_b32_limit(&x, x32)) return 0;
    fe_sqr(&x2, &x);
    fe_mul(&x3, &x2, &x);
    fe_set0(&seven);
    seven.d[0] = 7;
    fe_add(&y2, &x3, &seven);
    if (!fe_sqrt(&y, &y2)) return 0;
    if (fe_is_odd(&y)) fe_neg(&y, &y);
    gej_set_xy(r, &x, &y);
    return 1;
}

BLVM_HD int blvm_schnorr_verify_one(const uint8_t msg32[32], const uint8_t pkx[32],
                                    const uint8_t sig64[64]) {
    Sc s;
    if (sc_set_b32(&s, sig64 + 32)) return 0;
    if (sc_is0(&s)) return 0;
    Gej P, R;
    if (!ge_lift_x(&P, pkx)) return 0;
    if (!ge_lift_x(&R, sig64)) return 0;
    const char tag[] = "BIP0340/challenge";
    uint8_t taghash[32];
    blvm_sha256_fixed((const uint8_t*)tag, 17, taghash);
    uint8_t buf[160];
    for (int i = 0; i < 32; ++i) {
        buf[i] = taghash[i];
        buf[32 + i] = taghash[i];
        buf[64 + i] = sig64[i];
        buf[96 + i] = pkx[i];
        buf[128 + i] = msg32[i];
    }
    uint8_t ehash[32];
    blvm_sha256_fixed(buf, 160, ehash);
    Sc e;
    sc_set_b32(&e, ehash);
    Gej sg, ep, rhs, G;
    gej_set_g(&G);
    gej_mul(&sg, &G, &s);
    gej_mul(&ep, &P, &e);
    gej_add(&rhs, &R, &ep);
    if (sg.infinity || rhs.infinity) return 0;
    Fe x1, x2, y1, y2, zi, zi2, zi3;
    fe_inv(&zi, &sg.z);
    fe_sqr(&zi2, &zi);
    fe_mul(&zi3, &zi2, &zi);
    fe_mul(&x1, &sg.x, &zi2);
    fe_mul(&y1, &sg.y, &zi3);
    fe_inv(&zi, &rhs.z);
    fe_sqr(&zi2, &zi);
    fe_mul(&zi3, &zi2, &zi);
    fe_mul(&x2, &rhs.x, &zi2);
    fe_mul(&y2, &rhs.y, &zi3);
    fe_normalize(&x1);
    fe_normalize(&x2);
    fe_normalize(&y1);
    fe_normalize(&y2);
    return fe_eq(&x1, &x2) && fe_eq(&y1, &y2);
}

#ifdef __CUDACC__
__global__ void blvm_ecdsa_verify_batch_kernel(const uint8_t* msgs, const uint8_t* pks,
                                               const uint8_t* sigs, uint8_t* out, size_t n,
                                               const GeStorage* pre_g,
                                               const GeStorage* pre_g_128) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = (uint8_t)blvm_ecdsa_verify_one(msgs + i * 32, pks + i * 33, sigs + i * 64, pre_g,
                                            pre_g_128);
}

__global__ void blvm_schnorr_verify_batch_kernel(const uint8_t* msgs, const uint8_t* pks,
                                                 const uint8_t* sigs, uint8_t* out, size_t n) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = (uint8_t)blvm_schnorr_verify_one(msgs + i * 32, pks + i * 32, sigs + i * 64);
}
#endif
