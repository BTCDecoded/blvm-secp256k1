/* Strauss WNAF + GLV endomorphism for device ECDSA/Schnorr verify.
 * Phase 2 (default): Shamir-GLV 16-combo lincomb for u1*G+u2*Q (ECDSA verify).
 * Fallback: define BLVM_ECMULT_USE_STRAUSS for WINDOW_G PRE_G Strauss path.
 * Mirrors src/ecmult.rs ecmult_strauss:
 *   Q/A side: split_lambda + WINDOW_A=5 odd-multiples (affine + global Z)
 *   G side:   split_128 + WINDOW_G=15 + PRE_G / PRE_G_128 (Strauss only)
 * Fallback (null tables): G also uses WINDOW_A odd-multiples of G (host_kat / pre-upload).
 * WINDOW_A=6/TABLE_A=16 dens (H10) was REVERT — flat Msig/s, wall −2.3% vs H6. */
#pragma once
#include "blvm_point.cuh"

/* Default Shamir-GLV (16-combo). Force Strauss/PRE_G with -DBLVM_ECMULT_USE_STRAUSS.
 * Was held off for 497110 false-reject (nested P2WPKH gate bug — fixed 2026-07-30). */
#ifndef BLVM_ECMULT_USE_STRAUSS
#define BLVM_ECMULT_SHAMIR_GLV 1
#endif

enum {
    BLVM_WINDOW_A = 5,
    BLVM_TABLE_A = 8,
    BLVM_WINDOW_G = 15,
    BLVM_TABLE_G = 8192,
    BLVM_WNAF_SIZE = 129,
    BLVM_GE_STORAGE_BYTES = 64
};

/* Affine compact storage: Fe.d[4] LE limbs == Rust FeStorage.n[4]. */
struct GeStorage {
    Fe x;
    Fe y;
};

/* Device getters (defined in device.cu / host_kat stubs). Host path never calls these. */
#if defined(__CUDA_ARCH__)
__device__ const GeStorage* blvm_pre_g_get(void);
__device__ const GeStorage* blvm_pre_g_128_get(void);
#endif

/* beta: nontrivial cube root of 1 in Fp. lambda*(X:Y:Z) = (beta*X : Y : Z). */
BLVM_HD_INLINE void fe_const_beta(Fe* r) {
    r->d[0] = 0xC1396C28719501EEULL;
    r->d[1] = 0x9CF0497512F58995ULL;
    r->d[2] = 0x6E64479EAC3434E9ULL;
    r->d[3] = 0x7AE96A2B657C0710ULL;
}

/* Host ecmult_odd_multiples_table: pre[i] holds (2*i+1)*a in a shared-Z system;
 * zr[i] = Z ratios; *z = deferred global Z for the accumulator. */
BLVM_HD void ecmult_odd_multiples_table(Ge pre[BLVM_TABLE_A], Fe zr[BLVM_TABLE_A], Fe* z,
                                       const Gej* a) {
    Gej d;
    gej_double(&d, a);
    Ge d_ge;
    ge_set_xy(&d_ge, &d.x, &d.y);

    ge_set_gej_zinv(&pre[0], a, &d.z);
    Gej ai;
    gej_set_ge(&ai, &pre[0]);
    ai.z = a->z;
    zr[0] = d.z;

    for (int i = 1; i < BLVM_TABLE_A; ++i) {
        Gej sum;
        gej_add_ge_var(&sum, &ai, &d_ge, &zr[i]);
        ai = sum;
        ge_set_xy(&pre[i], &ai.x, &ai.y);
    }
    fe_mul(z, &ai.z, &d.z);
}

/* Bring pre[0..n) to the same global Z using zr ratios (host ge_table_set_globalz). */
BLVM_HD void ge_table_set_globalz(int n, Ge pre[BLVM_TABLE_A], const Fe zr[BLVM_TABLE_A]) {
    if (n <= 0) return;
    int i = n - 1;
    fe_normalize(&pre[i].y);
    Fe zs = zr[i];
    while (i > 0) {
        if (i != n - 1) {
            Fe t;
            fe_mul(&t, &zs, &zr[i]);
            zs = t;
        }
        --i;
        Ge ai = pre[i];
        ge_set_ge_zinv(&pre[i], &ai, &zs);
    }
}

/* Fill Jacobian pre[i] = (2*i+1)*a (fallback G table / legacy). */
BLVM_HD void ecmult_odd_multiples(Gej pre[BLVM_TABLE_A], const Gej* a) {
    pre[0] = *a;
    Gej d;
    gej_double(&d, a);
    for (int i = 1; i < BLVM_TABLE_A; ++i) {
        gej_add(&pre[i], &pre[i - 1], &d);
    }
}

/* WNAF: sum_i 2^i * wnaf[i] == a (mod n). Returns bits (= last_set_bit+1).
 * int16_t: WINDOW_G=15 digits reach ±(2^14-1). */
BLVM_HD int ecmult_wnaf(int16_t wnaf[BLVM_WNAF_SIZE], const Sc* a, int w) {
    Sc s = *a;
    int sign = 1;
    if (sc_get_bits(&s, 255, 1) != 0) {
        Sc tmp;
        sc_neg(&tmp, &s);
        s = tmp;
        sign = -1;
    }
    for (int i = 0; i < BLVM_WNAF_SIZE; ++i) wnaf[i] = 0;
    int bit = 0;
    int last_set = -1;
    int carry = 0;
    while (bit < BLVM_WNAF_SIZE) {
        if ((int)sc_get_bits(&s, (unsigned)bit, 1) == carry) {
            ++bit;
            continue;
        }
        int now = w;
        if (now > BLVM_WNAF_SIZE - bit) now = BLVM_WNAF_SIZE - bit;
        int word = (int)sc_get_bits(&s, (unsigned)bit, (unsigned)now) + carry;
        carry = (word >> (w - 1)) & 1;
        int wval = word - (carry << w);
        wnaf[bit] = (int16_t)(sign * wval);
        last_set = bit;
        bit += now;
    }
    return last_set + 1;
}

BLVM_HD_INLINE void ecmult_table_get_ge(Ge* r, const Ge* pre, int n) {
    int idx = ((n > 0 ? n : -n) - 1) / 2;
    *r = pre[idx];
    if (n < 0) fe_neg(&r->y, &r->y);
}

BLVM_HD_INLINE void ecmult_table_get_ge_lambda(Ge* r, const Ge* pre, const Fe* x_beta, int n) {
    int idx = ((n > 0 ? n : -n) - 1) / 2;
    ge_set_xy(r, &x_beta[idx], &pre[idx].y);
    if (n < 0) fe_neg(&r->y, &r->y);
}

BLVM_HD_INLINE void ecmult_table_get_ge_storage(Ge* r, const GeStorage* pre, int n) {
    int idx = ((n > 0 ? n : -n) - 1) / 2;
    ge_set_xy(r, &pre[idx].x, &pre[idx].y);
    if (n < 0) fe_neg(&r->y, &r->y);
}

BLVM_HD_INLINE void ecmult_table_get(Gej* r, const Gej* pre, int n) {
    int idx = ((n > 0 ? n : -n) - 1) / 2;
    *r = pre[idx];
    if (n < 0) fe_neg(&r->y, &r->y);
}

BLVM_HD_INLINE void ecmult_table_get_lambda(Gej* r, const Gej* pre, const Fe* beta, int n) {
    ecmult_table_get(r, pre, n);
    Fe bx;
    fe_mul(&bx, &r->x, beta);
    r->x = bx;
}

BLVM_HD_INLINE void ecmult_add_wnaf_ge(Gej* r, const Ge* pre, const Fe* x_beta, int16_t digit,
                                      int use_lambda) {
    if (digit == 0) return;
    Ge t;
    if (use_lambda)
        ecmult_table_get_ge_lambda(&t, pre, x_beta, (int)digit);
    else
        ecmult_table_get_ge(&t, pre, (int)digit);
    Gej sum;
    gej_add_ge_var(&sum, r, &t, nullptr);
    *r = sum;
}

/* Jacobian-table WNAF add (fallback G path without PRE_G). */
BLVM_HD_INLINE void ecmult_add_wnaf(Gej* r, const Gej* pre, const Fe* beta, int16_t digit,
                                    int use_lambda) {
    if (digit == 0) return;
    Gej t;
    if (use_lambda)
        ecmult_table_get_lambda(&t, pre, beta, (int)digit);
    else
        ecmult_table_get(&t, pre, (int)digit);
    Gej sum;
    gej_add(&sum, r, &t);
    *r = sum;
}

BLVM_HD_INLINE void ecmult_add_wnaf_storage(Gej* r, const GeStorage* pre, int16_t digit) {
    if (digit == 0) return;
    Ge t;
    ecmult_table_get_ge_storage(&t, pre, (int)digit);
    Gej sum;
    gej_add_ge_var(&sum, r, &t, nullptr);
    *r = sum;
}

/* PRE_G add into a deferred-Z accumulator (host add_zinv_var). */
BLVM_HD_INLINE void ecmult_add_wnaf_storage_zinv(Gej* r, const GeStorage* pre, int16_t digit,
                                                const Fe* z) {
    if (digit == 0) return;
    Ge t;
    ecmult_table_get_ge_storage(&t, pre, (int)digit);
    Gej sum;
    gej_add_zinv_var(&sum, r, &t, z);
    *r = sum;
}

BLVM_HD_INLINE int ecmult_have_pre_g(const GeStorage* pre_g, const GeStorage* pre_g_128) {
    return pre_g != nullptr && pre_g_128 != nullptr;
}

/* Build A-side affine odd-multiples + global Z + beta*x aux (host Strauss setup). */
BLVM_HD void ecmult_build_pre_a(Ge pre_a[BLVM_TABLE_A], Fe x_beta[BLVM_TABLE_A], Fe* z,
                                const Gej* a) {
    Fe zr[BLVM_TABLE_A];
    ecmult_odd_multiples_table(pre_a, zr, z, a);
    ge_table_set_globalz(BLVM_TABLE_A, pre_a, zr);
    Fe beta;
    fe_const_beta(&beta);
    for (int i = 0; i < BLVM_TABLE_A; ++i) {
        fe_mul(&x_beta[i], &pre_a[i].x, &beta);
    }
}

/* When mixing absolute-Z Jacobian G adds (no PRE_G), materialize true-affine A table. */
BLVM_HD void ecmult_pre_a_to_affine(Ge pre_a[BLVM_TABLE_A], Fe x_beta[BLVM_TABLE_A], Fe* z) {
    Fe zi, beta;
    fe_inv(&zi, z);
    for (int i = 0; i < BLVM_TABLE_A; ++i) {
        Ge ai = pre_a[i];
        ge_set_ge_zinv(&pre_a[i], &ai, &zi);
    }
    fe_set1(z);
    fe_const_beta(&beta);
    for (int i = 0; i < BLVM_TABLE_A; ++i) {
        fe_mul(&x_beta[i], &pre_a[i].x, &beta);
    }
}

BLVM_HD_INLINE void ecmult_apply_global_z(Gej* r, const Fe* z) {
    if (!r->infinity) {
        Fe t;
        fe_mul(&t, &r->z, z);
        r->z = t;
    }
}

/* r = na * a  via GLV + WNAF (no generator term). */
BLVM_HD void ecmult_glv(Gej* r, const Gej* a, const Sc* na) {
    if (a->infinity || sc_is0(na)) {
        gej_set_inf(r);
        return;
    }
    Sc na1, na_lam;
    sc_split_lambda(&na1, &na_lam, na);
    int16_t wnaf1[BLVM_WNAF_SIZE];
    int16_t wnaf_lam[BLVM_WNAF_SIZE];
    int bits1 = ecmult_wnaf(wnaf1, &na1, BLVM_WINDOW_A);
    int bits_lam = ecmult_wnaf(wnaf_lam, &na_lam, BLVM_WINDOW_A);
    int bits = bits1 > bits_lam ? bits1 : bits_lam;

    Ge pre[BLVM_TABLE_A];
    Fe x_beta[BLVM_TABLE_A];
    Fe z;
    ecmult_build_pre_a(pre, x_beta, &z, a);

    gej_set_inf(r);
    for (int i = bits - 1; i >= 0; --i) {
        Gej dbl;
        gej_double(&dbl, r);
        *r = dbl;
        if (i < bits1) ecmult_add_wnaf_ge(r, pre, x_beta, wnaf1[i], 0);
        if (i < bits_lam) ecmult_add_wnaf_ge(r, pre, x_beta, wnaf_lam[i], 1);
    }
    ecmult_apply_global_z(r, &z);
}

/* r = ng * G via split_128 + PRE_G / PRE_G_128 (host ecmult_gen_strauss). */
BLVM_HD void ecmult_gen_pre_g(Gej* r, const Sc* ng, const GeStorage* pre_g,
                              const GeStorage* pre_g_128) {
    Sc ng1, ng128;
    sc_split_128(&ng1, &ng128, ng);
    int16_t wnaf1[BLVM_WNAF_SIZE];
    int16_t wnaf128[BLVM_WNAF_SIZE];
    int bits1 = ecmult_wnaf(wnaf1, &ng1, BLVM_WINDOW_G);
    int bits128 = ecmult_wnaf(wnaf128, &ng128, BLVM_WINDOW_G);
    int bits = bits1 > bits128 ? bits1 : bits128;

    gej_set_inf(r);
    for (int i = bits - 1; i >= 0; --i) {
        Gej dbl;
        gej_double(&dbl, r);
        *r = dbl;
        if (i < bits1) ecmult_add_wnaf_storage(r, pre_g, wnaf1[i]);
        if (i < bits128) ecmult_add_wnaf_storage(r, pre_g_128, wnaf128[i]);
    }
}

/* r = ng * G + na * a  (ECDSA: ng=u1, na=u2, a=Q).
 * Device IBD path passes ctx PRE_G pointers (kernel args) — no device-symbol sync. */
BLVM_HD void ecmult_strauss(Gej* r, const Gej* a, const Sc* na, const Sc* ng,
                            const GeStorage* pre_g, const GeStorage* pre_g_128) {
    int have_a = !sc_is0(na) && !a->infinity;
    int have_g = ng != nullptr && !sc_is0(ng);
    int use_pre_g = ecmult_have_pre_g(pre_g, pre_g_128);

    if (!have_a && !have_g) {
        gej_set_inf(r);
        return;
    }
    if (!have_a) {
        if (use_pre_g) {
            ecmult_gen_pre_g(r, ng, pre_g, pre_g_128);
        } else {
            Gej g;
            gej_set_g(&g);
            ecmult_glv(r, &g, ng);
        }
        return;
    }
    if (!have_g) {
        ecmult_glv(r, a, na);
        return;
    }

    Sc na1, na_lam;
    sc_split_lambda(&na1, &na_lam, na);
    int16_t wnaf_na1[BLVM_WNAF_SIZE];
    int16_t wnaf_na_lam[BLVM_WNAF_SIZE];
    int bits_na1 = ecmult_wnaf(wnaf_na1, &na1, BLVM_WINDOW_A);
    int bits_na_lam = ecmult_wnaf(wnaf_na_lam, &na_lam, BLVM_WINDOW_A);
    int bits = bits_na1 > bits_na_lam ? bits_na1 : bits_na_lam;

    int16_t wnaf_ng1[BLVM_WNAF_SIZE];
    int16_t wnaf_ng_hi[BLVM_WNAF_SIZE];
    int bits_ng1 = 0;
    int bits_ng_hi = 0;
    Sc ng1, ng_hi;
    if (use_pre_g) {
        sc_split_128(&ng1, &ng_hi, ng);
        bits_ng1 = ecmult_wnaf(wnaf_ng1, &ng1, BLVM_WINDOW_G);
        bits_ng_hi = ecmult_wnaf(wnaf_ng_hi, &ng_hi, BLVM_WINDOW_G);
    } else {
        sc_split_lambda(&ng1, &ng_hi, ng);
        bits_ng1 = ecmult_wnaf(wnaf_ng1, &ng1, BLVM_WINDOW_A);
        bits_ng_hi = ecmult_wnaf(wnaf_ng_hi, &ng_hi, BLVM_WINDOW_A);
    }
    if (bits_ng1 > bits) bits = bits_ng1;
    if (bits_ng_hi > bits) bits = bits_ng_hi;

    Ge pre_a[BLVM_TABLE_A];
    Fe x_beta[BLVM_TABLE_A];
    Fe z;
    ecmult_build_pre_a(pre_a, x_beta, &z, a);

    Gej pre_g_fallback[BLVM_TABLE_A];
    if (!use_pre_g) {
        /* Absolute-Z Jacobian G adds: materialize true-affine A (one inv; host_kat only). */
        ecmult_pre_a_to_affine(pre_a, x_beta, &z);
        Gej g;
        gej_set_g(&g);
        ecmult_odd_multiples(pre_g_fallback, &g);
    }

    Fe beta;
    fe_const_beta(&beta);

    gej_set_inf(r);
    for (int i = bits - 1; i >= 0; --i) {
        Gej dbl;
        gej_double(&dbl, r);
        *r = dbl;
        if (i < bits_na1) ecmult_add_wnaf_ge(r, pre_a, x_beta, wnaf_na1[i], 0);
        if (i < bits_na_lam) ecmult_add_wnaf_ge(r, pre_a, x_beta, wnaf_na_lam[i], 1);
        if (use_pre_g) {
            if (i < bits_ng1) ecmult_add_wnaf_storage_zinv(r, pre_g, wnaf_ng1[i], &z);
            if (i < bits_ng_hi) ecmult_add_wnaf_storage_zinv(r, pre_g_128, wnaf_ng_hi[i], &z);
        } else {
            if (i < bits_ng1) ecmult_add_wnaf(r, pre_g_fallback, &beta, wnaf_ng1[i], 0);
            if (i < bits_ng_hi) ecmult_add_wnaf(r, pre_g_fallback, &beta, wnaf_ng_hi[i], 1);
        }
    }
    ecmult_apply_global_z(r, &z);
}

#if defined(BLVM_ECMULT_SHAMIR_GLV)

/* Mixed add r = a + b (affine); never in-place — gej_add_ge_var aliasing breaks doubles. */
#if defined(__CUDACC__)
__attribute__((noinline)) __attribute__((optnone))
#endif
BLVM_HD void ecmult_add_ge(Gej* r, const Gej* a, const Ge* b) {
    Gej sum;
    gej_add_ge_var(&sum, a, b, nullptr);
    *r = sum;
}

/* Affine add: Jacobian(Z=1) a + affine b → r. */
#if defined(__CUDACC__)
__attribute__((noinline)) __attribute__((optnone))
#endif
BLVM_HD void ecmult_jadd_ge(Gej* r, const Ge* a, const Ge* b) {
    Gej j;
    gej_set_ge(&j, a);
    ecmult_add_ge(r, &j, b);
}

/* r = ng*G + na*a via Shamir 4-way GLV interleave + 16-entry affine combo table.
 * P=G (generator), Q=a (pubkey). ~128 doublings + one mixed add per nonzero nibble. */
#if defined(__CUDACC__)
__attribute__((noinline)) __attribute__((optnone))
#endif
BLVM_HD void ecmult_shamir_glv(Gej* r, const Gej* q, const Sc* na, const Sc* ng) {
    int have_a = !sc_is0(na) && !q->infinity;
    int have_g = ng != nullptr && !sc_is0(ng);
    if (!have_a && !have_g) {
        gej_set_inf(r);
        return;
    }
    if (!have_a) {
        Gej g;
        gej_set_g(&g);
        ecmult_glv(r, &g, ng);
        return;
    }
    if (!have_g) {
        ecmult_glv(r, q, na);
        return;
    }

    Gej g_j;
    gej_set_g(&g_j);
    ScGlv da, db;
    sc_glv_decompose(&da, ng);
    sc_glv_decompose(&db, na);

    Ge aff_g, aff_q;
    gej_to_affine(&aff_g, &g_j);
    gej_to_affine(&aff_q, q);

    Fe beta;
    fe_const_beta(&beta);

    Ge pts[4];
    pts[0] = aff_g;
    if (da.k1_neg) fe_neg(&pts[0].y, &pts[0].y);
    fe_mul(&pts[1].x, &aff_g.x, &beta);
    pts[1].y = aff_g.y;
    if (da.k2_neg) fe_neg(&pts[1].y, &pts[1].y);
    pts[2] = aff_q;
    if (db.k1_neg) fe_neg(&pts[2].y, &pts[2].y);
    fe_mul(&pts[3].x, &aff_q.x, &beta);
    pts[3].y = aff_q.y;
    if (db.k2_neg) fe_neg(&pts[3].y, &pts[3].y);

    Ge table[16];
    table[1] = pts[0];
    table[2] = pts[1];
    table[4] = pts[2];
    table[8] = pts[3];

    Gej jc[11];
    ecmult_jadd_ge(&jc[0], &pts[0], &pts[1]); /* P1+P2  → table[3] */
    ecmult_jadd_ge(&jc[1], &pts[0], &pts[2]); /* P1+Q1  → table[5] */
    ecmult_jadd_ge(&jc[2], &pts[1], &pts[2]); /* P2+Q1  → table[6] */
    ecmult_jadd_ge(&jc[3], &pts[0], &pts[3]); /* P1+Q2  → table[9] */
    ecmult_jadd_ge(&jc[4], &pts[1], &pts[3]); /* P2+Q2  → table[10] */
    ecmult_jadd_ge(&jc[5], &pts[2], &pts[3]); /* Q1+Q2  → table[12] */
    ecmult_add_ge(&jc[6], &jc[0], &pts[2]); /* P1+P2+Q1  → table[7] */
    ecmult_add_ge(&jc[7], &jc[0], &pts[3]); /* P1+P2+Q2  → table[11] */
    ecmult_add_ge(&jc[8], &jc[1], &pts[3]); /* P1+Q1+Q2  → table[13] */
    ecmult_add_ge(&jc[9], &jc[2], &pts[3]); /* P2+Q1+Q2  → table[14] */
    ecmult_add_ge(&jc[10], &jc[6], &pts[3]); /* P1+P2+Q1+Q2 → table[15] */

    int degen = 0;
    for (int i = 0; i < 11; ++i) {
        if (jc[i].infinity) {
            degen = 1;
            break;
        }
    }

    if (degen) {
        int len1 = sc_bitlen(&da.k1);
        int len2 = sc_bitlen(&da.k2);
        int len3 = sc_bitlen(&db.k1);
        int len4 = sc_bitlen(&db.k2);
        int bits = len1;
        if (len2 > bits) bits = len2;
        if (len3 > bits) bits = len3;
        if (len4 > bits) bits = len4;
        gej_set_inf(r);
        for (int i = bits - 1; i >= 0; --i) {
            if (!r->infinity) {
                Gej dbl;
                gej_double(&dbl, r);
                *r = dbl;
            }
            for (int j = 0; j < 4; ++j) {
                const Sc* scs[4] = {&da.k1, &da.k2, &db.k1, &db.k2};
                if (!sc_bit(scs[j], i)) continue;
                if (r->infinity) {
                    gej_set_ge(r, &pts[j]);
                } else {
                    ecmult_add_ge(r, r, &pts[j]);
                }
            }
        }
        return;
    }

    Fe prefix[11];
    prefix[0] = jc[0].z;
    for (int i = 1; i < 11; ++i) fe_mul(&prefix[i], &prefix[i - 1], &jc[i].z);

    Fe inv_prod;
    fe_inv(&inv_prod, &prefix[10]);

    Fe z_inv[11];
    for (int i = 10; i > 0; --i) {
        fe_mul(&z_inv[i], &inv_prod, &prefix[i - 1]);
        fe_mul(&inv_prod, &inv_prod, &jc[i].z);
    }
    z_inv[0] = inv_prod;

    enum { T03 = 3, T05 = 5, T06 = 6, T09 = 9, T10 = 10, T12 = 12,
           T07 = 7, T11 = 11, T13 = 13, T14 = 14, T15 = 15 };
    const int tbl_map[11] = {T03, T05, T06, T09, T10, T12, T07, T11, T13, T14, T15};
    for (int i = 0; i < 11; ++i) {
        Fe zi2, zi3;
        fe_sqr(&zi2, &z_inv[i]);
        fe_mul(&zi3, &zi2, &z_inv[i]);
        fe_mul(&table[tbl_map[i]].x, &jc[i].x, &zi2);
        fe_mul(&table[tbl_map[i]].y, &jc[i].y, &zi3);
    }

    int len1 = sc_bitlen(&da.k1);
    int len2 = sc_bitlen(&da.k2);
    int len3 = sc_bitlen(&db.k1);
    int len4 = sc_bitlen(&db.k2);
    int bits = len1;
    if (len2 > bits) bits = len2;
    if (len3 > bits) bits = len3;
    if (len4 > bits) bits = len4;

    gej_set_inf(r);
    for (int i = bits - 1; i >= 0; --i) {
        if (!r->infinity) {
            Gej dbl;
            gej_double(&dbl, r);
            *r = dbl;
        }
        int idx = sc_bit(&da.k1, i) | (sc_bit(&da.k2, i) << 1) | (sc_bit(&db.k1, i) << 2) |
                  (sc_bit(&db.k2, i) << 3);
        if (idx == 0) continue;
        if (r->infinity) {
            gej_set_ge(r, &table[idx]);
        } else {
            ecmult_add_ge(r, r, &table[idx]);
        }
    }
}

#endif /* BLVM_ECMULT_SHAMIR_GLV */

BLVM_HD void gej_mul(Gej* r, const Gej* a, const Sc* k) { ecmult_glv(r, a, k); }

/* r = u1*G + u2*Q — device: Shamir-GLV 16-combo; host KAT: Strauss (nvcc -O3 HD quirk). */
BLVM_HD void gej_lincomb(Gej* r, const Sc* u1, const Sc* u2, const Gej* q,
                         const GeStorage* pre_g, const GeStorage* pre_g_128) {
#if defined(BLVM_ECMULT_SHAMIR_GLV)
#if defined(__CUDA_ARCH__)
    (void)pre_g;
    (void)pre_g_128;
    ecmult_shamir_glv(r, q, u2, u1);
#else
    ecmult_strauss(r, q, u2, u1, pre_g, pre_g_128);
#endif
#else
    ecmult_strauss(r, q, u2, u1, pre_g, pre_g_128);
#endif
}
