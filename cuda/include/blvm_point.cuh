/* BLVM Jacobian points (secp256k1 a=0) + scalar mult. */
#pragma once
#include "blvm_field.cuh"
#include "blvm_scalar.cuh"

/* Affine point (x, y) or infinity. */
struct Ge {
    Fe x, y;
    int infinity;
};

struct Gej {
    Fe x, y, z;
    int infinity;
};

BLVM_HD_INLINE void fe_set1(Fe* r) {
    fe_set0(r);
    r->d[0] = 1;
}

BLVM_HD_INLINE void ge_set_xy(Ge* r, const Fe* x, const Fe* y) {
    r->x = *x;
    r->y = *y;
    r->infinity = 0;
}

BLVM_HD_INLINE void ge_set_inf(Ge* r) {
    fe_set0(&r->x);
    fe_set0(&r->y);
    r->infinity = 1;
}

/* r = (a.x*zi², a.y*zi³). a must not be infinity. */
BLVM_HD_INLINE void ge_set_gej_zinv(Ge* r, const Gej* a, const Fe* zi) {
    Fe zi2, zi3;
    fe_sqr(&zi2, zi);
    fe_mul(&zi3, &zi2, zi);
    fe_mul(&r->x, &a->x, &zi2);
    fe_mul(&r->y, &a->y, &zi3);
    r->infinity = 0;
}

/* r = (a.x*zi², a.y*zi³). a must not be infinity. */
BLVM_HD_INLINE void ge_set_ge_zinv(Ge* r, const Ge* a, const Fe* zi) {
    Fe zi2, zi3;
    fe_sqr(&zi2, zi);
    fe_mul(&zi3, &zi2, zi);
    fe_mul(&r->x, &a->x, &zi2);
    fe_mul(&r->y, &a->y, &zi3);
    r->infinity = 0;
}

BLVM_HD_INLINE void gej_set_inf(Gej* r) {
    fe_set0(&r->x);
    fe_set0(&r->y);
    fe_set0(&r->z);
    r->infinity = 1;
}

BLVM_HD_INLINE void gej_set_xy(Gej* r, const Fe* x, const Fe* y) {
    r->x = *x;
    r->y = *y;
    fe_set1(&r->z);
    r->infinity = 0;
}

BLVM_HD_INLINE void gej_set_ge(Gej* r, const Ge* a) {
    if (a->infinity) {
        gej_set_inf(r);
        return;
    }
    r->x = a->x;
    r->y = a->y;
    fe_set1(&r->z);
    r->infinity = 0;
}

BLVM_HD void gej_double(Gej* r, const Gej* p) {
    if (p->infinity || fe_is0(&p->y)) {
        gej_set_inf(r);
        return;
    }
    Fe A, B, C, D, E, F, X3, Y3, Z3;
    fe_sqr(&A, &p->x);
    fe_sqr(&B, &p->y);
    fe_sqr(&C, &B);
    Fe xb;
    fe_add(&xb, &p->x, &B);
    fe_sqr(&D, &xb);
    fe_sub(&D, &D, &A);
    fe_sub(&D, &D, &C);
    fe_add(&D, &D, &D);
    fe_add(&E, &A, &A);
    fe_add(&E, &E, &A);
    fe_sqr(&F, &E);
    Fe D2;
    fe_add(&D2, &D, &D);
    fe_sub(&X3, &F, &D2);
    Fe DX;
    fe_sub(&DX, &D, &X3);
    fe_mul(&Y3, &E, &DX);
    Fe C8 = C;
    fe_add(&C8, &C8, &C8);
    fe_add(&C8, &C8, &C8);
    fe_add(&C8, &C8, &C8);
    fe_sub(&Y3, &Y3, &C8);
    fe_mul(&Z3, &p->y, &p->z);
    fe_add(&Z3, &Z3, &Z3);
    r->x = X3;
    r->y = Y3;
    r->z = Z3;
    r->infinity = 0;
}

BLVM_HD void gej_add(Gej* r, const Gej* p, const Gej* q) {
    if (p->infinity) {
        *r = *q;
        return;
    }
    if (q->infinity) {
        *r = *p;
        return;
    }
    Fe Z1Z1, Z2Z2, U1, U2, S1, S2, H, HH, HHH, rfe, V, X3, Y3, Z3;
    fe_sqr(&Z1Z1, &p->z);
    fe_sqr(&Z2Z2, &q->z);
    fe_mul(&U1, &p->x, &Z2Z2);
    fe_mul(&U2, &q->x, &Z1Z1);
    Fe Z1cub, Z2cub;
    fe_mul(&Z1cub, &p->z, &Z1Z1);
    fe_mul(&Z2cub, &q->z, &Z2Z2);
    fe_mul(&S1, &p->y, &Z2cub);
    fe_mul(&S2, &q->y, &Z1cub);
    fe_sub(&H, &U2, &U1);
    fe_sub(&rfe, &S2, &S1);
    if (fe_is0(&H)) {
        if (fe_is0(&rfe)) {
            gej_double(r, p);
            return;
        }
        gej_set_inf(r);
        return;
    }
    fe_sqr(&HH, &H);
    fe_mul(&HHH, &HH, &H);
    fe_mul(&V, &U1, &HH);
    fe_sqr(&X3, &rfe);
    fe_sub(&X3, &X3, &HHH);
    Fe V2;
    fe_add(&V2, &V, &V);
    fe_sub(&X3, &X3, &V2);
    fe_sub(&Y3, &V, &X3);
    fe_mul(&Y3, &rfe, &Y3);
    Fe S1HHH;
    fe_mul(&S1HHH, &S1, &HHH);
    fe_sub(&Y3, &Y3, &S1HHH);
    /* Z3 = Z1 * Z2 * H */
    fe_mul(&Z3, &p->z, &q->z);
    fe_mul(&Z3, &Z3, &H);
    r->x = X3;
    r->y = Y3;
    r->z = Z3;
    r->infinity = 0;
}

/* Mixed add: r = a + b with b affine (host gej_add_ge_var / add_ge_var_rzr).
 * If rzr != nullptr, *rzr = Z(r)/Z(a) (= H), or 1/0 on infinity shortcuts. */
BLVM_HD void gej_add_ge_var(Gej* r, const Gej* a, const Ge* b, Fe* rzr) {
    if (a->infinity) {
        if (rzr) fe_set1(rzr);
        gej_set_ge(r, b);
        return;
    }
    if (b->infinity) {
        if (rzr) fe_set0(rzr);
        *r = *a;
        return;
    }
    Fe z12, u2, s2, h, i, h2, h3, t, tmp;
    fe_sqr(&z12, &a->z);
    fe_mul(&u2, &b->x, &z12);
    fe_mul(&tmp, &b->y, &z12);
    fe_mul(&s2, &tmp, &a->z);
    fe_sub(&h, &u2, &a->x);
    fe_sub(&i, &a->y, &s2);
    if (fe_is0(&h)) {
        if (rzr) fe_set0(rzr);
        if (fe_is0(&i)) {
            gej_double(r, a);
        } else {
            gej_set_inf(r);
        }
        return;
    }
    if (rzr) *rzr = h;
    fe_mul(&r->z, &a->z, &h);
    fe_sqr(&h2, &h);
    fe_neg(&h2, &h2);
    fe_mul(&h3, &h2, &h);
    fe_mul(&t, &a->x, &h2);
    fe_sqr(&r->x, &i);
    fe_add(&r->x, &r->x, &h3);
    fe_add(&r->x, &r->x, &t);
    fe_add(&r->x, &r->x, &t);
    fe_add(&t, &t, &r->x);
    fe_mul(&r->y, &t, &i);
    fe_mul(&tmp, &h3, &a->y);
    fe_add(&r->y, &r->y, &tmp);
    r->infinity = 0;
}

/* r = a + b with b affine and bzinv = 1/Z_b in the deferred-Z Strauss system
 * (host gej_add_zinv_var). When bzinv=1 this equals gej_add_ge_var. */
BLVM_HD void gej_add_zinv_var(Gej* r, const Gej* a, const Ge* b, const Fe* bzinv) {
    if (a->infinity) {
        if (b->infinity) {
            gej_set_inf(r);
            return;
        }
        Fe bzinv2, bzinv3;
        fe_sqr(&bzinv2, bzinv);
        fe_mul(&bzinv3, &bzinv2, bzinv);
        fe_mul(&r->x, &b->x, &bzinv2);
        fe_mul(&r->y, &b->y, &bzinv3);
        fe_set1(&r->z);
        r->infinity = 0;
        return;
    }
    if (b->infinity) {
        *r = *a;
        return;
    }
    Fe az, z12, u2, s2, h, i, h2, h3, t, tmp;
    fe_mul(&az, &a->z, bzinv);
    fe_sqr(&z12, &az);
    fe_mul(&u2, &b->x, &z12);
    fe_mul(&tmp, &b->y, &z12);
    fe_mul(&s2, &tmp, &az);
    fe_sub(&h, &u2, &a->x);
    fe_sub(&i, &a->y, &s2);
    if (fe_is0(&h)) {
        if (fe_is0(&i)) {
            gej_double(r, a);
        } else {
            gej_set_inf(r);
        }
        return;
    }
    fe_mul(&r->z, &a->z, &h);
    fe_sqr(&h2, &h);
    fe_neg(&h2, &h2);
    fe_mul(&h3, &h2, &h);
    fe_mul(&t, &a->x, &h2);
    fe_sqr(&r->x, &i);
    fe_add(&r->x, &r->x, &h3);
    fe_add(&r->x, &r->x, &t);
    fe_add(&r->x, &r->x, &t);
    fe_add(&t, &t, &r->x);
    fe_mul(&r->y, &t, &i);
    fe_mul(&tmp, &h3, &a->y);
    fe_add(&r->y, &r->y, &tmp);
    r->infinity = 0;
}

BLVM_HD int ge_from_compressed(Gej* r, const uint8_t pk33[33]) {
    if (pk33[0] != 2 && pk33[0] != 3) return 0;
    Fe x, x2, x3, y2, y, seven;
    if (!fe_set_b32_limit(&x, pk33 + 1)) return 0;
    fe_sqr(&x2, &x);
    fe_mul(&x3, &x2, &x);
    fe_set0(&seven);
    seven.d[0] = 7;
    fe_add(&y2, &x3, &seven);
    if (!fe_sqrt(&y, &y2)) return 0;
    if (((pk33[0] == 3) ? 1 : 0) != fe_is_odd(&y)) fe_neg(&y, &y);
    gej_set_xy(r, &x, &y);
    return 1;
}

BLVM_HD void gej_get_x(Fe* xout, const Gej* a) {
    if (a->infinity) {
        fe_set0(xout);
        return;
    }
    Fe zi, zi2;
    fe_inv(&zi, &a->z);
    fe_sqr(&zi2, &zi);
    fe_mul(xout, &a->x, &zi2);
    fe_normalize(xout);
}

BLVM_HD void gej_to_affine(Ge* r, const Gej* a) {
    if (a->infinity) {
        ge_set_inf(r);
        return;
    }
    if (a->z.d[0] == 1 && a->z.d[1] == 0 && a->z.d[2] == 0 && a->z.d[3] == 0) {
        r->x = a->x;
        r->y = a->y;
        r->infinity = 0;
        return;
    }
    Fe zi;
    fe_inv(&zi, &a->z);
    ge_set_gej_zinv(r, a, &zi);
}

BLVM_HD void gej_set_g(Gej* r) {
    const uint8_t gx[32] = {0x79, 0xBE, 0x66, 0x7E, 0xF9, 0xDC, 0xBB, 0xAC, 0x55, 0xA0, 0x62, 0x95,
                            0xCE, 0x87, 0x0B, 0x07, 0x02, 0x9B, 0xFC, 0xDB, 0x2D, 0xCE, 0x28, 0xD9,
                            0x59, 0xF2, 0x81, 0x5B, 0x16, 0xF8, 0x17, 0x98};
    const uint8_t gy[32] = {0x48, 0x3A, 0xDA, 0x77, 0x26, 0xA3, 0xC4, 0x65, 0x5D, 0xA4, 0xFB, 0xFC,
                            0x0E, 0x11, 0x08, 0xA8, 0xFD, 0x17, 0xB4, 0x48, 0xA6, 0x85, 0x54, 0x19,
                            0x9C, 0x47, 0xD0, 0x8F, 0xFB, 0x10, 0xD4, 0xB8};
    Fe x, y;
    fe_from_be32(&x, gx);
    fe_from_be32(&y, gy);
    gej_set_xy(r, &x, &y);
}
