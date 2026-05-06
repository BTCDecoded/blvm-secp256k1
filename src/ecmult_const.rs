//! Constant-time point multiplication: `Q = q*A` (libsecp256k1 `ecmult_const` algorithm).
//! Port of [ecmult_const_impl.h](https://github.com/bitcoin-core/secp256k1) (5-bit group, GLV, effective affine).

use std::sync::OnceLock;
use subtle::{Choice, ConstantTimeEq};

use crate::ecmult::ecmult_odd_multiples_table;
use crate::field::FieldElement;
use crate::group::{ge_table_set_globalz, generator_g, Ge, Gej};
use crate::scalar::Scalar;

const ECMULT_CONST_GROUP_SIZE: u32 = 5;
const ECMULT_CONST_TABLE_SIZE: usize = 1 << (ECMULT_CONST_GROUP_SIZE as usize - 1);
const ECMULT_CONST_GROUPS: i32 = 26;
// (129 + 4) / 5 = 26 groups; 26 * 5 = 130 = ECMULT_CONST_BITS
const _ECMULT_CONST_BITS: u32 = ECMULT_CONST_GROUPS as u32 * ECMULT_CONST_GROUP_SIZE;

// K for ECMULT_CONST_BITS=130, GROUP_SIZE=5. Packed like scalar_4x64.h
// `SECP256K1_SCALAR_CONST(d7,…,d0)` => (d1<<32|d0, d3<<32|d2, d5<<32|d4, d7<<32|d6).
const ECMULT_K: Scalar = Scalar {
    d: [
        0xb5c2c1dcde9798d9u64,
        0x589ae84826ba29e4u64,
        0xc2bdd6bf7c118d6bu64,
        0xa4e88a7dcb13034eu64,
    ],
};

// 2^128 (S_OFFSET in upstream) as 4x64 little-endian scalar limbs.
const S_OFFSET: Scalar = Scalar {
    d: [0, 0, 1, 0],
};

#[allow(clippy::needless_range_loop)] // index drives constant-time cmove selection
fn ecmult_const_table_get_ge(r: &mut Ge, pre: &[Ge; ECMULT_CONST_TABLE_SIZE], n: u32) {
    let neg = ((n >> (ECMULT_CONST_GROUP_SIZE - 1)) & 1) ^ 1;
    let negm = 0u32.wrapping_sub(neg);
    let index = (negm ^ n) & ((1u32 << (ECMULT_CONST_GROUP_SIZE - 1)) - 1);
    *r = pre[0];
    for m in 1..ECMULT_CONST_TABLE_SIZE {
        let eq = (m as u32).ct_eq(&index);
        r.x.cmov(&pre[m].x, eq);
        r.y.cmov(&pre[m].y, eq);
    }
    r.infinity = false;
    let mut neg_y = FieldElement::zero();
    neg_y.negate(&r.y, 1);
    r.y.cmov(&neg_y, Choice::from(neg as u8));
}

/// Decompose `q` into the two GLV halves `(v1, v2)` plus the `S_OFFSET` shift used by the
/// scalar-walk loop. Same algebra as the upstream prologue.
#[inline(always)]
fn ecmult_const_decompose(q: &Scalar) -> (Scalar, Scalar) {
    let mut s = Scalar::zero();
    s.add(q, &ECMULT_K);
    let ssum = s;
    s.half_modn(&ssum);

    let mut v1 = Scalar::zero();
    let mut v2 = Scalar::zero();
    Scalar::split_lambda(&mut v1, &mut v2, &s);

    let v1o = v1;
    v1.add(&v1o, &S_OFFSET);
    let v2o = v2;
    v2.add(&v2o, &S_OFFSET);
    (v1, v2)
}

/// Constant-time scalar walk over already-prepared **affine** tables (`pre_a` and `pre_a_lam`)
/// with a shared `global_z`. Public-data layout of the tables; secret-dependent indexing is
/// resolved via cmov-only reads (`ecmult_const_table_get_ge`).
#[inline(always)]
fn ecmult_const_walk(
    r: &mut Gej,
    pre_a: &[Ge; ECMULT_CONST_TABLE_SIZE],
    pre_a_lam: &[Ge; ECMULT_CONST_TABLE_SIZE],
    global_z: &FieldElement,
    v1: &Scalar,
    v2: &Scalar,
) {
    let mut first = true;
    for g in (0..ECMULT_CONST_GROUPS).rev() {
        let off = (g as u32) * ECMULT_CONST_GROUP_SIZE;
        let bits1 = v1.get_bits_var(off, ECMULT_CONST_GROUP_SIZE);
        let bits2 = v2.get_bits_var(off, ECMULT_CONST_GROUP_SIZE);
        let mut t = Ge::default();
        ecmult_const_table_get_ge(&mut t, pre_a, bits1);
        if first {
            r.set_ge(&t);
            first = false;
        } else {
            for _ in 0..ECMULT_CONST_GROUP_SIZE {
                let rj = *r;
                r.double_ct(&rj);
            }
            let rj = *r;
            r.add_ge(&rj, &t);
        }
        ecmult_const_table_get_ge(&mut t, pre_a_lam, bits2);
        let rj = *r;
        r.add_ge(&rj, &t);
    }
    let r_z = r.z;
    r.z.mul(&r_z, global_z);
}

/// R = `q*A` (constant-time w.r.t. `q`); public point `A` is allowed to use variable precomputation.
/// If `A` is infinity, `r` is set to infinity.
pub fn ecmult_const(r: &mut Gej, a: &Ge, q: &Scalar) {
    if a.infinity {
        r.set_infinity();
        return;
    }

    let (v1, v2) = ecmult_const_decompose(q);

    let mut tj = Gej::default();
    tj.set_ge(a);
    let mut pre_a = [Ge::default(); ECMULT_CONST_TABLE_SIZE];
    let mut zr = [FieldElement::zero(); ECMULT_CONST_TABLE_SIZE];
    let mut global_z = FieldElement::one();
    ecmult_odd_multiples_table(
        ECMULT_CONST_TABLE_SIZE,
        &mut pre_a,
        &mut zr,
        &mut global_z,
        &tj,
    );
    ge_table_set_globalz(ECMULT_CONST_TABLE_SIZE, &mut pre_a, &zr);
    let mut pre_a_lam = [Ge::default(); ECMULT_CONST_TABLE_SIZE];
    for i in 0..ECMULT_CONST_TABLE_SIZE {
        pre_a_lam[i].mul_lambda(&pre_a[i]);
    }

    ecmult_const_walk(r, &pre_a, &pre_a_lam, &global_z, &v1, &v2);
}

/// Cached `pre_a` / `pre_a_lam` tables and `global_z` for the **generator G**, built once.
/// All entries are public data; only the *index* into them (per scalar limb group) is secret.
struct GenConstTables {
    pre_a: [Ge; ECMULT_CONST_TABLE_SIZE],
    pre_a_lam: [Ge; ECMULT_CONST_TABLE_SIZE],
    global_z: FieldElement,
}

fn gen_const_tables() -> &'static GenConstTables {
    static T: OnceLock<GenConstTables> = OnceLock::new();
    T.get_or_init(|| {
        let g = generator_g();
        let mut tj = Gej::default();
        tj.set_ge(&g);

        let mut pre_a = [Ge::default(); ECMULT_CONST_TABLE_SIZE];
        let mut zr = [FieldElement::zero(); ECMULT_CONST_TABLE_SIZE];
        let mut global_z = FieldElement::one();
        ecmult_odd_multiples_table(
            ECMULT_CONST_TABLE_SIZE,
            &mut pre_a,
            &mut zr,
            &mut global_z,
            &tj,
        );
        ge_table_set_globalz(ECMULT_CONST_TABLE_SIZE, &mut pre_a, &zr);

        let mut pre_a_lam = [Ge::default(); ECMULT_CONST_TABLE_SIZE];
        for i in 0..ECMULT_CONST_TABLE_SIZE {
            pre_a_lam[i].mul_lambda(&pre_a[i]);
        }

        for entry in pre_a.iter_mut() {
            entry.x.normalize();
            entry.y.normalize();
        }
        for entry in pre_a_lam.iter_mut() {
            entry.x.normalize();
            entry.y.normalize();
        }
        global_z.normalize();

        GenConstTables {
            pre_a,
            pre_a_lam,
            global_z,
        }
    })
}

/// Legacy GLV-based constant-time `k*G`. Retained for A/B comparison against the
/// multi-comb path (`ecmult_gen_comb`). The production signing path uses the comb;
/// this routine remains crate-visible so tests (or an external harness) can compare both.
#[allow(dead_code)]
pub(crate) fn ecmult_gen_const_glv(r: &mut Gej, k: &Scalar) {
    let (v1, v2) = ecmult_const_decompose(k);
    let t = gen_const_tables();
    ecmult_const_walk(r, &t.pre_a, &t.pre_a_lam, &t.global_z, &v1, &v2);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecmult;
    use sha2::{Digest, Sha256};

    fn scalar_small(v: u64) -> Scalar {
        let mut s = Scalar::zero();
        s.set_int(v as u32);
        s
    }

    /// Mix add 2*G+G=3*G in Jacobian must match ecmult_gen(3).
    #[test]
    fn add_ge_2g_plus_g() {
        let g = generator_g();
        let two = scalar_small(2);
        let three = scalar_small(3);
        let mut t2g = Gej::default();
        ecmult::ecmult_gen(&mut t2g, &two);
        let mut r = Gej::default();
        r.add_ge(&t2g, &g);
        let mut e3 = Gej::default();
        ecmult::ecmult_gen(&mut e3, &three);
        let mut got = Ge::default();
        got.set_gej_var(&r);
        let mut want = Ge::default();
        want.set_gej_var(&e3);
        got.x.normalize();
        want.x.normalize();
        assert_eq!(got.x, want.x, "add_ge(2G,G) x");
    }

    #[test]
    fn ecmult_const_vs_ecmult_small_int() {
        let g = generator_g();
        let mut gj = Gej::default();
        gj.set_ge(&g);
        for v in 1u32..50u32 {
            let mut s = Scalar::zero();
            s.set_int(v);
            let mut a = Gej::default();
            ecmult::ecmult(&mut a, &gj, &s, None);
            let mut c = Gej::default();
            ecmult_const(&mut c, &g, &s);
            let mut ag = Ge::default();
            ag.set_gej_var(&a);
            let mut cg = Ge::default();
            cg.set_gej_var(&c);
            ag.x.normalize();
            cg.x.normalize();
            ag.y.normalize();
            cg.y.normalize();
            assert_eq!(ag.x, cg.x, "ecmult vs ecmult_const x for k={v}");
            assert_eq!(ag.y, cg.y, "ecmult vs ecmult_const y for k={v}");
        }
    }

    /// Deterministic "random" scalars: `ecmult` vs `ecmult_const` on generator must agree.
    #[test]
    fn ecmult_const_vs_ecmult_hashed() {
        let g = generator_g();
        let mut gj = Gej::default();
        gj.set_ge(&g);
        for i in 0u32..256 {
            let h = Sha256::digest((i as u64).to_le_bytes());
            let b32: [u8; 32] = h.into();
            let mut s = Scalar::zero();
            let _ = s.set_b32(&b32);
            if s.is_zero() {
                continue;
            }
            let mut a = Gej::default();
            ecmult::ecmult(&mut a, &gj, &s, None);
            let mut c = Gej::default();
            ecmult_const(&mut c, &g, &s);
            let mut ag = Ge::default();
            ag.set_gej_var(&a);
            let mut cg = Ge::default();
            cg.set_gej_var(&c);
            ag.x.normalize();
            cg.x.normalize();
            ag.y.normalize();
            cg.y.normalize();
            assert_eq!(ag.x, cg.x, "x i={i}");
            assert_eq!(ag.y, cg.y, "y i={i}");
        }
    }
}
