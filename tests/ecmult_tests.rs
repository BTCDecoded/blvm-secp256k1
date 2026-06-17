//! Ecmult tests: scalar multiplication R = k*G and R = na*A + ng*G.

use blvm_secp256k1::ecmult::{ecmult, ecmult_gen};
use blvm_secp256k1::group::{Ge, Gej, generator_g};
use blvm_secp256k1::scalar::Scalar;

/// ecmult_gen(1) = G
#[test]
fn test_ecmult_gen_one() {
    let g = generator_g();
    let mut one = Scalar::zero();
    one.set_int(1);

    let mut r = Gej::default();
    ecmult_gen(&mut r, &one);

    assert!(!r.is_infinity());
    let mut aff = Ge::default();
    aff.set_gej_var(&r);
    aff.x.normalize();
    aff.y.normalize();
    let mut gx = g.x;
    let mut gy = g.y;
    gx.normalize();
    gy.normalize();
    assert_eq!(aff.x, gx, "ecmult_gen(1) must equal G");
    assert_eq!(aff.y, gy, "ecmult_gen(1) must equal G");
}

/// ecmult_gen(2) = 2*G (cross-validate with G+G)
#[test]
fn test_ecmult_gen_two() {
    let g = generator_g();
    let mut gej = Gej::default();
    gej.set_ge(&g);
    let mut sum = Gej::default();
    sum.add_ge_var(&gej, &g);
    let mut two_g_aff = Ge::default();
    two_g_aff.set_gej_var(&sum);
    two_g_aff.x.normalize();
    two_g_aff.y.normalize();

    let mut two = Scalar::zero();
    two.set_int(2);
    let mut r = Gej::default();
    ecmult_gen(&mut r, &two);

    assert!(!r.is_infinity());
    let mut aff = Ge::default();
    aff.set_gej_var(&r);
    aff.x.normalize();
    aff.y.normalize();
    assert_eq!(aff.x, two_g_aff.x, "ecmult_gen(2) must equal 2G");
    assert_eq!(aff.y, two_g_aff.y, "ecmult_gen(2) must equal 2G");
}

/// ecmult_gen(0) = infinity
#[test]
fn test_ecmult_gen_zero() {
    let zero = Scalar::zero();
    let mut r = Gej::default();
    ecmult_gen(&mut r, &zero);
    assert!(r.is_infinity(), "ecmult_gen(0) must be infinity");
}

/// ecmult with na=0, ng=k gives k*G (same as ecmult_gen)
#[test]
fn test_ecmult_gen_via_ecmult() {
    let mut k = Scalar::zero();
    k.set_int(42);

    let mut r_gen = Gej::default();
    ecmult_gen(&mut r_gen, &k);

    let mut inf = Gej::default();
    inf.set_infinity();
    let zero = Scalar::zero();
    let mut r_ecmult = Gej::default();
    ecmult(&mut r_ecmult, &inf, &zero, Some(&k));

    let mut aff_gen = Ge::default();
    aff_gen.set_gej_var(&r_gen);
    let mut aff_ecmult = Ge::default();
    aff_ecmult.set_gej_var(&r_ecmult);
    aff_gen.x.normalize();
    aff_gen.y.normalize();
    aff_ecmult.x.normalize();
    aff_ecmult.y.normalize();
    assert_eq!(
        aff_gen.x, aff_ecmult.x,
        "ecmult(inf,0,k) must equal ecmult_gen(k)"
    );
    assert_eq!(
        aff_gen.y, aff_ecmult.y,
        "ecmult(inf,0,k) must equal ecmult_gen(k)"
    );
}
