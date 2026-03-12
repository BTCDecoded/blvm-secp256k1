//! Group operation tests.

use blvm_secp256k1::field::FieldElement;
use blvm_secp256k1::group::{ge_set_all_gej_var, generator_g, Ge, Gej};

/// Generator G matches secp256k1 curve (SEC2 2.7.1)
#[test]
fn test_generator_matches_curve() {
    let g = generator_g();
    let mut gx = g.x;
    let mut gy = g.y;
    gx.normalize();
    gy.normalize();
    let mut our_x = [0u8; 32];
    let mut our_y = [0u8; 32];
    gx.get_b32(&mut our_x);
    gy.get_b32(&mut our_y);
    // SEC2 secp256k1 generator G (32 bytes = 64 hex chars each)
    let expected_x =
        hex::decode("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798").unwrap();
    let expected_y =
        hex::decode("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8").unwrap();
    assert_eq!(our_x, expected_x.as_slice(), "G.x mismatch");
    assert_eq!(our_y, expected_y.as_slice(), "G.y mismatch");
}

#[test]
fn test_ge_generator_on_curve() {
    let g = generator_g();
    assert!(!g.is_infinity());

    // Verify y^2 = x^3 + 7 (secp256k1 curve equation)
    let mut y2 = FieldElement::zero();
    y2.sqr(&g.y);
    y2.normalize();
    let mut x2 = FieldElement::zero();
    x2.sqr(&g.x);
    let mut x3 = FieldElement::zero();
    x3.mul(&x2, &g.x);
    x3.add_int(7);
    x3.normalize();
    assert_eq!(y2, x3, "G must satisfy y^2 = x^3 + 7");
}

#[test]
fn test_gej_set_ge_roundtrip() {
    let g = generator_g();
    let mut gej = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    gej.set_ge(&g);
    assert!(!gej.is_infinity());

    let mut ge2 = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    ge2.set_gej_var(&gej);
    assert!(!ge2.is_infinity());
    ge2.x.normalize();
    ge2.y.normalize();
    let mut gx = g.x;
    let mut gy = g.y;
    gx.normalize();
    gy.normalize();
    assert_eq!(ge2.x, gx);
    assert_eq!(ge2.y, gy);
}

#[test]
fn test_gej_double() {
    let g = generator_g();
    let mut gej = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    gej.set_ge(&g);

    let mut dbl = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    dbl.double(&gej);
    assert!(!dbl.is_infinity());

    let mut dbl_aff = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    dbl_aff.set_gej_var(&dbl);
    assert!(!dbl_aff.is_infinity());

    // 2G should be on curve
    let mut y2 = FieldElement::zero();
    y2.sqr(&dbl_aff.y);
    let mut x2 = FieldElement::zero();
    x2.sqr(&dbl_aff.x);
    let mut x3 = FieldElement::zero();
    x3.mul(&x2, &dbl_aff.x);
    x3.add_int(7);
    y2.normalize();
    x3.normalize();
    assert_eq!(y2, x3, "2G must be on curve");
}

/// G + G = 2G via add_var (same-point doubling path)
#[test]
fn test_gej_add_var_same_point() {
    let g = generator_g();
    let mut gej = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    gej.set_ge(&g);

    let mut sum = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    sum.add_var(&gej, &gej);
    assert!(!sum.is_infinity(), "G + G should not be infinity");

    let mut dbl = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    dbl.double(&gej);

    let mut sum_aff = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    sum_aff.set_gej_var(&sum);
    let mut dbl_aff = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    dbl_aff.set_gej_var(&dbl);
    sum_aff.x.normalize();
    sum_aff.y.normalize();
    dbl_aff.x.normalize();
    dbl_aff.y.normalize();
    assert_eq!(sum_aff.x, dbl_aff.x, "G+G must equal 2G");
    assert_eq!(sum_aff.y, dbl_aff.y);
}

/// G + G = 2G via add_ge_var (same-point doubling path)
#[test]
fn test_gej_add_ge_var_same_point() {
    let g = generator_g();
    let mut gej = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    gej.set_ge(&g);

    let mut sum = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    sum.add_ge_var(&gej, &g);
    assert!(
        !sum.is_infinity(),
        "G + G (add_ge_var) should not be infinity"
    );

    let mut dbl = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    dbl.double(&gej);

    let mut sum_aff = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    sum_aff.set_gej_var(&sum);
    let mut dbl_aff = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    dbl_aff.set_gej_var(&dbl);
    sum_aff.x.normalize();
    sum_aff.y.normalize();
    dbl_aff.x.normalize();
    dbl_aff.y.normalize();
    assert_eq!(sum_aff.x, dbl_aff.x, "G+G (add_ge_var) must equal 2G");
    assert_eq!(sum_aff.y, dbl_aff.y);
}

/// G + 2G = 3G via add (different points)
#[test]
fn test_gej_add_var() {
    let g = generator_g();
    let mut gej_g = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    gej_g.set_ge(&g);
    let mut gej_2g = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    gej_2g.double(&gej_g);

    let mut sum = Gej {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        z: FieldElement::zero(),
        infinity: false,
    };
    sum.add_var(&gej_g, &gej_2g);
    assert!(!sum.is_infinity(), "G + 2G should not be infinity");

    let mut sum_aff = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    sum_aff.set_gej_var(&sum);
    assert!(!sum_aff.is_infinity());

    // Verify 3G is on curve
    let mut y2 = FieldElement::zero();
    y2.sqr(&sum_aff.y);
    y2.normalize();
    let mut x2 = FieldElement::zero();
    x2.sqr(&sum_aff.x);
    let mut x3 = FieldElement::zero();
    x3.mul(&x2, &sum_aff.x);
    x3.add_int(7);
    x3.normalize();
    assert_eq!(y2, x3, "3G must be on curve");
}

#[test]
fn test_ge_neg() {
    let g = generator_g();
    let mut neg = Ge {
        x: FieldElement::zero(),
        y: FieldElement::zero(),
        infinity: false,
    };
    neg.neg(&g);
    assert!(!neg.is_infinity());
    neg.x.normalize();
    neg.y.normalize();
    let mut gx = g.x;
    gx.normalize();
    assert_eq!(neg.x, gx);
    let mut expected_neg_y = FieldElement::zero();
    expected_neg_y.negate(&g.y, 1);
    expected_neg_y.normalize();
    assert_eq!(neg.y, expected_neg_y);
}

/// ge_set_all_gej_var: batch Gej->Ge matches individual set_gej_var
#[test]
fn test_ge_set_all_gej_var() {
    let g = generator_g();
    let mut gej = Gej::default();
    gej.set_ge(&g);

    // Build [G, 2G, 3G] in Jacobian
    let mut gejs = vec![Gej::default(); 3];
    gejs[0] = gej;
    let mut sum = Gej::default();
    sum.add_ge_var(&gej, &g);
    gejs[1] = sum;
    let sum_in = sum;
    sum.add_ge_var(&sum_in, &g);
    gejs[2] = sum;

    // Batch convert
    let mut batch_aff = vec![Ge::default(); 3];
    ge_set_all_gej_var(&mut batch_aff, &gejs);

    // Individual convert and compare
    for i in 0..3 {
        let mut ind_aff = Ge::default();
        ind_aff.set_gej_var(&gejs[i]);
        batch_aff[i].x.normalize();
        batch_aff[i].y.normalize();
        ind_aff.x.normalize();
        ind_aff.y.normalize();
        assert_eq!(batch_aff[i].x, ind_aff.x, "batch[{}].x", i);
        assert_eq!(batch_aff[i].y, ind_aff.y, "batch[{}].y", i);
    }
}
