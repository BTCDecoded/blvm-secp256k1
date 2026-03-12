//! Field arithmetic tests and cross-validation.

use blvm_secp256k1::field::FieldElement;

#[test]
fn test_field_zero_one() {
    let z = FieldElement::zero();
    assert!(z.is_zero());

    let o = FieldElement::one();
    assert!(!o.is_zero());

    let mut buf = [0u8; 32];
    o.get_b32(&mut buf);
    assert_eq!(buf[31], 1);
    assert!(buf[..31].iter().all(|&b| b == 0));
}

#[test]
fn test_field_set_get_roundtrip() {
    // Use a value < p: 0x1234... (32 bytes, big-endian)
    let bytes =
        hex::decode("0000000000000000000000000000000000000000000000000000000000001234").unwrap();
    let mut arr = [0u8; 32];
    arr.copy_from_slice(&bytes);

    let mut fe = FieldElement::zero();
    let ok = fe.set_b32_limit(&arr);
    assert!(ok, "set_b32_limit should accept this value");
    fe.normalize();
    let mut out = [0u8; 32];
    fe.get_b32(&mut out);
    assert_eq!(arr, out, "roundtrip failed");
}

#[test]
fn test_field_add() {
    let mut a = FieldElement::one();
    let b = FieldElement::one();
    a.add_assign(&b);
    a.normalize();

    let mut expected = FieldElement::zero();
    expected.set_b32_mod(&[0u8; 32]);
    expected.n[0] = 2;
    expected.normalize();

    assert_eq!(a.n[0], 2);
    assert!(a.n[1..].iter().all(|&x| x == 0));
}

#[test]
fn test_field_mul_sqr() {
    // a = 2, b = 3 => a*b = 6
    let mut a = FieldElement::zero();
    a.n[0] = 2;
    let mut b = FieldElement::zero();
    b.n[0] = 3;

    let mut ab = FieldElement::zero();
    ab.mul(&a, &b);
    ab.normalize();
    assert_eq!(ab.n[0], 6);
    assert!(ab.n[1..].iter().all(|&x| x == 0));

    // a^2 = 4
    let mut a2 = FieldElement::zero();
    a2.sqr(&a);
    a2.normalize();
    assert_eq!(a2.n[0], 4);
}

#[test]
fn test_field_inv() {
    // inv(2) * 2 = 1 mod p.
    let mut a = FieldElement::zero();
    a.n[0] = 2;

    let mut inv_a = FieldElement::zero();
    inv_a.inv(&a);

    let mut prod = FieldElement::zero();
    prod.mul(&a, &inv_a);
    prod.normalize();

    let one = FieldElement::one();
    assert_eq!(prod, one, "a * inv(a) should equal 1");
}

#[test]
fn test_field_inv_one() {
    // inv(1) = 1
    let a = FieldElement::one();
    let mut inv_a = FieldElement::zero();
    inv_a.inv(&a);
    assert_eq!(inv_a, a);
}

#[test]
fn test_field_negate() {
    // Test that a + (-a) = 0 mod p. Use small a = 5.
    let mut a = FieldElement::zero();
    a.n[0] = 5;

    let mut neg = FieldElement::zero();
    neg.negate(&a, 1);

    let mut sum = FieldElement::zero();
    sum.add_assign(&a);
    sum.add_assign(&neg);
    sum.normalize();

    assert!(sum.is_zero(), "a + (-a) should be 0");
}

#[test]
fn test_field_normalizes_to_zero_var() {
    // -a + a should normalize to 0 (used in gej_add when P1 = P2)
    // Use secp256k1 generator G x-coordinate (well-known constant)
    // secp256k1 generator G x-coordinate (32 bytes = 64 hex chars)
    let gx =
        hex::decode("79be667ef9dcbbac55a06295ce870b0729bfcdb2dce28d959f2815b16f817998").unwrap();
    let mut gx_arr = [0u8; 32];
    gx_arr.copy_from_slice(&gx);
    let mut a = FieldElement::zero();
    a.set_b32_mod(&gx_arr);

    let mut h = FieldElement::zero();
    h.negate(&a, 1);
    h.add_assign(&a);
    assert!(
        h.normalizes_to_zero_var(),
        "-a + a must normalize to zero (P+P doubling path)"
    );
}
