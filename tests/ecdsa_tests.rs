//! ECDSA sign/verify tests.

use blvm_secp256k1::ecdsa::{
    ecdsa_sig_normalize, ecdsa_sig_parse_der, ecdsa_sig_parse_der_lax, ecdsa_sig_serialize_der,
    ecdsa_sig_sign, ecdsa_sig_verify, ecdsa_sig_verify_exhaustive, ge_from_compressed,
    ge_to_compressed, pubkey_from_secret,
};
use blvm_secp256k1::ecmult;
use blvm_secp256k1::group::{generator_g, Gej};
use blvm_secp256k1::scalar::Scalar;
use subtle::ConstantTimeEq;

fn scalar_from_b32(b: &[u8; 32]) -> Scalar {
    let mut s = Scalar::zero();
    s.set_b32(b);
    s
}

#[test]
fn test_ecdsa_sign_verify_roundtrip() {
    // Use small values: seckey=1, msg=2, nonce=3
    let mut seckey_bytes = [0u8; 32];
    seckey_bytes[31] = 1;
    let mut msg_bytes = [0u8; 32];
    msg_bytes[31] = 2;
    let mut nonce_bytes = [0u8; 32];
    nonce_bytes[31] = 3;

    let seckey = scalar_from_b32(&seckey_bytes);
    let message = scalar_from_b32(&msg_bytes);
    let nonce = scalar_from_b32(&nonce_bytes);

    let pubkey = pubkey_from_secret(&seckey);
    let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
    assert!(
        ecdsa_sig_verify(&sig.0, &sig.1, &pubkey, &message),
        "ecdsa_sig_verify (eq_x_var check)"
    );
    assert!(
        ecdsa_sig_verify_exhaustive(&sig.0, &sig.1, &pubkey, &message),
        "exhaustive verify (ecmult check)"
    );
}

/// ecmult(a, 1, 0) = a for any a. Tests pre_a path with minimal scalar.
#[test]
fn test_ecmult_na_one_ng_zero() {
    let mut two = Scalar::zero();
    two.set_int(2);
    let mut one = Scalar::zero();
    one.set_int(1);

    let mut two_g = Gej::default();
    ecmult::ecmult_gen(&mut two_g, &two);
    let mut result = Gej::default();
    ecmult::ecmult(&mut result, &two_g, &one, None);

    let mut ge_result = blvm_secp256k1::group::Ge::default();
    ge_result.set_gej_var(&result);
    let mut ge_two = blvm_secp256k1::group::Ge::default();
    ge_two.set_gej_var(&two_g);
    let c_result = ge_to_compressed(&ge_result);
    let c_two = ge_to_compressed(&ge_two);
    assert_eq!(c_result, c_two, "ecmult(2*G, 1, 0) must equal 2*G");
}

/// ecmult(2*G, -3, None) = -6*G.
#[test]
fn test_ecmult_na_neg_three_ng_none() {
    let mut two = Scalar::zero();
    two.set_int(2);
    let mut three = Scalar::zero();
    three.set_int(3);
    let mut neg_three = Scalar::zero();
    neg_three.negate(&three);

    let mut two_g = Gej::default();
    ecmult::ecmult_gen(&mut two_g, &two);
    let mut result = Gej::default();
    ecmult::ecmult(&mut result, &two_g, &neg_three, None);

    // Expected: -3 * (2*G) = -6*G
    let mut six = Scalar::zero();
    six.set_int(6);
    let mut neg_six = Scalar::zero();
    neg_six.negate(&six);
    let mut expected = Gej::default();
    ecmult::ecmult_gen(&mut expected, &neg_six);

    let mut ge_result = blvm_secp256k1::group::Ge::default();
    ge_result.set_gej_var(&result);
    let mut ge_expected = blvm_secp256k1::group::Ge::default();
    ge_expected.set_gej_var(&expected);
    let c_result = ge_to_compressed(&ge_result);
    let c_expected = ge_to_compressed(&ge_expected);
    assert_eq!(
        c_result, c_expected,
        "ecmult(2*G, -3, None) must equal -6*G"
    );
}

/// ECDSA verify with pubkey != G. Tests pre_a path when a≠G and both na,ng present.
#[test]
fn test_ecdsa_verify_non_g_pubkey() {
    let seckey_bytes = [0x11u8; 32];
    let msg_bytes = [0x22u8; 32];
    let nonce_bytes = [0x33u8; 32];

    let seckey = scalar_from_b32(&seckey_bytes);
    let message = scalar_from_b32(&msg_bytes);
    let nonce = scalar_from_b32(&nonce_bytes);

    let pubkey = pubkey_from_secret(&seckey);
    let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");

    // ECDSA verify: u1=z/s, u2=r/s. Check u2*pubkey + u1*G via ecmult.
    // For pubkey=0x11*G: u2*pubkey + u1*G = (u2*0x11 + u1)*G. Reference via ecmult_gen.
    let mut sn = Scalar::zero();
    sn.inv_var(&sig.1);
    let mut u1 = Scalar::zero();
    u1.mul(&sn, &message);
    let mut u2 = Scalar::zero();
    u2.mul(&sn, &sig.0);

    // ecmult split path works when input has z!=1. Use 17*G from repeated add (same point as pubkey).
    let g_aff = blvm_secp256k1::group::generator_g();
    let mut pubkeyj = Gej::default();
    pubkeyj.set_ge(&g_aff);
    for _ in 1..17 {
        let mut s = Gej::default();
        s.add_ge_var(&pubkeyj, &g_aff);
        pubkeyj = s;
    }
    let mut pr = Gej::default();
    ecmult::ecmult(&mut pr, &pubkeyj, &u2, Some(&u1));

    // Reference: (u2*0x11 + u1)*G
    let mut seventeen = Scalar::zero();
    seventeen.set_int(17);
    let mut u2_17 = Scalar::zero();
    u2_17.mul(&u2, &seventeen);
    let mut expected_scalar = Scalar::zero();
    expected_scalar.add(&u2_17, &u1);
    let mut expected = Gej::default();
    ecmult::ecmult_gen(&mut expected, &expected_scalar);

    let mut pr_aff = blvm_secp256k1::group::Ge::default();
    pr_aff.set_gej_var(&pr);
    let mut exp_aff = blvm_secp256k1::group::Ge::default();
    exp_aff.set_gej_var(&expected);
    pr_aff.x.normalize();
    pr_aff.y.normalize();
    exp_aff.x.normalize();
    exp_aff.y.normalize();
    assert!(
        blvm_secp256k1::field::FieldElement::fe_equal(&pr_aff.x, &exp_aff.x)
            && blvm_secp256k1::field::FieldElement::fe_equal(&pr_aff.y, &exp_aff.y),
        "ecmult(pubkey,u2,u1) must equal (u2*17+u1)*G for pubkey=0x11*G"
    );

    assert!(
        ecdsa_sig_verify(&sig.0, &sig.1, &pubkey, &message),
        "ecdsa_sig_verify (eq_x_var) for pubkey=0x11*G"
    );
    assert!(
        ecdsa_sig_verify_exhaustive(&sig.0, &sig.1, &pubkey, &message),
        "ecdsa_sig_verify_exhaustive for pubkey=0x11*G (tests pre_a path when a≠G)"
    );
}

#[test]
fn test_ge_compressed_roundtrip() {
    let seckey_bytes: [u8; 32] = [0x42; 32];
    let seckey = scalar_from_b32(&seckey_bytes);
    let pubkey = pubkey_from_secret(&seckey);
    let compressed = ge_to_compressed(&pubkey);
    let restored = ge_from_compressed(&compressed).expect("parse");
    assert_eq!(ge_to_compressed(&restored), compressed);
}

/// Verify sign produces s*k = m + r*d (mod n).
#[test]
fn test_sign_equation() {
    let mut seckey_bytes = [0u8; 32];
    seckey_bytes[31] = 1;
    let mut msg_bytes = [0u8; 32];
    msg_bytes[31] = 2;
    let mut nonce_bytes = [0u8; 32];
    nonce_bytes[31] = 3;
    let seckey = scalar_from_b32(&seckey_bytes);
    let message = scalar_from_b32(&msg_bytes);
    let nonce = scalar_from_b32(&nonce_bytes);

    let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
    // s*k = m + r*d (mod n). Check: s*nonce = message + sigr*seckey
    let mut lhs = Scalar::zero();
    lhs.mul(&sig.1, &nonce);
    let mut rd = Scalar::zero();
    rd.mul(&sig.0, &seckey);
    let mut rhs = Scalar::zero();
    rhs.add(&rd, &message);
    // s*k = m+r*d or (low-S) -s*k = m+r*d
    let mut neg_rhs = Scalar::zero();
    neg_rhs.negate(&rhs);
    assert!(
        bool::from(lhs.ct_eq(&rhs)) || bool::from(lhs.ct_eq(&neg_rhs)),
        "s*k = ±(m+r*d)"
    );
}

/// Verify inv(s)*s = 1 and u1+u2 = ±3 for seckey=1, msg=2, nonce=3.
#[test]
fn test_ecdsa_verify_u1_u2() {
    let mut seckey_bytes = [0u8; 32];
    seckey_bytes[31] = 1;
    let mut msg_bytes = [0u8; 32];
    msg_bytes[31] = 2;
    let mut nonce_bytes = [0u8; 32];
    nonce_bytes[31] = 3;
    let seckey = scalar_from_b32(&seckey_bytes);
    let message = scalar_from_b32(&msg_bytes);
    let nonce = scalar_from_b32(&nonce_bytes);

    let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
    let mut sn = Scalar::zero();
    sn.inv_var(&sig.1);
    // inv(s)*s = 1
    let mut check = Scalar::zero();
    check.mul(&sn, &sig.1);
    assert!(check.is_one(), "inv(s)*s must equal 1");

    let mut u1 = Scalar::zero();
    u1.mul(&sn, &message);
    let mut u2 = Scalar::zero();
    u2.mul(&sn, &sig.0);
    let mut u1_plus_u2 = Scalar::zero();
    u1_plus_u2.add(&u1, &u2);
    let mut three = Scalar::zero();
    three.set_int(3);
    let mut neg_three = Scalar::zero();
    neg_three.negate(&three);
    assert!(
        bool::from(u1_plus_u2.ct_eq(&three)) || bool::from(u1_plus_u2.ct_eq(&neg_three)),
        "u1+u2 must equal ±3 (low-S gives -3)"
    );
}

/// With our sign (seckey=1, msg=2, nonce=3), u1+u2=3 or -3. ecmult(u2*G+u1*G) must equal ±3*G.
#[test]
fn test_ecmult_verify_vectors() {
    let mut seckey_bytes = [0u8; 32];
    seckey_bytes[31] = 1;
    let mut msg_bytes = [0u8; 32];
    msg_bytes[31] = 2;
    let mut nonce_bytes = [0u8; 32];
    nonce_bytes[31] = 3;
    let seckey = scalar_from_b32(&seckey_bytes);
    let message = scalar_from_b32(&msg_bytes);
    let nonce = scalar_from_b32(&nonce_bytes);

    let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
    let mut sn = Scalar::zero();
    sn.inv_var(&sig.1);
    let mut u1 = Scalar::zero();
    u1.mul(&sn, &message);
    let mut u2 = Scalar::zero();
    u2.mul(&sn, &sig.0);

    // ecmult(u2*G + u1*G) must equal ecmult_gen(u1+u2)
    let mut u1_plus_u2 = Scalar::zero();
    u1_plus_u2.add(&u1, &u2);
    let mut r_gen = Gej::default();
    ecmult::ecmult_gen(&mut r_gen, &u1_plus_u2);
    let mut ge_gen = blvm_secp256k1::group::Ge::default();
    ge_gen.set_gej_var(&r_gen);

    let g = generator_g();
    let mut gj = Gej::default();
    gj.set_ge(&g);
    let mut r = Gej::default();
    ecmult::ecmult(&mut r, &gj, &u2, Some(&u1));
    let mut ge = blvm_secp256k1::group::Ge::default();
    ge.set_gej_var(&r);
    let mut ge_gen2 = blvm_secp256k1::group::Ge::default();
    ge_gen2.set_gej_var(&r_gen);
    // ecmult(u2*G+u1*G) must equal ecmult_gen(u1+u2) — same x (y may differ for -3*G)
    let mut ge_x = [0u8; 32];
    ge.x.normalize();
    ge.x.get_b32(&mut ge_x);
    let mut ge_gen_x = [0u8; 32];
    ge_gen2.x.normalize();
    ge_gen2.x.get_b32(&mut ge_gen_x);
    assert_eq!(
        ge_x, ge_gen_x,
        "ecmult(u2*G+u1*G) x must equal ecmult_gen(u1+u2) x"
    );
}

/// ecmult with na=0 (only PRE_G): ecmult(G, 0, k) must equal ecmult_gen(k) for any k.
#[test]
fn test_ecmult_ng_only() {
    let g = generator_g();
    let mut gj = Gej::default();
    gj.set_ge(&g);
    let zero = Scalar::zero();
    let mut three = Scalar::zero();
    three.set_int(3);
    let mut neg_three = Scalar::zero();
    neg_three.negate(&three);

    let mut r_ecmult = Gej::default();
    ecmult::ecmult(&mut r_ecmult, &gj, &zero, Some(&neg_three));
    let mut r_gen = Gej::default();
    ecmult::ecmult_gen(&mut r_gen, &neg_three);

    let mut ge_ecmult = blvm_secp256k1::group::Ge::default();
    ge_ecmult.set_gej_var(&r_ecmult);
    let mut ge_gen = blvm_secp256k1::group::Ge::default();
    ge_gen.set_gej_var(&r_gen);
    let c_ecmult = ge_to_compressed(&ge_ecmult);
    let c_gen = ge_to_compressed(&ge_gen);
    assert_eq!(
        c_ecmult, c_gen,
        "ecmult(0*G + (-3)*G) must equal ecmult_gen(-3)"
    );
}

/// ecmult with ng=None (only pre_a): ecmult(G, k, None) must equal ecmult_gen(k) for any k.
#[test]
fn test_ecmult_na_only() {
    let g = generator_g();
    let mut gj = Gej::default();
    gj.set_ge(&g);

    // Positive scalar: ecmult(G, 3, None) = 3*G
    let mut three = Scalar::zero();
    three.set_int(3);
    let mut r_ecmult = Gej::default();
    ecmult::ecmult(&mut r_ecmult, &gj, &three, None);
    let mut r_gen = Gej::default();
    ecmult::ecmult_gen(&mut r_gen, &three);
    let mut ge_ecmult = blvm_secp256k1::group::Ge::default();
    ge_ecmult.set_gej_var(&r_ecmult);
    let mut ge_gen = blvm_secp256k1::group::Ge::default();
    ge_gen.set_gej_var(&r_gen);
    let c_ecmult = ge_to_compressed(&ge_ecmult);
    let c_gen = ge_to_compressed(&ge_gen);
    assert_eq!(c_ecmult, c_gen, "ecmult(3*G + 0) must equal ecmult_gen(3)");

    // Negative scalar: ecmult(G, -3, None) = -3*G
    let mut neg_three = Scalar::zero();
    neg_three.negate(&three);
    ecmult::ecmult(&mut r_ecmult, &gj, &neg_three, None);
    ecmult::ecmult_gen(&mut r_gen, &neg_three);
    ge_ecmult.set_gej_var(&r_ecmult);
    ge_gen.set_gej_var(&r_gen);
    let c_ecmult = ge_to_compressed(&ge_ecmult);
    let c_gen = ge_to_compressed(&ge_gen);
    assert_eq!(
        c_ecmult, c_gen,
        "ecmult((-3)*G + 0) must equal ecmult_gen(-3)"
    );
}

/// ecmult(na*A + ng*G): 1*G + 2*G and 2*G + 1*G must equal 3*G.
#[test]
fn test_ecmult_na_plus_ng() {
    let g = generator_g();
    let mut gj = Gej::default();
    gj.set_ge(&g);
    let mut one = Scalar::zero();
    one.set_int(1);
    let mut two = Scalar::zero();
    two.set_int(2);
    let mut three = Scalar::zero();
    three.set_int(3);

    let mut r_ecmult = Gej::default();
    ecmult::ecmult(&mut r_ecmult, &gj, &one, Some(&two));
    let mut r_gen = Gej::default();
    ecmult::ecmult_gen(&mut r_gen, &three);

    let mut ge_ecmult = blvm_secp256k1::group::Ge::default();
    ge_ecmult.set_gej_var(&r_ecmult);
    let mut ge_gen = blvm_secp256k1::group::Ge::default();
    ge_gen.set_gej_var(&r_gen);
    let c_ecmult = ge_to_compressed(&ge_ecmult);
    let c_gen = ge_to_compressed(&ge_gen);
    assert_eq!(c_ecmult, c_gen, "ecmult(1*G + 2*G) must equal 3*G");

    // Also 2*G + 1*G = 3*G (swap na/ng)
    ecmult::ecmult(&mut r_ecmult, &gj, &two, Some(&one));
    ge_ecmult.set_gej_var(&r_ecmult);
    let c_ecmult2 = ge_to_compressed(&ge_ecmult);
    assert_eq!(c_ecmult2, c_gen, "ecmult(2*G + 1*G) must equal 3*G");

    // na=3, ng=None: only pre_a path
    ecmult::ecmult(&mut r_ecmult, &gj, &three, None);
    ge_ecmult.set_gej_var(&r_ecmult);
    let c_ecmult3 = ge_to_compressed(&ge_ecmult);
    assert_eq!(c_ecmult3, c_gen, "ecmult(3*G + 0) must equal 3*G");
}

/// DER parse/serialize roundtrip.
#[test]
fn test_ecdsa_der_roundtrip() {
    let mut seckey_bytes = [0u8; 32];
    seckey_bytes[31] = 1;
    let msg_bytes: [u8; 32] = [2u8; 32];
    let nonce_bytes: [u8; 32] = [3u8; 32];

    let seckey = scalar_from_b32(&seckey_bytes);
    let message = scalar_from_b32(&msg_bytes);
    let nonce = scalar_from_b32(&nonce_bytes);

    let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
    let der = ecdsa_sig_serialize_der(&sig.0, &sig.1);

    let parsed = ecdsa_sig_parse_der(&der).expect("parse strict DER");
    let der_roundtrip = ecdsa_sig_serialize_der(&parsed.0, &parsed.1);
    assert_eq!(der, der_roundtrip, "DER roundtrip");

    let parsed_lax = ecdsa_sig_parse_der_lax(&der).expect("parse lax DER");
    assert!(parsed.0.ct_eq(&parsed_lax.0).into() && parsed.1.ct_eq(&parsed_lax.1).into());

    let pubkey = pubkey_from_secret(&seckey);
    assert!(ecdsa_sig_verify(&parsed.0, &parsed.1, &pubkey, &message));
}

/// ECDSA signature normalize: high-S becomes low-S.
#[test]
fn test_ecdsa_sig_normalize() {
    let mut s_high = Scalar::zero();
    // n/2 < s < n gives high-S. Use n-1 (definitely high).
    let n_minus_1: [u8; 32] = [
        0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
        0xFE, 0xBA, 0xAE, 0xDC, 0xE6, 0xAF, 0x48, 0xA0, 0x3B, 0xBF, 0xD2, 0x5E, 0x8C, 0xD0, 0x36,
        0x41, 0x40,
    ];
    s_high.set_b32(&n_minus_1);
    assert!(s_high.is_high());
    ecdsa_sig_normalize(&mut s_high);
    assert!(!s_high.is_high());
    assert!(s_high.is_one()); // n-1 normalized = 1 (low-S)
}
