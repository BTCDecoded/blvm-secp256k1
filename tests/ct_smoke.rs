//! Differential and platform smoke checks for constant-time–oriented code paths.
//! This is not a substitute for dudect or hardware trace analysis.

use blvm_secp256k1::ecdsa::{ge_from_compressed, pubkey_from_secret};
use blvm_secp256k1::ecdh;
use blvm_secp256k1::schnorr::{schnorr_sign, schnorr_verify};
use blvm_secp256k1::scalar::Scalar;
use sha2::{Digest, Sha256};
use subtle::ConstantTimeEq;

/// Secret-path signing and ECDH stay consistent for fixed and hashed “random” keys.
#[test]
fn secret_paths_fixed_vs_hashed_scalars() {
    let msg = [0xabu8; 32];
    let aux = [0u8; 32];
    for i in 0u32..32 {
        let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
        let mut sk = [0u8; 32];
        sk.copy_from_slice(&h);
        let mut d = Scalar::zero();
        if d.set_b32(&sk) || d.is_zero() {
            continue;
        }
        let sig = schnorr_sign(&sk, &msg, &aux).expect("schnorr sign");
        let pk = {
            let mut ge = pubkey_from_secret(&d);
            let mut c = [0u8; 33];
            c[0] = if ge.y.is_odd() { 0x03 } else { 0x02 };
            let mut xb = [0u8; 32];
            ge.x.normalize();
            ge.x.get_b32(&mut xb);
            c[1..33].copy_from_slice(&xb);
            c
        };
        let p = ge_from_compressed(&pk.try_into().unwrap()).unwrap();
        let mut xo = [0u8; 32];
        let mut xx = p.x;
        xx.normalize();
        xx.get_b32(&mut xo);
        assert!(schnorr_verify(&sig, &msg, &xo));

        let shared = ecdh::ecdh(&p, &d).expect("ecdh");
        assert_ne!(shared, [0u8; 32]);
    }
}

/// On supported targets, scalar inversion is the safegcd path (see `scalar` module).
#[test]
fn scalar_inv_is_nonzero_on_small_input() {
    let mut a = Scalar::zero();
    a.set_int(0x7fff_ffffu32);
    let mut i = Scalar::zero();
    i.inv(&a);
    let mut prod = Scalar::zero();
    prod.mul(&a, &i);
    let mut one = Scalar::zero();
    one.set_int(1);
    assert!(bool::from(prod.ct_eq(&one)), "a * a^-1 = 1");
}
