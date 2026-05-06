//! ECDH: Elliptic Curve Diffie-Hellman shared secret.
//!
//! Computes shared_point = pubkey * seckey, then hashes (x,y) to output.
//! Uses SHA256(0x02|0x03 || x32) as default (compressed point hash).
//!
//! **Timing:** `seckey` is multiplied with constant-time `ecmult_const` (libsecp256k1-style).

use crate::ecdsa::ge_from_compressed;
use crate::ecmult_const;
use crate::group::{Ge, Gej};
use crate::scalar::Scalar;
use sha2::{Digest, Sha256};

/// ECDH hash: SHA256(version || x32) where version = (y_odd as u8) | 0x02.
/// Matches libsecp256k1 ecdh_hash_function_sha256.
pub fn ecdh_hash_sha256(output: &mut [u8; 32], x: &[u8; 32], y: &[u8; 32]) {
    let version = (y[31] & 0x01) | 0x02;
    let mut hasher = Sha256::new();
    hasher.update([version]);
    hasher.update(x);
    hasher.finalize_into((&mut output[..]).into());
}

/// Compute ECDH shared secret: output = hash(pubkey * seckey).
/// Returns None if seckey invalid (zero or overflow).
/// Uses default SHA256 hash (compressed point).
pub fn ecdh(pubkey: &Ge, seckey: &Scalar) -> Option<[u8; 32]> {
    if seckey.is_zero() {
        return None;
    }
    if pubkey.infinity {
        return None;
    }
    let mut pubkeyj = Gej::default();
    pubkeyj.set_ge(pubkey);
    let mut res = Gej::default();
    ecmult_const(&mut res, pubkey, seckey);
    let mut pt = Ge::default();
    pt.set_gej_var(&res);
    if pt.infinity {
        return None;
    }
    pt.x.normalize();
    pt.y.normalize();
    let mut x = [0u8; 32];
    let mut y = [0u8; 32];
    pt.x.get_b32(&mut x);
    pt.y.get_b32(&mut y);
    let mut output = [0u8; 32];
    ecdh_hash_sha256(&mut output, &x, &y);
    Some(output)
}

/// ECDH from compressed pubkey bytes. Returns None if invalid.
pub fn ecdh_compressed(pubkey_bytes: &[u8; 33], seckey_bytes: &[u8; 32]) -> Option<[u8; 32]> {
    let pubkey = ge_from_compressed(pubkey_bytes)?;
    let mut seckey = Scalar::zero();
    if seckey.set_b32(seckey_bytes) {
        return None;
    }
    if seckey.is_zero() {
        return None;
    }
    ecdh(&pubkey, &seckey)
}
