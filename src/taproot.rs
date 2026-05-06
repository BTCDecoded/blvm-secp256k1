//! BIP 341 Taproot utilities.
//!
//! Standalone functions for Taproot key path and script path validation:
//! - Tagged hashes (TapTweak, TapLeaf, TapBranch)
//! - x-only pubkey tweak (internal key + merkle root → output key)
//! - pubkey_combine for OP_CHECKSIGADD and multi-key scripts
//!
//! **Timing:** `t*G` in `xonly_pubkey_tweak_add` uses constant-time `ecmult_gen_const` so the
//! tweak scalar does not go through the fast variable-time generator ladder.

use crate::ecdsa::{ge_from_compressed, ge_to_compressed};
use crate::ecmult_gen_const;
use crate::group::{Ge, Gej};
use crate::scalar::Scalar;
use crate::schnorr::{lift_x, tagged_hash};

// BIP 341 tagged hash tags
const TAG_TAP_TWEAK: &[u8] = b"TapTweak";
const TAG_TAP_LEAF: &[u8] = b"TapLeaf";
const TAG_TAP_BRANCH: &[u8] = b"TapBranch";
const TAG_TAP_SIGHASH: &[u8] = b"TapSighash";

/// BIP 341 TapTweak hash: tag "TapTweak", data = internal_key || merkle_root.
pub fn tap_tweak_hash(internal_key: &[u8; 32], merkle_root: &[u8; 32]) -> [u8; 32] {
    let mut data = Vec::with_capacity(64);
    data.extend_from_slice(internal_key);
    data.extend_from_slice(merkle_root);
    tagged_hash(TAG_TAP_TWEAK, &data)
}

/// BIP 341 TapLeaf hash: tag "TapLeaf", data = leaf_version || compact_size(script_len) || script.
pub fn tap_leaf_hash(leaf_version: u8, script: &[u8]) -> [u8; 32] {
    let mut data = Vec::with_capacity(1 + compact_size_len(script.len()) + script.len());
    data.push(leaf_version);
    data.extend(compact_size_encode(script.len()));
    data.extend_from_slice(script);
    tagged_hash(TAG_TAP_LEAF, &data)
}

/// BIP 341 TapBranch hash: tag "TapBranch", data = left || right (caller must sort).
pub fn tap_branch_hash(left: &[u8; 32], right: &[u8; 32]) -> [u8; 32] {
    let mut data = Vec::with_capacity(64);
    data.extend_from_slice(left);
    data.extend_from_slice(right);
    tagged_hash(TAG_TAP_BRANCH, &data)
}

/// BIP 341 TapSighash: tag "TapSighash", hashes the SigMsg for Taproot/Tapscript verification.
pub fn tap_sighash_hash(data: &[u8]) -> [u8; 32] {
    tagged_hash(TAG_TAP_SIGHASH, data)
}

fn compact_size_len(n: usize) -> usize {
    if n < 0xfd {
        1
    } else if n <= 0xffff {
        3
    } else if n <= 0xffff_ffff {
        5
    } else {
        9
    }
}

fn compact_size_encode(n: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(9);
    if n < 0xfd {
        out.push(n as u8);
    } else if n <= 0xffff {
        out.push(0xfd);
        out.extend_from_slice(&(n as u16).to_le_bytes());
    } else if n <= 0xffff_ffff {
        out.push(0xfe);
        out.extend_from_slice(&(n as u32).to_le_bytes());
    } else {
        out.push(0xff);
        out.extend_from_slice(&(n as u64).to_le_bytes());
    }
    out
}

/// x-only pubkey from affine point. Returns (x_bytes, parity) where parity is 0 for even y, 1 for odd.
pub fn xonly_from_point(ge: &Ge) -> ([u8; 32], u8) {
    if ge.infinity {
        return ([0u8; 32], 0);
    }
    let mut x = ge.x;
    x.normalize();
    let mut y = ge.y;
    y.normalize();
    let mut out = [0u8; 32];
    x.get_b32(&mut out);
    let parity = if y.is_odd() { 1 } else { 0 };
    (out, parity)
}

/// Standalone x-only pubkey tweak: Q = P + t*G (or -P + t*G for x-only).
/// Returns (output_xonly, parity). parity = 0 if Q has even y, 1 if odd.
pub fn xonly_pubkey_tweak_add(xonly: &[u8; 32], tweak: &[u8; 32]) -> Option<([u8; 32], u8)> {
    let p = lift_x(xonly)?;
    let mut t = Scalar::zero();
    if t.set_b32(tweak) {
        return None; // overflow: t >= n
    }

    let mut tj = Gej::default();
    ecmult_gen_const(&mut tj, &t);
    let mut pj = Gej::default();
    pj.set_ge(&p);
    let pj_old = pj;
    pj.add_var(&pj_old, &tj);

    if pj.is_infinity() {
        return None;
    }

    let mut q = Ge::default();
    q.set_gej_var(&pj);
    q.x.normalize();
    q.y.normalize();

    Some(xonly_from_point(&q))
}

/// Compute Taproot output key from internal key and merkle root.
/// Returns output_key (32 bytes) or None if invalid.
pub fn taproot_output_key(internal_key: &[u8; 32], merkle_root: &[u8; 32]) -> Option<[u8; 32]> {
    let tweak = tap_tweak_hash(internal_key, merkle_root);
    let (output, _parity) = xonly_pubkey_tweak_add(internal_key, &tweak)?;
    Some(output)
}

/// Combine multiple compressed pubkeys via EC point addition.
/// Returns None if any pubkey is invalid or result is infinity.
pub fn pubkey_combine(pubkeys: &[[u8; 33]]) -> Option<[u8; 33]> {
    if pubkeys.is_empty() {
        return None;
    }
    let mut sum = ge_from_compressed(&pubkeys[0])?;
    for pk in pubkeys.iter().skip(1) {
        let p = ge_from_compressed(pk)?;
        let mut sumj = Gej::default();
        sumj.set_ge(&sum);
        let mut pj = Gej::default();
        pj.set_ge(&p);
        let sumj_old = sumj;
        sumj.add_var(&sumj_old, &pj);
        if sumj.is_infinity() {
            return None;
        }
        sum.set_gej_var(&sumj);
    }
    sum.x.normalize();
    sum.y.normalize();
    Some(ge_to_compressed(&sum))
}
