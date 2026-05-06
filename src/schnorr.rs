//! BIP 340 Schnorr signatures for secp256k1.
//!
//! 32-byte x-only pubkeys, 64-byte signatures (r || s).
//! Tagged hashes: SHA256(SHA256(tag) || SHA256(tag) || x).
//!
//! **Timing:** `schnorr_sign` and `xonly_pubkey_from_secret` use constant-time `k*G` and pubkey
//! derivation. `schnorr_verify` / batch verify use public-data variable-time `ecmult`.

use sha2::{Digest, Sha256};
use subtle::Choice;

use crate::ecmult;
use crate::ecmult_gen_const;
use crate::field::FieldElement;
use crate::group::{Ge, Gej};
use crate::scalar::Scalar;

/// BIP 340 tagged hash: hash_name(x) = SHA256(SHA256(tag) || SHA256(tag) || x).
/// Pub(crate) for schnorr_half_agg.
pub(crate) fn tagged_hash(tag: &[u8], data: &[u8]) -> [u8; 32] {
    let midstate = tagged_hash_midstate(tag);
    tagged_hash_from_midstate(&midstate, data)
}

/// Precompute SHA256 state after processing SHA256(tag) || SHA256(tag) (one 64-byte block).
/// Clone the returned hasher and feed data to avoid recomputing the tag hash.
fn tagged_hash_midstate(tag: &[u8]) -> Sha256 {
    let tag_hash = Sha256::digest(tag);
    let mut hasher = Sha256::new();
    hasher.update(tag_hash);
    hasher.update(tag_hash);
    hasher
}

#[inline(always)]
fn tagged_hash_from_midstate(midstate: &Sha256, data: &[u8]) -> [u8; 32] {
    let mut h = midstate.clone();
    h.update(data);
    h.finalize().into()
}

use std::sync::OnceLock;

fn challenge_midstate() -> &'static Sha256 {
    static MID: OnceLock<Sha256> = OnceLock::new();
    MID.get_or_init(|| tagged_hash_midstate(b"BIP0340/challenge"))
}

fn nonce_midstate() -> &'static Sha256 {
    static MID: OnceLock<Sha256> = OnceLock::new();
    MID.get_or_init(|| tagged_hash_midstate(b"BIP0340/nonce"))
}

fn aux_midstate() -> &'static Sha256 {
    static MID: OnceLock<Sha256> = OnceLock::new();
    MID.get_or_init(|| tagged_hash_midstate(b"BIP0340/aux"))
}

/// BIP0340/challenge tagged hash with precomputed midstate. Pub(crate) for musig.
#[inline(always)]
pub(crate) fn tagged_hash_challenge(data: &[u8]) -> [u8; 32] {
    tagged_hash_from_midstate(challenge_midstate(), data)
}

/// lift_x: Given 32-byte x (field element), return point P with even y or None.
/// Fails if x is not on curve or x >= p.
#[inline(always)]
pub fn lift_x(x_bytes: &[u8; 32]) -> Option<Ge> {
    let mut x = FieldElement::zero();
    if !x.set_b32_limit(x_bytes) {
        return None;
    }
    let mut p = Ge::default();
    if !p.set_xo_var(&x, false) {
        return None;
    }
    Some(p)
}

/// x-only pubkey bytes from affine point (32 bytes = x coordinate).
/// Caller must ensure P has even y.
fn ge_to_xonly(p: &Ge) -> [u8; 32] {
    let mut x = p.x;
    x.normalize();
    let mut out = [0u8; 32];
    x.get_b32(&mut out);
    out
}

/// BIP 340 x-only pubkey from secret key.
/// Returns 32-byte pubkey (x coordinate of d'·G with even y).
pub fn xonly_pubkey_from_secret(seckey: &[u8; 32]) -> Option<[u8; 32]> {
    let mut d = Scalar::zero();
    if d.set_b32(seckey) {
        return None;
    }
    if d.is_zero() {
        return None;
    }
    let mut pj = Gej::default();
    ecmult_gen_const(&mut pj, &d);
    let mut p = Ge::default();
    p.set_gej_var(&pj);
    p.x.normalize();
    // x-only encoding depends only on x (same for (x,y) and (x,-y)).
    let mut out = [0u8; 32];
    p.x.get_b32(&mut out);
    Some(out)
}

/// BIP 340 Schnorr verify.
/// sig64: 64-byte signature (r || s)
/// msg: message bytes (arbitrary length)
/// pubkey_x32: 32-byte x-only pubkey
#[inline(always)]
pub fn schnorr_verify(sig64: &[u8; 64], msg: &[u8], pubkey_x32: &[u8; 32]) -> bool {
    let r_bytes: [u8; 32] = sig64[0..32].try_into().unwrap();
    let s_bytes: [u8; 32] = sig64[32..64].try_into().unwrap();

    let mut r_fe = FieldElement::zero();
    if !r_fe.set_b32_limit(&r_bytes) {
        return false; // r >= p
    }

    let mut s = Scalar::zero();
    if s.set_b32(&s_bytes) {
        return false; // s >= n
    }

    let p = match lift_x(pubkey_x32) {
        Some(p) => p,
        None => return false,
    };

    let e_hash = {
        let mut h = challenge_midstate().clone();
        h.update(r_bytes);
        h.update(pubkey_x32);
        h.update(msg);
        let result: [u8; 32] = h.finalize().into();
        result
    };

    let mut e = Scalar::zero();
    e.set_b32(&e_hash);

    let mut neg_e = Scalar::zero();
    neg_e.negate(&e);

    let mut pj = Gej::default();
    pj.set_ge(&p);

    // R = s*G + (-e)*P
    let mut rj = Gej::default();
    ecmult::ecmult(&mut rj, &pj, &neg_e, Some(&s));

    if rj.infinity {
        return false;
    }

    let mut r_computed = Ge::default();
    r_computed.set_gej_var(&rj);
    r_computed.x.normalize();
    r_computed.y.normalize();

    if r_computed.y.is_odd() {
        return false;
    }

    let mut r_computed_bytes = [0u8; 32];
    r_computed.x.get_b32(&mut r_computed_bytes);

    r_computed_bytes == r_bytes
}

/// BIP 340 Schnorr batch verify.
/// Returns true if all signatures are valid, false if any is invalid.
/// Uses ecmult_multi for both sum_left and sum_right.
pub fn schnorr_verify_batch(sigs: &[[u8; 64]], msgs: &[&[u8]], pubkeys: &[[u8; 32]]) -> bool {
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n == 0 {
        return true;
    }
    if n == 1 {
        return schnorr_verify(&sigs[0], msgs[0], &pubkeys[0]);
    }

    // Try batch algorithm: sum_left == sum_right
    match schnorr_verify_batch_inner(sigs, msgs, pubkeys, n) {
        Some(ok) => ok,
        None => (0..n).all(|i| schnorr_verify(&sigs[i], msgs[i], &pubkeys[i])),
    }
}

/// Inner batch verify. Returns Some(true/false) on success, None on parse failure.
/// Uses stack allocation for n <= 16 (all buffers) or n <= 64 (all_scalars/all_points only).
fn schnorr_verify_batch_inner(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
) -> Option<bool> {
    if n <= 16 {
        schnorr_verify_batch_stack_16(sigs, msgs, pubkeys, n)
    } else if n <= 64 {
        schnorr_verify_batch_stack_64(sigs, msgs, pubkeys, n)
    } else {
        schnorr_verify_batch_heap(sigs, msgs, pubkeys, n)
    }
}

struct BatchParseState<'a> {
    s_scalars: &'a mut [Scalar],
    e_scalars: &'a mut [Scalar],
    pk_points: &'a mut [Ge],
    r_points: &'a mut [Ge],
    randoms: &'a mut [Scalar],
}

fn schnorr_verify_batch_parse(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
    state: &mut BatchParseState<'_>,
) -> Option<()> {
    for i in 0..n {
        let r_bytes: [u8; 32] = sigs[i][0..32].try_into().unwrap();
        let s_bytes: [u8; 32] = sigs[i][32..64].try_into().unwrap();

        let mut r_fe = FieldElement::zero();
        if !r_fe.set_b32_limit(&r_bytes) {
            return None;
        }

        let mut s = Scalar::zero();
        if s.set_b32(&s_bytes) {
            return None;
        }
        state.s_scalars[i] = s;

        state.pk_points[i] = lift_x(&pubkeys[i])?;

        let e_hash = {
            let mut h = challenge_midstate().clone();
            h.update(r_bytes);
            h.update(pubkeys[i]);
            h.update(msgs[i]);
            let result: [u8; 32] = h.finalize().into();
            result
        };

        let mut e = Scalar::zero();
        e.set_b32(&e_hash);
        state.e_scalars[i] = e;

        state.r_points[i] = lift_x(&r_bytes)?;

        let mut random = Scalar::zero();
        let _ = random.set_b32(&r_bytes);
        let mut random_tmp = Scalar::zero();
        random_tmp.add(&random, &e);
        state.randoms[i] = random_tmp;
    }
    Some(())
}

fn schnorr_verify_batch_stack_16(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
) -> Option<bool> {
    let mut s_scalars = [Scalar::zero(); 16];
    let mut e_scalars = [Scalar::zero(); 16];
    let mut pk_points = [Ge::default(); 16];
    let mut r_points = [Ge::default(); 16];
    let mut randoms = [Scalar::zero(); 16];
    schnorr_verify_batch_parse(
        sigs,
        msgs,
        pubkeys,
        n,
        &mut BatchParseState {
            s_scalars: &mut s_scalars,
            e_scalars: &mut e_scalars,
            pk_points: &mut pk_points,
            r_points: &mut r_points,
            randoms: &mut randoms,
        },
    )?;

    let mut g_scalar = Scalar::zero();
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&randoms[i], &s_scalars[i]);
        let mut g_new = Scalar::zero();
        g_new.add(&g_scalar, &term);
        g_scalar = g_new;
    }

    let mut all_scalars = [Scalar::zero(); 32];
    let mut all_points = [Ge::default(); 32];
    for i in 0..n {
        let mut neg_e = Scalar::zero();
        neg_e.negate(&e_scalars[i]);
        let mut z_neg_e = Scalar::zero();
        z_neg_e.mul(&randoms[i], &neg_e);
        all_scalars[i] = z_neg_e;
        all_points[i] = pk_points[i];
    }
    for i in 0..n {
        let mut neg_rand = Scalar::zero();
        neg_rand.negate(&randoms[i]);
        all_scalars[n + i] = neg_rand;
        all_points[n + i] = r_points[i];
    }

    let mut result = Gej::default();
    ecmult::ecmult_multi(
        &mut result,
        &g_scalar,
        &all_scalars[..2 * n],
        &all_points[..2 * n],
    );
    Some(result.is_infinity())
}

fn schnorr_verify_batch_stack_64(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
) -> Option<bool> {
    let mut s_scalars = Vec::with_capacity(n);
    let mut e_scalars = Vec::with_capacity(n);
    let mut pk_points = Vec::with_capacity(n);
    let mut r_points = Vec::with_capacity(n);
    let mut randoms = Vec::with_capacity(n);
    s_scalars.resize(n, Scalar::zero());
    e_scalars.resize(n, Scalar::zero());
    pk_points.resize(n, Ge::default());
    r_points.resize(n, Ge::default());
    randoms.resize(n, Scalar::zero());
    schnorr_verify_batch_parse(
        sigs,
        msgs,
        pubkeys,
        n,
        &mut BatchParseState {
            s_scalars: &mut s_scalars,
            e_scalars: &mut e_scalars,
            pk_points: &mut pk_points,
            r_points: &mut r_points,
            randoms: &mut randoms,
        },
    )?;

    let mut g_scalar = Scalar::zero();
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&randoms[i], &s_scalars[i]);
        let mut g_new = Scalar::zero();
        g_new.add(&g_scalar, &term);
        g_scalar = g_new;
    }

    let mut all_scalars = [Scalar::zero(); 128];
    let mut all_points = [Ge::default(); 128];
    for i in 0..n {
        let mut neg_e = Scalar::zero();
        neg_e.negate(&e_scalars[i]);
        let mut z_neg_e = Scalar::zero();
        z_neg_e.mul(&randoms[i], &neg_e);
        all_scalars[i] = z_neg_e;
        all_points[i] = pk_points[i];
    }
    for i in 0..n {
        let mut neg_rand = Scalar::zero();
        neg_rand.negate(&randoms[i]);
        all_scalars[n + i] = neg_rand;
        all_points[n + i] = r_points[i];
    }

    let mut result = Gej::default();
    ecmult::ecmult_multi(
        &mut result,
        &g_scalar,
        &all_scalars[..2 * n],
        &all_points[..2 * n],
    );
    Some(result.is_infinity())
}

fn schnorr_verify_batch_heap(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
) -> Option<bool> {
    let mut s_scalars = Vec::with_capacity(n);
    let mut e_scalars = Vec::with_capacity(n);
    let mut pk_points = Vec::with_capacity(n);
    let mut r_points = Vec::with_capacity(n);
    let mut randoms = Vec::with_capacity(n);

    for i in 0..n {
        let r_bytes: [u8; 32] = sigs[i][0..32].try_into().unwrap();
        let s_bytes: [u8; 32] = sigs[i][32..64].try_into().unwrap();

        let mut r_fe = FieldElement::zero();
        if !r_fe.set_b32_limit(&r_bytes) {
            return None;
        }

        let mut s = Scalar::zero();
        if s.set_b32(&s_bytes) {
            return None;
        }
        s_scalars.push(s);

        let p = lift_x(&pubkeys[i])?;
        pk_points.push(p);

        let e_hash = {
            let mut h = challenge_midstate().clone();
            h.update(r_bytes);
            h.update(pubkeys[i]);
            h.update(msgs[i]);
            let result: [u8; 32] = h.finalize().into();
            result
        };

        let mut e = Scalar::zero();
        e.set_b32(&e_hash);
        e_scalars.push(e);

        let r_ge = lift_x(&r_bytes)?;
        r_points.push(r_ge);

        let mut random = Scalar::zero();
        let _ = random.set_b32(&r_bytes);
        let mut random_tmp = Scalar::zero();
        random_tmp.add(&random, &e);
        randoms.push(random_tmp);
    }

    let mut g_scalar = Scalar::zero();
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&randoms[i], &s_scalars[i]);
        let mut g_new = Scalar::zero();
        g_new.add(&g_scalar, &term);
        g_scalar = g_new;
    }

    let mut all_scalars = Vec::with_capacity(2 * n);
    let mut all_points = Vec::with_capacity(2 * n);

    for i in 0..n {
        let mut neg_e = Scalar::zero();
        neg_e.negate(&e_scalars[i]);
        let mut z_neg_e = Scalar::zero();
        z_neg_e.mul(&randoms[i], &neg_e);
        all_scalars.push(z_neg_e);
        all_points.push(pk_points[i]);
    }
    for i in 0..n {
        let mut neg_rand = Scalar::zero();
        neg_rand.negate(&randoms[i]);
        all_scalars.push(neg_rand);
        all_points.push(r_points[i]);
    }

    let mut result = Gej::default();
    ecmult::ecmult_multi(&mut result, &g_scalar, &all_scalars, &all_points);
    Some(result.is_infinity())
}

/// BIP 340 Schnorr sign.
/// seckey: 32-byte secret key
/// msg: message bytes (arbitrary length)
/// aux_rand32: 32-byte auxiliary randomness (use zeros if no randomness)
pub fn schnorr_sign(seckey: &[u8; 32], msg: &[u8], aux_rand32: &[u8; 32]) -> Option<[u8; 64]> {
    let mut d = Scalar::zero();
    if d.set_b32(seckey) {
        return None;
    }
    if d.is_zero() {
        return None;
    }

    // P = d*G — d is secret, must use the constant-time path. (Previous implementation
    // mistakenly used the var-time `ecmult::ecmult_gen`, which leaked timing on `d`.)
    let mut pj = Gej::default();
    ecmult_gen_const(&mut pj, &d);
    let mut p_ge = Ge::default();
    p_ge.set_gej_var(&pj);
    p_ge.x.normalize();
    p_ge.y.normalize();
    let pk_parity_odd = p_ge.y.is_odd();
    let pk = ge_to_xonly(&p_ge);

    schnorr_sign_inner(seckey, &d, pk_parity_odd, &pk, msg, aux_rand32)
}

/// BIP 340 Schnorr keypair: seckey + cached x-only pubkey (with parity already resolved).
/// Equivalent to libsecp256k1's `secp256k1_keypair` and lets callers reuse one CT `d*G`
/// across many signatures.
#[derive(Clone, Copy, Debug)]
pub struct Keypair {
    /// **Original** secret key bytes (NOT y-flipped). The y-flip needed to produce a
    /// signing-domain `d_adj` is decided by `pk_parity_odd` and re-applied per sign call.
    /// Storing the original bytes avoids exposing the y-flip transformation to callers.
    pub(crate) seckey: [u8; 32],
    /// X-only pubkey (32 bytes) matching `d*G` lifted to the even-y representative.
    pub(crate) pubkey_xonly: [u8; 32],
    /// Whether `(d*G).y` was odd before the BIP340 even-y normalization (so signing
    /// must use `-d` instead of `d`). One bit, public; carries no secret data.
    pub(crate) pk_parity_odd: bool,
}

impl Keypair {
    /// Derive `Keypair` from a 32-byte secret key. Performs **one** CT `d*G` and caches
    /// the pubkey + parity so subsequent `schnorr_sign_with_keypair` calls only need the
    /// `k*G` for the nonce (matches libsecp's `keypair_create` + `schnorrsig_sign32`).
    pub fn from_seckey(seckey32: &[u8; 32]) -> Option<Self> {
        let mut d = Scalar::zero();
        if d.set_b32(seckey32) {
            return None;
        }
        if d.is_zero() {
            return None;
        }
        let mut pj = Gej::default();
        ecmult_gen_const(&mut pj, &d);
        let mut p = Ge::default();
        p.set_gej_var(&pj);
        p.x.normalize();
        p.y.normalize();
        let pk_parity_odd = p.y.is_odd();
        let pubkey_xonly = ge_to_xonly(&p);
        Some(Keypair {
            seckey: *seckey32,
            pubkey_xonly,
            pk_parity_odd,
        })
    }

    /// X-only pubkey (32 bytes) from this keypair.
    #[inline]
    pub fn xonly_pubkey(&self) -> [u8; 32] {
        self.pubkey_xonly
    }
}

/// BIP 340 Schnorr sign using a precomputed [`Keypair`]. Per-call cost is **one** CT
/// `k*G` (the nonce point) plus hashing, mirroring libsecp's keypair-based sign loop.
///
/// Constant-time w.r.t. the secret key inside `keypair`.
pub fn schnorr_sign_with_keypair(
    keypair: &Keypair,
    msg: &[u8],
    aux_rand32: &[u8; 32],
) -> Option<[u8; 64]> {
    let mut d = Scalar::zero();
    if d.set_b32(&keypair.seckey) {
        return None;
    }
    if d.is_zero() {
        return None;
    }
    schnorr_sign_inner(
        &keypair.seckey,
        &d,
        keypair.pk_parity_odd,
        &keypair.pubkey_xonly,
        msg,
        aux_rand32,
    )
}

#[inline]
fn schnorr_sign_inner(
    seckey: &[u8; 32],
    d: &Scalar,
    pk_parity_odd: bool,
    pk: &[u8; 32],
    msg: &[u8],
    aux_rand32: &[u8; 32],
) -> Option<[u8; 64]> {
    let mut d_adj = *d;
    if pk_parity_odd {
        d_adj.negate(d);
    }

    let aux_hash = tagged_hash_from_midstate(aux_midstate(), aux_rand32);
    let mut masked_key = [0u8; 32];
    for i in 0..32 {
        masked_key[i] = seckey[i] ^ aux_hash[i];
    }

    let k_hash = {
        let mut h = nonce_midstate().clone();
        h.update(masked_key);
        h.update(*pk);
        h.update(msg);
        let result: [u8; 32] = h.finalize().into();
        result
    };

    let mut k = Scalar::zero();
    if k.set_b32(&k_hash) {
        return None;
    }
    if k.is_zero() {
        return None;
    }

    // Single constant-time ecmult for nonce point R = k*G
    let mut rj = Gej::default();
    ecmult_gen_const(&mut rj, &k);
    let mut r_ge = Ge::default();
    r_ge.set_gej_var(&rj);
    r_ge.x.normalize();
    r_ge.y.normalize();

    // Branchless: negate k and r.y if r.y is odd (BIP340 requires even-y nonce point).
    // Checking/branching on r_ge.y.is_odd() would leak nonce-dependent timing.
    let parity = r_ge.y.is_odd() as i32;
    k.cond_negate(parity);
    let mut neg_ry = FieldElement::zero();
    neg_ry.negate(&r_ge.y, 1);
    neg_ry.normalize();
    r_ge.y.normalize();
    r_ge.y.cmov(&neg_ry, Choice::from(parity as u8));

    let mut r_bytes = [0u8; 32];
    r_ge.x.get_b32(&mut r_bytes);

    let e_hash = {
        let mut h = challenge_midstate().clone();
        h.update(r_bytes);
        h.update(*pk);
        h.update(msg);
        let result: [u8; 32] = h.finalize().into();
        result
    };

    let mut e = Scalar::zero();
    e.set_b32(&e_hash);

    // s = k + e*d_adj
    let mut ed = Scalar::zero();
    ed.mul(&e, &d_adj);
    let mut s = Scalar::zero();
    s.add(&k, &ed);

    let mut s_bytes = [0u8; 32];
    s.get_b32(&mut s_bytes);

    let mut sig = [0u8; 64];
    sig[0..32].copy_from_slice(&r_bytes);
    sig[32..64].copy_from_slice(&s_bytes);

    Some(sig)
}
