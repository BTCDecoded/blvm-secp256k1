//! ECDSA signatures for secp256k1.
//!
//! Verify and sign (r,s) against message hash.
//! Low-S enforcement; compact (64-byte) and DER parsing.

use crate::ecmult;
use crate::field::FieldElement;
use crate::group::{generator_g, Ge, Gej};
use crate::scalar::Scalar;
use subtle::ConstantTimeEq;

use std::sync::OnceLock;

/// Group order n as field element (for ECDSA xr + n check).
fn ecdsa_const_order_as_fe() -> &'static FieldElement {
    static FE: OnceLock<FieldElement> = OnceLock::new();
    FE.get_or_init(|| {
        let n: [u8; 32] = [
            0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
            0xFF, 0xFE, 0xBA, 0xAE, 0xDC, 0xE6, 0xAF, 0x48, 0xA0, 0x3B, 0xBF, 0xD2, 0x5E, 0x8C,
            0xD0, 0x36, 0x41, 0x41,
        ];
        let mut r = FieldElement::zero();
        r.set_b32_mod(&n);
        r
    })
}

/// p - n (for ECDSA: if xr >= p-n then xr + n >= p, skip second case).
fn ecdsa_const_p_minus_order() -> &'static FieldElement {
    static FE: OnceLock<FieldElement> = OnceLock::new();
    FE.get_or_init(|| {
        let p_minus_n: [u8; 32] = [
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x01, 0x45, 0x51, 0x23, 0x19, 0x50, 0xB7, 0x5F, 0xC4, 0x40, 0x2D, 0xA1, 0x72,
            0x2F, 0xC9, 0xBA, 0xEE,
        ];
        let mut r = FieldElement::zero();
        r.set_b32_mod(&p_minus_n);
        r
    })
}

/// Batch verify ECDSA signatures. Returns true if all valid, false if any invalid.
/// Uses random linear combination: sum z_i*(u1_i*G + u2_i*P_i) == sum z_i*R_i.
/// z_i = r_i + s_i (deterministic from sig). Rejects high-S.
pub fn ecdsa_verify_batch(sigs: &[[u8; 64]], msgs: &[[u8; 32]], pubkeys: &[[u8; 33]]) -> bool {
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n == 0 {
        return true;
    }
    if n == 1 {
        if let (Some((sigr, sigs_scalar)), Some(pk), Some(msg)) = (
            ecdsa_sig_parse_compact(&sigs[0]),
            ge_from_compressed(&pubkeys[0]),
            {
                let mut m = Scalar::zero();
                if m.set_b32(&msgs[0]) {
                    None
                } else {
                    Some(m)
                }
            },
        ) {
            return ecdsa_sig_verify(&sigr, &sigs_scalar, &pk, &msg);
        }
        return false;
    }

    if n <= 64 {
        ecdsa_verify_batch_stack(sigs, msgs, pubkeys, n)
    } else {
        ecdsa_verify_batch_heap(sigs, msgs, pubkeys, n)
    }
}

fn ecdsa_verify_batch_stack(
    sigs: &[[u8; 64]],
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
    n: usize,
) -> bool {
    use crate::group::Gej;

    let mut u1_vec = [Scalar::zero(); 64];
    let mut u2_vec = [Scalar::zero(); 64];
    let mut pk_vec = [Ge::default(); 64];
    let mut r_pt_vec = [Ge::default(); 64];
    let mut z_vec = [Scalar::zero(); 64];

    for i in 0..n {
        let (sigr, sigs_scalar) = match ecdsa_sig_parse_compact(&sigs[i]) {
            Some(s) => s,
            None => return false,
        };
        if sigr.is_zero() || sigs_scalar.is_zero() || sigs_scalar.is_high() {
            return false;
        }
        let pk = match ge_from_compressed(&pubkeys[i]) {
            Some(p) => p,
            None => return false,
        };
        let mut m = Scalar::zero();
        if m.set_b32(&msgs[i]) {
            return false;
        }

        let mut sn = Scalar::zero();
        sn.inv_var(&sigs_scalar);
        let mut u1 = Scalar::zero();
        u1.mul(&sn, &m);
        let mut u2 = Scalar::zero();
        u2.mul(&sn, &sigr);
        u1_vec[i] = u1;
        u2_vec[i] = u2;
        pk_vec[i] = pk;

        let mut r_buf = [0u8; 32];
        sigr.get_b32(&mut r_buf);
        let mut rx_fe = FieldElement::zero();
        if !rx_fe.set_b32_limit(&r_buf) {
            return false;
        }
        let mut r_ge = Ge::default();
        if !r_ge.set_xo_var(&rx_fe, false) {
            return false;
        }
        r_pt_vec[i] = r_ge;

        let mut z = Scalar::zero();
        z.add(&sigr, &sigs_scalar);
        z_vec[i] = z;
    }

    let mut g_scalar = Scalar::zero();
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&z_vec[i], &u1_vec[i]);
        let mut g_new = Scalar::zero();
        g_new.add(&g_scalar, &term);
        g_scalar = g_new;
    }

    let mut left_scalars = [Scalar::zero(); 64];
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&z_vec[i], &u2_vec[i]);
        left_scalars[i] = term;
    }

    let mut sum_left = Gej::default();
    ecmult::ecmult_multi(&mut sum_left, &g_scalar, &left_scalars[..n], &pk_vec[..n]);

    let mut sum_right = Gej::default();
    ecmult::ecmult_multi(&mut sum_right, &Scalar::zero(), &z_vec[..n], &r_pt_vec[..n]);

    let mut neg_right = Gej::default();
    neg_right.neg(&sum_right);
    let mut diff = Gej::default();
    diff.add_var(&sum_left, &neg_right);
    diff.is_infinity()
}

fn ecdsa_verify_batch_heap(
    sigs: &[[u8; 64]],
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
    n: usize,
) -> bool {
    use crate::group::Gej;

    let mut u1_vec = Vec::with_capacity(n);
    let mut u2_vec = Vec::with_capacity(n);
    let mut pk_vec = Vec::with_capacity(n);
    let mut r_pt_vec = Vec::with_capacity(n);
    let mut z_vec = Vec::with_capacity(n);

    for i in 0..n {
        let (sigr, sigs_scalar) = match ecdsa_sig_parse_compact(&sigs[i]) {
            Some(s) => s,
            None => return false,
        };
        if sigr.is_zero() || sigs_scalar.is_zero() || sigs_scalar.is_high() {
            return false;
        }
        let pk = match ge_from_compressed(&pubkeys[i]) {
            Some(p) => p,
            None => return false,
        };
        let mut m = Scalar::zero();
        if m.set_b32(&msgs[i]) {
            return false;
        }

        let mut sn = Scalar::zero();
        sn.inv_var(&sigs_scalar);
        let mut u1 = Scalar::zero();
        u1.mul(&sn, &m);
        let mut u2 = Scalar::zero();
        u2.mul(&sn, &sigr);
        u1_vec.push(u1);
        u2_vec.push(u2);
        pk_vec.push(pk);

        let mut r_buf = [0u8; 32];
        sigr.get_b32(&mut r_buf);
        let mut rx_fe = FieldElement::zero();
        if !rx_fe.set_b32_limit(&r_buf) {
            return false;
        }
        let mut r_ge = Ge::default();
        if !r_ge.set_xo_var(&rx_fe, false) {
            return false;
        }
        r_pt_vec.push(r_ge);

        let mut z = Scalar::zero();
        z.add(&sigr, &sigs_scalar);
        z_vec.push(z);
    }

    let mut g_scalar = Scalar::zero();
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&z_vec[i], &u1_vec[i]);
        let mut g_new = Scalar::zero();
        g_new.add(&g_scalar, &term);
        g_scalar = g_new;
    }

    let mut left_scalars = Vec::with_capacity(n);
    for i in 0..n {
        let mut term = Scalar::zero();
        term.mul(&z_vec[i], &u2_vec[i]);
        left_scalars.push(term);
    }

    let mut sum_left = Gej::default();
    ecmult::ecmult_multi(&mut sum_left, &g_scalar, &left_scalars, &pk_vec);

    let mut sum_right = Gej::default();
    ecmult::ecmult_multi(&mut sum_right, &Scalar::zero(), &z_vec, &r_pt_vec);

    let mut neg_right = Gej::default();
    neg_right.neg(&sum_right);
    let mut diff = Gej::default();
    diff.add_var(&sum_left, &neg_right);
    diff.is_infinity()
}

/// Verify ECDSA signature (r,s) against pubkey and message hash.
/// Returns true if valid. Does not enforce low-S (accepts both).
pub fn ecdsa_sig_verify(sigr: &Scalar, sigs: &Scalar, pubkey: &Ge, message: &Scalar) -> bool {
    if sigr.is_zero() || sigs.is_zero() {
        return false;
    }

    let mut sn = Scalar::zero();
    sn.inv_var(sigs);
    let mut u1 = Scalar::zero();
    u1.mul(&sn, message);
    let mut u2 = Scalar::zero();
    u2.mul(&sn, sigr);

    let mut pubkeyj = Gej::default();
    pubkeyj.set_ge(pubkey);
    let mut pr = Gej::default();
    ecmult::ecmult(&mut pr, &pubkeyj, &u2, Some(&u1));

    if pr.infinity {
        return false;
    }

    // sigr bytes -> xr (field element). libsecp256k1: scalar_get_b32 then fe_set_b32_limit.
    let mut c = [0u8; 32];
    sigr.get_b32(&mut c);
    let mut xr = FieldElement::zero();
    if !xr.set_b32_limit(&c) {
        return false;
    }

    // Check xr == X(pr) mod n. Since 2*n > p, h in {0,1}.
    // Case 1: xr == pr.x/pr.z^2  <=>  xr*z^2 == pr.x
    if pr.eq_x_var(&xr) {
        return true;
    }
    // Case 2: xr + n == pr.x/pr.z^2  (only if xr + n < p)
    if FieldElement::cmp_var(&xr, ecdsa_const_p_minus_order()) >= 0 {
        return false;
    }
    let order = ecdsa_const_order_as_fe();
    let mut xr_plus_n = FieldElement::zero();
    xr_plus_n.add(&xr, order);
    pr.eq_x_var(&xr_plus_n)
}

/// Exhaustive verification: convert pr to affine, get x, reduce mod n, compare.
/// Used for debugging; matches libsecp256k1 EXHAUSTIVE_TEST_ORDER path.
pub fn ecdsa_sig_verify_exhaustive(
    sigr: &Scalar,
    sigs: &Scalar,
    pubkey: &Ge,
    message: &Scalar,
) -> bool {
    if sigr.is_zero() || sigs.is_zero() {
        return false;
    }
    let mut sn = Scalar::zero();
    sn.inv_var(sigs);
    let mut u1 = Scalar::zero();
    u1.mul(&sn, message);
    let mut u2 = Scalar::zero();
    u2.mul(&sn, sigr);
    let mut pubkeyj = Gej::default();
    pubkeyj.set_ge(pubkey);
    let one = {
        let mut f = FieldElement::zero();
        f.set_int(1);
        f
    };
    if !pubkeyj.infinity && FieldElement::fe_equal(&pubkeyj.z, &one) {
        let g = generator_g();
        let mut tmp = Gej::default();
        tmp.add_ge_var(&pubkeyj, &g);
        let mut minus_g = Ge::default();
        minus_g.neg(&g);
        let mut minus_gj = Gej::default();
        minus_gj.set_ge(&minus_g);
        pubkeyj.add_var(&tmp, &minus_gj);
    }
    let mut pr = Gej::default();
    ecmult::ecmult(&mut pr, &pubkeyj, &u2, Some(&u1));
    if pr.infinity {
        return false;
    }
    let mut pr_ge = Ge::default();
    pr_ge.set_gej_var(&pr);
    pr_ge.x.normalize();
    let mut c = [0u8; 32];
    pr_ge.x.get_b32(&mut c);
    let mut computed_r = Scalar::zero();
    computed_r.set_b32(&c);
    let ok = bool::from(sigr.ct_eq(&computed_r));
    if !ok {
        let mut sr = [0u8; 32];
        let mut cr = [0u8; 32];
        sigr.get_b32(&mut sr);
        computed_r.get_b32(&mut cr);
        eprintln!("ecdsa_sig_verify_exhaustive: sigr={sr:02x?} computed_r={cr:02x?}");
    }
    ok
}

/// Parse compact 64-byte signature (r || s, big-endian).
/// Returns (sigr, sigs) or None if invalid.
pub fn ecdsa_sig_parse_compact(compact: &[u8; 64]) -> Option<(Scalar, Scalar)> {
    let mut sigr = Scalar::zero();
    let mut sigs = Scalar::zero();
    let mut r_bytes = [0u8; 32];
    let mut s_bytes = [0u8; 32];
    r_bytes.copy_from_slice(&compact[0..32]);
    s_bytes.copy_from_slice(&compact[32..64]);
    if sigr.set_b32(&r_bytes) || sigs.set_b32(&s_bytes) {
        return None; // overflow: r or s >= n
    }
    if sigr.is_zero() || sigs.is_zero() {
        return None;
    }
    Some((sigr, sigs))
}

/// Read DER length. Returns (consumed, length) or None on error.
fn der_read_len(sig: &[u8]) -> Option<(usize, usize)> {
    if sig.is_empty() {
        return None;
    }
    let b1 = sig[0];
    if b1 == 0xFF {
        return None; // X.690: 0xFF shall not be used
    }
    if (b1 & 0x80) == 0 {
        return Some((1, b1 as usize)); // short form
    }
    if b1 == 0x80 {
        return None; // indefinite length not allowed in DER
    }
    let lenleft = (b1 & 0x7F) as usize;
    if lenleft > sig.len() - 1 {
        return None;
    }
    if sig[1] == 0 {
        return None; // not shortest encoding
    }
    if lenleft > 8 {
        return None; // would exceed size_t
    }
    let mut len = 0usize;
    for i in 0..lenleft {
        len = (len << 8) | sig[1 + i] as usize;
    }
    if len < 128 {
        return None; // not shortest encoding
    }
    if len > sig.len() - 1 - lenleft {
        return None;
    }
    Some((1 + lenleft, len))
}

/// Parse DER integer into scalar. Returns (consumed, scalar) or None. Rejects negative.
fn der_parse_integer_strict(sig: &[u8]) -> Option<(usize, Scalar)> {
    if sig.is_empty() || sig[0] != 0x02 {
        return None;
    }
    let (len_consumed, rlen) = der_read_len(&sig[1..])?;
    let content_start = 1 + len_consumed;
    if rlen == 0 || content_start + rlen > sig.len() {
        return None;
    }
    let content = &sig[content_start..content_start + rlen];
    if content[0] == 0x00 && rlen > 1 && (content[1] & 0x80) == 0 {
        return None; // excessive 0x00 padding
    }
    if content[0] == 0xFF && rlen > 1 && (content[1] & 0x80) != 0 {
        return None; // excessive 0xFF padding
    }
    if (content[0] & 0x80) != 0 {
        return None; // negative
    }
    let mut skip = 0;
    if rlen > 0 && content[0] == 0 {
        skip = 1; // skip leading zero
    }
    let actual_len = rlen - skip;
    if actual_len > 32 {
        return None;
    }
    let mut ra = [0u8; 32];
    ra[32 - actual_len..].copy_from_slice(&content[skip..]);
    let mut r = Scalar::zero();
    if r.set_b32(&ra) {
        return None; // overflow
    }
    Some((content_start + rlen, r))
}

/// Parse ECDSA signature from strict DER. BIP66-compliant.
pub fn ecdsa_sig_parse_der(sig: &[u8]) -> Option<(Scalar, Scalar)> {
    if sig.is_empty() || sig[0] != 0x30 {
        return None;
    }
    let (len_consumed, rlen) = der_read_len(&sig[1..])?;
    let content_start = 1 + len_consumed;
    if rlen != sig.len() - content_start {
        return None; // garbage after tuple
    }
    let mut pos = content_start;
    let (r_consumed, sigr) = der_parse_integer_strict(&sig[pos..])?;
    pos += r_consumed;
    let (s_consumed, sigs) = der_parse_integer_strict(&sig[pos..])?;
    pos += s_consumed;
    if pos != sig.len() {
        return None;
    }
    if sigr.is_zero() || sigs.is_zero() {
        return None;
    }
    Some((sigr, sigs))
}

/// Parse ECDSA signature from lax DER. Handles pre-BIP66 non-standard encodings.
/// Matches libsecp256k1 contrib/lax_der_parsing.c (returns None on overflow; C returns invalid sig).
pub fn ecdsa_sig_parse_der_lax(sig: &[u8]) -> Option<(Scalar, Scalar)> {
    if sig.is_empty() || sig[0] != 0x30 {
        return None;
    }
    let mut pos = 1;
    if pos >= sig.len() {
        return None;
    }
    let mut lenbyte = sig[pos];
    pos += 1;
    if (lenbyte & 0x80) != 0 {
        let lenlen = (lenbyte & 0x7F) as usize;
        if lenlen > sig.len() - pos {
            return None;
        }
        pos += lenlen; // skip sequence length bytes
    }
    if pos >= sig.len() || sig[pos] != 0x02 {
        return None;
    }
    pos += 1;
    if pos >= sig.len() {
        return None;
    }
    lenbyte = sig[pos];
    pos += 1;
    let (rlen, rpos) = if (lenbyte & 0x80) != 0 {
        let mut lenlen = (lenbyte & 0x7F) as usize;
        if lenlen > sig.len() - pos {
            return None;
        }
        while lenlen > 0 && pos < sig.len() && sig[pos] == 0 {
            pos += 1;
            lenlen -= 1;
        }
        if lenlen >= 8 {
            return None;
        }
        let mut rlen = 0usize;
        while lenlen > 0 {
            if pos >= sig.len() {
                return None;
            }
            rlen = (rlen << 8) | sig[pos] as usize;
            pos += 1;
            lenlen -= 1;
        }
        if rlen > sig.len() - pos {
            return None;
        }
        let rpos = pos;
        pos += rlen;
        (rlen, rpos)
    } else {
        let rlen = lenbyte as usize;
        if rlen > sig.len() - pos {
            return None;
        }
        let rpos = pos;
        pos += rlen;
        (rlen, rpos)
    };
    if pos >= sig.len() || sig[pos] != 0x02 {
        return None;
    }
    pos += 1;
    if pos >= sig.len() {
        return None;
    }
    lenbyte = sig[pos];
    pos += 1;
    let (slen, spos) = if (lenbyte & 0x80) != 0 {
        let mut lenlen = (lenbyte & 0x7F) as usize;
        if lenlen > sig.len() - pos {
            return None;
        }
        while lenlen > 0 && pos < sig.len() && sig[pos] == 0 {
            pos += 1;
            lenlen -= 1;
        }
        if lenlen >= 8 {
            return None;
        }
        let mut slen = 0usize;
        while lenlen > 0 {
            if pos >= sig.len() {
                return None;
            }
            slen = (slen << 8) | sig[pos] as usize;
            pos += 1;
            lenlen -= 1;
        }
        if slen > sig.len() - pos {
            return None;
        }
        let spos = pos;
        (slen, spos)
    } else {
        let slen = lenbyte as usize;
        if slen > sig.len() - pos {
            return None;
        }
        let spos = pos;
        (slen, spos)
    };
    let mut rlen_actual = rlen;
    let mut rpos_actual = rpos;
    while rlen_actual > 0 && sig[rpos_actual] == 0 {
        rlen_actual -= 1;
        rpos_actual += 1;
    }
    let mut slen_actual = slen;
    let mut spos_actual = spos;
    while slen_actual > 0 && sig[spos_actual] == 0 {
        slen_actual -= 1;
        spos_actual += 1;
    }
    let mut tmpsig = [0u8; 64];
    if rlen_actual > 32 || slen_actual > 32 {
        return None;
    }
    if rlen_actual > 0 {
        tmpsig[32 - rlen_actual..32].copy_from_slice(&sig[rpos_actual..rpos_actual + rlen_actual]);
    }
    if slen_actual > 0 {
        tmpsig[64 - slen_actual..64].copy_from_slice(&sig[spos_actual..spos_actual + slen_actual]);
    }
    ecdsa_sig_parse_compact(&tmpsig)
}

/// Serialize (sigr, sigs) to DER. Minimal encoding (matches libsecp256k1).
pub fn ecdsa_sig_serialize_der(sigr: &Scalar, sigs: &Scalar) -> Vec<u8> {
    let mut r = [0u8; 33];
    let mut s = [0u8; 33];
    let mut r32 = [0u8; 32];
    let mut s32 = [0u8; 32];
    sigr.get_b32(&mut r32);
    sigs.get_b32(&mut s32);
    r[1..33].copy_from_slice(&r32);
    s[1..33].copy_from_slice(&s32);
    let mut rp = 0usize;
    let mut sp = 0usize;
    let mut len_r = 33usize;
    let mut len_s = 33usize;
    while len_r > 1 && r[rp] == 0 && r[rp + 1] < 0x80 {
        len_r -= 1;
        rp += 1;
    }
    while len_s > 1 && s[sp] == 0 && s[sp + 1] < 0x80 {
        len_s -= 1;
        sp += 1;
    }
    let total = 6 + len_r + len_s;
    let mut out = vec![0u8; total];
    out[0] = 0x30;
    out[1] = (4 + len_r + len_s) as u8;
    out[2] = 0x02;
    out[3] = len_r as u8;
    out[4..4 + len_r].copy_from_slice(&r[rp..rp + len_r]);
    out[4 + len_r] = 0x02;
    out[5 + len_r] = len_s as u8;
    out[6 + len_r..].copy_from_slice(&s[sp..sp + len_s]);
    out
}

/// Normalize signature to low-S. If s is high (>= n/2), replace with n - s.
pub fn ecdsa_sig_normalize(sigs: &mut Scalar) {
    if sigs.is_high() {
        let tmp = *sigs;
        sigs.negate(&tmp);
    }
}

/// Serialize Ge to 33-byte compressed pubkey (0x02/0x03 + x).
pub fn ge_to_compressed(ge: &Ge) -> [u8; 33] {
    let mut out = [0u8; 33];
    if ge.infinity {
        return out;
    }
    let mut x = ge.x;
    x.normalize();
    let mut y = ge.y;
    y.normalize();
    out[0] = if y.is_odd() { 0x03 } else { 0x02 };
    let mut x_bytes = [0u8; 32];
    x.get_b32(&mut x_bytes);
    out[1..33].copy_from_slice(&x_bytes);
    out
}

/// Parse 33-byte compressed pubkey to Ge. Returns None if invalid.
pub fn ge_from_compressed(bytes: &[u8; 33]) -> Option<Ge> {
    let prefix = bytes[0];
    if prefix != 0x02 && prefix != 0x03 {
        return None;
    }
    let mut x_bytes = [0u8; 32];
    x_bytes.copy_from_slice(&bytes[1..33]);
    let mut x = FieldElement::zero();
    if !x.set_b32_limit(&x_bytes) {
        return None;
    }
    let mut ge = Ge::default();
    if !ge.set_xo_var(&x, prefix == 0x03) {
        return None;
    }
    Some(ge)
}

/// Parse 65-byte uncompressed pubkey (0x04 || x || y) to Ge. Returns None if invalid.
pub fn ge_from_uncompressed(bytes: &[u8; 65]) -> Option<Ge> {
    if bytes[0] != 0x04 {
        return None;
    }
    let mut x_bytes = [0u8; 32];
    let mut y_bytes = [0u8; 32];
    x_bytes.copy_from_slice(&bytes[1..33]);
    y_bytes.copy_from_slice(&bytes[33..65]);
    let mut x = FieldElement::zero();
    let mut y = FieldElement::zero();
    if !x.set_b32_limit(&x_bytes) || !y.set_b32_limit(&y_bytes) {
        return None;
    }
    let mut ge = Ge {
        x,
        y,
        infinity: false,
    };
    ge.y.normalize();
    // Verify point is on curve: y^2 == x^3 + 7
    let mut x3 = FieldElement::zero();
    let mut x2 = FieldElement::zero();
    x2.sqr(&ge.x);
    x3.mul(&ge.x, &x2);
    x3.add_int(7); // secp256k1 curve constant B
    let mut y2 = FieldElement::zero();
    y2.sqr(&ge.y);
    if !FieldElement::fe_equal(&y2, &x3) {
        return None;
    }
    Some(ge)
}

/// Parse pubkey (33-byte compressed or 65-byte uncompressed) to Ge. Returns None if invalid.
pub fn ge_from_pubkey_bytes(bytes: &[u8]) -> Option<Ge> {
    if bytes.len() == 33 {
        let arr: [u8; 33] = bytes.try_into().unwrap();
        ge_from_compressed(&arr)
    } else if bytes.len() == 65 {
        let arr: [u8; 65] = bytes.try_into().unwrap();
        ge_from_uncompressed(&arr)
    } else {
        None
    }
}

/// Derive public key from secret. Returns affine point G * secret.
pub fn pubkey_from_secret(secret: &Scalar) -> Ge {
    let mut rj = Gej::default();
    ecmult::ecmult_gen(&mut rj, secret);
    let mut r = Ge::default();
    r.set_gej_var(&rj);
    r
}

/// Recover public key from ECDSA signature (r,s) with recovery id.
/// recid in 0..4: (recid & 2) indicates xr overflow (use xr+n), (recid & 1) indicates y parity.
pub fn ecdsa_sig_recover(sigr: &Scalar, sigs: &Scalar, message: &Scalar, recid: u8) -> Option<Ge> {
    if sigr.is_zero() || sigs.is_zero() || recid > 3 {
        return None;
    }
    let mut brx = [0u8; 32];
    sigr.get_b32(&mut brx);
    let mut fx = FieldElement::zero();
    if !fx.set_b32_limit(&brx) {
        return None;
    }
    if (recid & 2) != 0 {
        if FieldElement::cmp_var(&fx, ecdsa_const_p_minus_order()) >= 0 {
            return None;
        }
        let order = ecdsa_const_order_as_fe();
        let mut tmp = FieldElement::zero();
        tmp.add(&fx, order);
        fx = tmp;
    }
    let mut r_ge = Ge::default();
    if !r_ge.set_xo_var(&fx, (recid & 1) != 0) {
        return None;
    }
    let mut rj = Gej::default();
    rj.set_ge(&r_ge);
    let mut rn = Scalar::zero();
    rn.inv_var(sigr);
    let mut u1 = Scalar::zero();
    u1.mul(&rn, message);
    let u1_val = u1;
    u1.negate(&u1_val);
    let mut u2 = Scalar::zero();
    u2.mul(&rn, sigs);
    let mut qj = Gej::default();
    ecmult::ecmult(&mut qj, &rj, &u2, None);
    let mut u1_g = Gej::default();
    ecmult::ecmult_gen(&mut u1_g, &u1);
    let qj_in = qj;
    qj.add_var(&qj_in, &u1_g);
    if qj.infinity {
        return None;
    }
    let mut pubkey = Ge::default();
    pubkey.set_gej_var(&qj);
    Some(pubkey)
}

/// ECDSA sign with recovery id. Returns (sigr, sigs, recid) or None.
/// recid encodes (x_overflow << 1) | y_parity for pubkey recovery.
pub fn ecdsa_sig_sign_recoverable(
    seckey: &Scalar,
    message: &Scalar,
    nonce: &Scalar,
) -> Option<(Scalar, Scalar, u8)> {
    let mut rp = Gej::default();
    ecmult::ecmult_gen(&mut rp, nonce);
    let mut r_ge = Ge::default();
    r_ge.set_gej_var(&rp);
    r_ge.x.normalize();
    r_ge.y.normalize();

    let mut b = [0u8; 32];
    r_ge.x.get_b32(&mut b);
    let mut sigr = Scalar::zero();
    sigr.set_b32(&b);
    if sigr.is_zero() {
        return None;
    }

    let mut n = Scalar::zero();
    n.mul(&sigr, seckey);
    let n_mul = n;
    n.add(&n_mul, message);

    let mut sigs = Scalar::zero();
    sigs.inv_var(nonce);
    let inv = sigs;
    sigs.mul(&inv, &n);

    let high = sigs.is_high();
    sigs.cond_negate(high as i32);

    if sigs.is_zero() {
        return None;
    }

    let order = ecdsa_const_order_as_fe();
    let overflow = FieldElement::cmp_var(&r_ge.x, order) >= 0;
    let recid = (overflow as u8) << 1 | (r_ge.y.is_odd() as u8);
    Some((sigr, sigs, recid))
}

/// ECDSA sign. Returns (sigr, sigs) or None on failure.
/// Uses provided nonce (RFC 6979 or random). Enforces low-S.
pub fn ecdsa_sig_sign(
    seckey: &Scalar,
    message: &Scalar,
    nonce: &Scalar,
) -> Option<(Scalar, Scalar)> {
    let mut rp = Gej::default();
    ecmult::ecmult_gen(&mut rp, nonce);
    let mut r_ge = Ge::default();
    r_ge.set_gej_var(&rp);
    r_ge.x.normalize();
    r_ge.y.normalize();

    let mut b = [0u8; 32];
    r_ge.x.get_b32(&mut b);
    let mut sigr = Scalar::zero();
    sigr.set_b32(&b);
    if sigr.is_zero() {
        return None;
    }

    let mut n = Scalar::zero();
    n.mul(&sigr, seckey);
    let n_mul = n;
    n.add(&n_mul, message);

    let mut sigs = Scalar::zero();
    sigs.inv_var(nonce);
    let inv = sigs;
    sigs.mul(&inv, &n);

    let high = sigs.is_high();
    sigs.cond_negate(high as i32);

    if sigs.is_zero() {
        return None;
    }
    Some((sigr, sigs))
}

/// Verify ECDSA signature directly from DER bytes, pubkey bytes, and message hash.
/// Bypasses all intermediate serialization/deserialization.
/// `der_sig`: DER-encoded signature (WITHOUT sighash byte)
/// `pubkey_bytes`: 33-byte compressed or 65-byte uncompressed pubkey
/// `msg_hash`: 32-byte message hash
/// `strict_der`: if true, use strict DER parsing (BIP66); if false, use lax DER
/// `enforce_low_s`: if true, reject high-S signatures
/// Returns: Some(true) = valid, Some(false) = invalid sig, None = parse error
#[inline]
pub fn verify_ecdsa_direct(
    der_sig: &[u8],
    pubkey_bytes: &[u8],
    msg_hash: &[u8; 32],
    strict_der: bool,
    enforce_low_s: bool,
) -> Option<bool> {
    let (sigr, mut sigs) = if strict_der {
        ecdsa_sig_parse_der(der_sig)?
    } else {
        ecdsa_sig_parse_der_lax(der_sig)?
    };
    if enforce_low_s && sigs.is_high() {
        return Some(false);
    }
    if sigs.is_high() {
        let tmp = sigs;
        sigs.negate(&tmp);
    }
    let pk = ge_from_pubkey_bytes(pubkey_bytes)?;
    let mut msg = Scalar::zero();
    if msg.set_b32(msg_hash) {
        return None;
    }
    Some(ecdsa_sig_verify(&sigr, &sigs, &pk, &msg))
}

#[cfg(test)]
mod bip66_strict_der_regression {
    use super::*;

    /// Strict DER (BIP66) parse: mainnet-era production signature.
    /// Provenance: P2PKH spend from historic mainnet (height 363_726, tx index 2).
    #[test]
    fn strict_der_parses_bip66_epoch_mainnet_vector() {
        let der_hex = "304502206794ddbfbb9a36d9daecbab5b19ca523d3fa53530a7467475d812b6375204a85\
            0221009cab7076ce3319ab94a31df5c63b1b827b8d99b70822b0e9625a4aff2433a7d0";
        let der = hex::decode(der_hex).expect("hex");
        assert_eq!(der.len(), 71, "DER without sighash byte");
        let parsed = ecdsa_sig_parse_der(&der);
        assert!(
            parsed.is_some(),
            "strict DER parse failed for known-mainnet sig"
        );
    }
}
