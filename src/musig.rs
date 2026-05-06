//! MuSig2 for BIP340-compatible multi-signatures (BIP 327).
//!
//! Key aggregation, nonce gen/agg, nonce_process, partial_sign, partial_sig_verify, partial_sig_agg.
//! Cross-validated with secp256k1-fork (libsecp256k1 FFI).
//!
//! **Timing:** `nonce_gen` and tweak application use constant-time `k*G` for secret scalars.
//! Verification and key aggregation (public points) use variable-time `ecmult` / `ecmult_gen`.

use sha2::{Digest, Sha256};

use crate::ecdsa::{ge_from_compressed, ge_to_compressed, pubkey_from_secret};
use crate::ecmult;
use crate::ecmult_gen_const;
use crate::field::FieldElement;
use crate::group::{generator_g, Ge, Gej};
use crate::scalar::Scalar;
use crate::schnorr;

fn tagged_hash(tag: &[u8], data: &[u8]) -> [u8; 32] {
    let tag_hash = Sha256::digest(tag);
    let mut hasher = Sha256::new();
    hasher.update(tag_hash);
    hasher.update(tag_hash);
    hasher.update(data);
    hasher.finalize().into()
}

const TAG_KEYAGG_LIST: &[u8] = b"KeyAgg list";
const TAG_KEYAGG_COEFF: &[u8] = b"KeyAgg coefficient";
const TAG_MUSIG_AUX: &[u8] = b"MuSig/aux";
const TAG_MUSIG_NONCE: &[u8] = b"MuSig/nonce";
const TAG_MUSIG_NONCECOEF: &[u8] = b"MuSig/noncecoef";

/// Output of nonce_gen: ((ell64, pubkey33), msg66).
pub type NonceGenOutput = (([u8; 64], [u8; 33]), [u8; 66]);

/// Key aggregation cache. Holds Q (agg pubkey), second_pk, pks_hash, parity_acc, tweak.
#[derive(Clone)]
pub struct KeyAggCache {
    pub pk: Ge,
    pub second_pk: Ge,
    pub pks_hash: [u8; 32],
    pub parity_acc: u8,
    pub tweak: Scalar,
}

impl KeyAggCache {
    /// Aggregate public keys. Returns None if any pubkey invalid or Q is infinity.
    pub fn new(pubkeys: &[[u8; 33]]) -> Option<Self> {
        if pubkeys.is_empty() {
            return None;
        }
        let pks_hash = hash_keys(pubkeys)?;
        let second_pk = get_second_key(pubkeys)?;

        let mut qj = Gej::default();
        qj.set_infinity();
        for pk_bytes in pubkeys {
            let p = ge_from_compressed(pk_bytes)?;
            if p.infinity {
                return None;
            }
            let a = keyagg_coeff_internal(&pks_hash, &p, &second_pk);
            let mut pj = Gej::default();
            pj.set_ge(&p);
            let mut ap = Gej::default();
            ecmult::ecmult(&mut ap, &pj, &a, None);
            let qj_old = qj;
            qj.add_var(&qj_old, &ap);
        }

        let mut q = Ge::default();
        q.set_gej_var(&qj);
        if q.infinity {
            return None;
        }

        q.x.normalize();
        q.y.normalize();

        Some(KeyAggCache {
            pk: q,
            second_pk,
            pks_hash,
            parity_acc: 0,
            tweak: Scalar::zero(),
        })
    }

    /// X-only aggregate public key (32 bytes).
    pub fn agg_pk_xonly(&self) -> [u8; 32] {
        let mut q = self.pk;
        if q.y.is_odd() {
            let mut neg = Ge::default();
            neg.neg(&q);
            q = neg;
        }
        let mut out = [0u8; 32];
        q.x.normalize();
        q.x.get_b32(&mut out);
        out
    }

    /// Plain (compressed) aggregate public key.
    pub fn agg_pk_plain(&self) -> [u8; 33] {
        ge_to_compressed(&self.pk)
    }

    /// Apply plain EC tweak (BIP32). Returns None if invalid.
    pub fn pubkey_ec_tweak_add(&mut self, tweak32: &[u8; 32]) -> Option<()> {
        let mut t = Scalar::zero();
        if t.set_b32(tweak32) {
            return None; // overflow
        }
        self.apply_tweak(&t, false)
    }

    /// Apply x-only tweak (Taproot). Returns None if invalid.
    pub fn pubkey_xonly_tweak_add(&mut self, tweak32: &[u8; 32]) -> Option<()> {
        let mut t = Scalar::zero();
        if t.set_b32(tweak32) {
            return None;
        }
        self.apply_tweak(&t, true)
    }

    fn apply_tweak(&mut self, t: &Scalar, is_xonly: bool) -> Option<()> {
        let mut g = 1i32;
        if is_xonly && self.pk.y.is_odd() {
            g = -1;
            self.parity_acc ^= 1;
            let mut neg_t = Scalar::zero();
            neg_t.negate(&self.tweak);
            self.tweak = neg_t;
        }
        let mut t_acc = Scalar::zero();
        t_acc.add(&self.tweak, t);
        self.tweak = t_acc;

        let mut tj = Gej::default();
        ecmult_gen_const(&mut tj, t);
        let mut pj = Gej::default();
        pj.set_ge(&self.pk);
        if g < 0 {
            let mut neg_p = Ge::default();
            neg_p.neg(&self.pk);
            pj.set_ge(&neg_p);
        }
        let pj_old = pj;
        pj.add_var(&pj_old, &tj);
        self.pk.set_gej_var(&pj);
        if self.pk.infinity {
            return None;
        }
        self.pk.x.normalize();
        self.pk.y.normalize();
        Some(())
    }
}

fn hash_keys(pubkeys: &[[u8; 33]]) -> Option<[u8; 32]> {
    let mut data = Vec::with_capacity(pubkeys.len() * 33);
    for pk in pubkeys {
        data.extend_from_slice(pk);
    }
    Some(tagged_hash(TAG_KEYAGG_LIST, &data))
}

fn get_second_key(pubkeys: &[[u8; 33]]) -> Option<Ge> {
    let first = &pubkeys[0];
    for pk in pubkeys.iter().skip(1) {
        if pk != first {
            return ge_from_compressed(pk);
        }
    }
    let mut inf = Ge::default();
    inf.set_infinity();
    Some(inf)
}

fn keyagg_coeff_internal(pks_hash: &[u8; 32], pk: &Ge, second_pk: &Ge) -> Scalar {
    if !second_pk.infinity && ge_eq(pk, second_pk) {
        let mut one = Scalar::zero();
        one.set_int(1);
        return one;
    }
    let mut data = Vec::with_capacity(32 + 33);
    data.extend_from_slice(pks_hash);
    data.extend_from_slice(&ge_to_compressed(pk));
    let h = tagged_hash(TAG_KEYAGG_COEFF, &data);
    let mut r = Scalar::zero();
    r.set_b32(&h);
    r
}

fn ge_eq(a: &Ge, b: &Ge) -> bool {
    if a.infinity != b.infinity {
        return false;
    }
    if a.infinity {
        return true;
    }
    let mut ax = a.x;
    let mut bx = b.x;
    ax.normalize();
    bx.normalize();
    let mut ay = a.y;
    let mut by = b.y;
    ay.normalize();
    by.normalize();
    FieldElement::fe_equal(&ax, &bx) && FieldElement::fe_equal(&ay, &by)
}

/// Parse 66-byte pubnonce (33+33 compressed R1, R2). Returns None if invalid.
pub fn pubnonce_parse(in66: &[u8; 66]) -> Option<[u8; 66]> {
    let r1_bytes: [u8; 33] = in66[0..33].try_into().unwrap();
    let r2_bytes: [u8; 33] = in66[33..66].try_into().unwrap();
    let _r1 = ge_from_compressed(&r1_bytes)?;
    let _r2 = ge_from_compressed(&r2_bytes)?;
    if _r1.infinity || _r2.infinity {
        return None;
    }
    Some(*in66)
}

/// Parse 66-byte aggnonce. 33 zero bytes = infinity for that component.
pub fn aggnonce_parse(in66: &[u8; 66]) -> Option<(Ge, Ge)> {
    let r1_bytes: [u8; 33] = in66[0..33].try_into().unwrap();
    let r2_bytes: [u8; 33] = in66[33..66].try_into().unwrap();
    let r1 = ge_from_compressed_ext(&r1_bytes)?;
    let r2 = ge_from_compressed_ext(&r2_bytes)?;
    Some((r1, r2))
}

fn ge_from_compressed_ext(bytes: &[u8; 33]) -> Option<Ge> {
    if bytes.iter().all(|&b| b == 0) {
        let mut inf = Ge::default();
        inf.set_infinity();
        return Some(inf);
    }
    ge_from_compressed(bytes)
}

fn ge_to_compressed_ext(ge: &Ge) -> [u8; 33] {
    if ge.infinity {
        return [0u8; 33];
    }
    ge_to_compressed(ge)
}

/// Aggregate pubnonces: R1 = sum R1_i, R2 = sum R2_i.
pub fn nonce_agg(pubnonces: &[[u8; 66]]) -> Option<[u8; 66]> {
    if pubnonces.is_empty() {
        return None;
    }
    let mut r1j = Gej::default();
    let mut r2j = Gej::default();
    r1j.set_infinity();
    r2j.set_infinity();

    for pn in pubnonces {
        let (r1, r2) = aggnonce_parse(pn)?;
        let mut r1_gej = Gej::default();
        r1_gej.set_ge(&r1);
        let mut r2_gej = Gej::default();
        r2_gej.set_ge(&r2);
        let r1j_old = r1j;
        r1j.add_var(&r1j_old, &r1_gej);
        let r2j_old = r2j;
        r2j.add_var(&r2j_old, &r2_gej);
    }

    let mut r1_ge = Ge::default();
    let mut r2_ge = Ge::default();
    r1_ge.set_gej_var(&r1j);
    r2_ge.set_gej_var(&r2j);

    let mut out = [0u8; 66];
    out[0..33].copy_from_slice(&ge_to_compressed_ext(&r1_ge));
    out[33..66].copy_from_slice(&ge_to_compressed_ext(&r2_ge));
    Some(out)
}

/// Session context for signing. Created from aggnonce + msg + keyagg_cache.
#[derive(Clone)]
pub struct Session {
    pub fin_nonce_parity: u8,
    pub fin_nonce: [u8; 32],
    pub noncecoef: Scalar,
    pub challenge: Scalar,
    pub s_part: Scalar,
}

/// Create session from aggregated nonce, 32-byte message, and keyagg cache.
pub fn nonce_process(
    aggnonce: &[u8; 66],
    msg32: &[u8; 32],
    keyagg_cache: &KeyAggCache,
) -> Option<Session> {
    let (r1, r2) = aggnonce_parse(aggnonce)?;

    let mut agg_pk32 = [0u8; 32];
    let mut pk_x = keyagg_cache.pk.x;
    pk_x.normalize();
    pk_x.get_b32(&mut agg_pk32);

    let mut nonce_data = Vec::with_capacity(33 + 33 + 32 + 32);
    nonce_data.extend_from_slice(&ge_to_compressed_ext(&r1));
    nonce_data.extend_from_slice(&ge_to_compressed_ext(&r2));
    nonce_data.extend_from_slice(&agg_pk32);
    nonce_data.extend_from_slice(msg32);
    let noncehash = tagged_hash(TAG_MUSIG_NONCECOEF, &nonce_data);

    let mut b = Scalar::zero();
    b.set_b32(&noncehash);

    let mut fin_nonce_ptj = Gej::default();
    fin_nonce_ptj.set_ge(&r1);
    let mut r2j = Gej::default();
    r2j.set_ge(&r2);
    let mut b_r2 = Gej::default();
    ecmult::ecmult(&mut b_r2, &r2j, &b, None);
    let fin_old = fin_nonce_ptj;
    fin_nonce_ptj.add_var(&fin_old, &b_r2);

    let mut fin_nonce_pt = Ge::default();
    fin_nonce_pt.set_gej_var(&fin_nonce_ptj);

    if fin_nonce_pt.infinity {
        fin_nonce_pt = generator_g();
    }

    fin_nonce_pt.x.normalize();
    fin_nonce_pt.y.normalize();
    let fin_nonce_parity = if fin_nonce_pt.y.is_odd() { 1 } else { 0 };
    let mut fin_nonce = [0u8; 32];
    fin_nonce_pt.x.get_b32(&mut fin_nonce);

    let mut challenge_data = Vec::with_capacity(32 + 32 + 32);
    challenge_data.extend_from_slice(&fin_nonce);
    challenge_data.extend_from_slice(&agg_pk32);
    challenge_data.extend_from_slice(msg32);
    let e_hash = schnorr::tagged_hash_challenge(&challenge_data);
    let mut challenge = Scalar::zero();
    challenge.set_b32(&e_hash);

    let mut s_part = Scalar::zero();
    if !keyagg_cache.tweak.is_zero() {
        let mut e_tweak = Scalar::zero();
        e_tweak.mul(&challenge, &keyagg_cache.tweak);
        if keyagg_cache.pk.y.is_odd() {
            let e_val = e_tweak;
            e_tweak.negate(&e_val);
        }
        s_part = e_tweak;
    }

    Some(Session {
        fin_nonce_parity,
        fin_nonce,
        noncecoef: b,
        challenge,
        s_part,
    })
}

/// BIP 327 nonce input: bytes(1,len) = 1-byte length then data. bytes(4,len) = 3 zeros + 1-byte len.
fn nonce_helper_1(data: Option<&[u8]>) -> Vec<u8> {
    let mut out = Vec::new();
    if let Some(d) = data {
        out.push(d.len() as u8);
        out.extend_from_slice(d);
    } else {
        out.push(0);
    }
    out
}

fn nonce_helper_8(data: Option<&[u8]>) -> Vec<u8> {
    let mut out = vec![0u8; 7];
    if let Some(d) = data {
        out.push(d.len() as u8);
        out.extend_from_slice(d);
    } else {
        out.push(0);
    }
    out
}

fn nonce_helper_4(data: Option<&[u8]>) -> Vec<u8> {
    let mut out = vec![0u8; 3];
    if let Some(d) = data {
        out.push(d.len() as u8);
        out.extend_from_slice(d);
    } else {
        out.push(0);
    }
    out
}

/// Generate nonce pair. session_secrand32 must be unique per call; zeroed on success.
/// pubkey33: signer's compressed pubkey (required).
/// Returns (secnonce_k1_k2_pk, pubnonce) or None. secnonce is (k1[32], k2[32], pk[33]) = 97 bytes.
pub fn nonce_gen(
    session_secrand32: &mut [u8; 32],
    seckey: Option<&[u8; 32]>,
    pubkey33: &[u8; 33],
    msg32: Option<&[u8; 32]>,
    keyagg_cache: Option<&KeyAggCache>,
    extra_input32: Option<&[u8; 32]>,
) -> Option<NonceGenOutput> {
    if session_secrand32.iter().all(|&b| b == 0) {
        return None;
    }
    if let Some(sk) = seckey {
        let mut s = Scalar::zero();
        if s.set_b32(sk) || s.is_zero() {
            return None;
        }
    }

    let mut rand = [0u8; 32];
    if let Some(_sk) = seckey {
        let aux_hash = tagged_hash(TAG_MUSIG_AUX, session_secrand32);
        for i in 0..32 {
            rand[i] = session_secrand32[i] ^ aux_hash[i];
        }
    } else {
        rand.copy_from_slice(session_secrand32);
    }

    let aggpk = keyagg_cache.map(|c| c.agg_pk_xonly());
    let aggpk_slice: Option<&[u8]> = aggpk.as_ref().map(|a| a.as_slice());

    let mut data = Vec::new();
    data.extend_from_slice(&rand);
    data.extend(nonce_helper_1(Some(pubkey33)));
    data.extend(nonce_helper_1(aggpk_slice));
    data.push(if msg32.is_some() { 1 } else { 0 });
    if let Some(m) = msg32 {
        data.extend(nonce_helper_8(Some(m)));
    } else {
        data.extend(nonce_helper_8(None));
    }
    data.extend(nonce_helper_4(extra_input32.map(|e| e.as_slice())));

    let mut k1 = Scalar::zero();
    let mut k2 = Scalar::zero();
    {
        let mut d = data.clone();
        d.push(0);
        let h1 = tagged_hash(TAG_MUSIG_NONCE, &d);
        k1.set_b32(&h1);
        if k1.is_zero() {
            return None;
        }
    }
    {
        let mut d = data;
        d.push(1);
        let h2 = tagged_hash(TAG_MUSIG_NONCE, &d);
        k2.set_b32(&h2);
        if k2.is_zero() {
            return None;
        }
    }

    let mut r1j = Gej::default();
    let mut r2j = Gej::default();
    ecmult_gen_const(&mut r1j, &k1);
    ecmult_gen_const(&mut r2j, &k2);

    let mut r1 = Ge::default();
    let mut r2 = Ge::default();
    r1.set_gej_var(&r1j);
    r2.set_gej_var(&r2j);

    let mut pubnonce = [0u8; 66];
    pubnonce[0..33].copy_from_slice(&ge_to_compressed(&r1));
    pubnonce[33..66].copy_from_slice(&ge_to_compressed(&r2));

    let mut secnonce_k = [0u8; 64];
    k1.get_b32((&mut secnonce_k[0..32]).try_into().unwrap());
    k2.get_b32((&mut secnonce_k[32..64]).try_into().unwrap());

    session_secrand32.copy_from_slice(&[0u8; 32]);

    Some(((secnonce_k, *pubkey33), pubnonce))
}

/// Produce partial signature. Consumes secnonce (caller must not reuse).
/// seckey must match pubkey33. Returns 32-byte partial sig or None.
pub fn partial_sign(
    secnonce: &mut ([u8; 64], [u8; 33]),
    seckey: &[u8; 32],
    keyagg_cache: &KeyAggCache,
    session: &Session,
) -> Option<[u8; 32]> {
    if secnonce.0.iter().all(|&b| b == 0) {
        return None; // already used
    }
    let mut k1 = Scalar::zero();
    let mut k2 = Scalar::zero();
    k1.set_b32(secnonce.0[0..32].try_into().unwrap());
    k2.set_b32(secnonce.0[32..64].try_into().unwrap());
    if k1.is_zero() || k2.is_zero() {
        return None;
    }
    let pk_from_nonce = ge_from_compressed(&secnonce.1)?;
    let pk_from_sk = pubkey_from_secret(&{
        let mut d = Scalar::zero();
        d.set_b32(seckey);
        d
    });
    if !ge_eq(&pk_from_nonce, &pk_from_sk) {
        return None;
    }

    let mut sk = Scalar::zero();
    if sk.set_b32(seckey) || sk.is_zero() {
        return None;
    }

    // Negate sk if Q has odd y XOR parity_acc (BIP 327: d = g*gacc*d').
    // Condition is on PUBLIC aggregate key parity — not the secret key itself.
    // Use cond_negate for consistency with the library's CT style.
    let negate_sk =
        ((keyagg_cache.pk.y.is_odd() as i32) ^ (keyagg_cache.parity_acc as i32 & 1)) & 1;
    sk.cond_negate(negate_sk);

    let mu = keyagg_coeff_internal(
        &keyagg_cache.pks_hash,
        &pk_from_nonce,
        &keyagg_cache.second_pk,
    );
    let mut sk_mu = Scalar::zero();
    sk_mu.mul(&sk, &mu);

    let mut k1_adj = k1;
    let mut k2_adj = k2;
    // Negate both nonce halves if finalized nonce R has odd y (BIP 327).
    // fin_nonce_parity is a PUBLIC value (derived from public nonces),
    // but k1/k2 are SECRET — use cond_negate to avoid an observable branch.
    let parity = session.fin_nonce_parity as i32;
    k1_adj.cond_negate(parity);
    k2_adj.cond_negate(parity);

    let mut s = Scalar::zero();
    s.mul(&session.challenge, &sk_mu);
    let mut k2_b = Scalar::zero();
    k2_b.mul(&session.noncecoef, &k2_adj);
    let mut k_sum = Scalar::zero();
    k_sum.add(&k1_adj, &k2_b);
    let s_old = s;
    s.add(&s_old, &k_sum);

    secnonce.0.copy_from_slice(&[0u8; 64]);

    let mut out = [0u8; 32];
    s.get_b32(&mut out);
    Some(out)
}

/// Verify a partial signature.
pub fn partial_sig_verify(
    partial_sig: &[u8; 32],
    pubnonce: &[u8; 66],
    pubkey33: &[u8; 33],
    keyagg_cache: &KeyAggCache,
    session: &Session,
) -> bool {
    let (r1, r2) = match aggnonce_parse(pubnonce) {
        Some(x) => x,
        None => return false,
    };
    if r1.infinity || r2.infinity {
        return false; // pubnonce uses compressed, not ext
    }
    let pk = match ge_from_compressed(pubkey33) {
        Some(p) => p,
        None => return false,
    };

    let mut s = Scalar::zero();
    if s.set_b32(partial_sig) {
        return false;
    }

    let mut rj = Gej::default();
    rj.set_ge(&r1);
    let mut r2j = Gej::default();
    r2j.set_ge(&r2);
    let mut b_r2 = Gej::default();
    ecmult::ecmult(&mut b_r2, &r2j, &session.noncecoef, None);
    let rj_old = rj;
    rj.add_var(&rj_old, &b_r2);
    let mut rj_final = rj;
    if session.fin_nonce_parity != 0 {
        let rj_old = rj_final;
        rj_final.neg(&rj_old);
    }

    let mu = keyagg_coeff_internal(&keyagg_cache.pks_hash, &pk, &keyagg_cache.second_pk);
    let mut e_mu = Scalar::zero();
    e_mu.mul(&session.challenge, &mu);
    if keyagg_cache.pk.y.is_odd() != (keyagg_cache.parity_acc != 0) {
        let e_val = e_mu;
        e_mu.negate(&e_val);
    }

    let mut neg_s = Scalar::zero();
    neg_s.negate(&s);
    let mut pkj = Gej::default();
    pkj.set_ge(&pk);
    let mut tmp = Gej::default();
    ecmult::ecmult(&mut tmp, &pkj, &e_mu, Some(&neg_s));
    let tmp_old = tmp;
    tmp.add_var(&tmp_old, &rj_final);
    tmp.is_infinity()
}

/// Aggregate partial signatures into final 64-byte BIP 340 signature.
pub fn partial_sig_agg(session: &Session, partial_sigs: &[[u8; 32]]) -> Option<[u8; 64]> {
    if partial_sigs.is_empty() {
        return None;
    }
    let mut s_sum = session.s_part;
    for ps in partial_sigs {
        let mut term = Scalar::zero();
        if term.set_b32(ps) {
            return None;
        }
        let s_old = s_sum;
        s_sum.add(&s_old, &term);
    }
    let mut out = [0u8; 64];
    out[0..32].copy_from_slice(&session.fin_nonce);
    s_sum.get_b32((&mut out[32..64]).try_into().unwrap());
    Some(out)
}
