//! ElligatorSwift encoding for secp256k1 (BIP 324 v2 transport).
//!
//! 64-byte encodings (u || t) that are indistinguishable from random.
//! See https://eprint.iacr.org/2022/759 and libsecp256k1 doc/ellswift.md.
//!
//! Only built for x86_64/aarch64 (requires field is_square_var).

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
mod imp {
    use sha2::{Digest, Sha256};

    use crate::ecdsa::{ge_from_compressed, pubkey_from_secret};
    use crate::ecmult_const;
    use crate::field::FieldElement;
    use crate::group::{Ge, Gej};
    use crate::scalar::Scalar;

    const SECP256K1_B: u32 = 7;

    /// c1 = (sqrt(-3)-1)/2
    fn ellswift_c1() -> FieldElement {
        let b: [u8; 32] = [
            0x85, 0x16, 0x95, 0xd4, 0x9a, 0x83, 0xf8, 0xef, 0x91, 0x9b, 0xb8, 0x61, 0x53, 0xcb,
            0xcb, 0x16, 0x63, 0x0f, 0xb6, 0x8a, 0xed, 0x0a, 0x76, 0x6a, 0x3e, 0xc6, 0x93, 0xd6,
            0x8e, 0x6a, 0xfa, 0x40,
        ];
        let mut r = FieldElement::zero();
        r.set_b32_mod(&b);
        r
    }

    /// c2 = (-sqrt(-3)-1)/2
    fn ellswift_c2() -> FieldElement {
        let b: [u8; 32] = [
            0x7a, 0xe9, 0x6a, 0x2b, 0x65, 0x7c, 0x07, 0x10, 0x6e, 0x64, 0x47, 0x9e, 0xac, 0x34,
            0x34, 0xe9, 0x9c, 0xf0, 0x49, 0x75, 0x12, 0xf5, 0x89, 0x95, 0xc1, 0x39, 0x6c, 0x28,
            0x71, 0x95, 0x01, 0xee,
        ];
        let mut r = FieldElement::zero();
        r.set_b32_mod(&b);
        r
    }

    /// c3 = (-sqrt(-3)+1)/2 = -c1 = c2+1
    fn ellswift_c3() -> FieldElement {
        let b: [u8; 32] = [
            0x7a, 0xe9, 0x6a, 0x2b, 0x65, 0x7c, 0x07, 0x10, 0x6e, 0x64, 0x47, 0x9e, 0xac, 0x34,
            0x34, 0xe9, 0x9c, 0xf0, 0x49, 0x75, 0x12, 0xf5, 0x89, 0x95, 0xc1, 0x39, 0x6c, 0x28,
            0x71, 0x95, 0x01, 0xef,
        ];
        let mut r = FieldElement::zero();
        r.set_b32_mod(&b);
        r
    }

    /// c4 = (sqrt(-3)+1)/2 = -c2 = c1+1
    fn ellswift_c4() -> FieldElement {
        let b: [u8; 32] = [
            0x85, 0x16, 0x95, 0xd4, 0x9a, 0x83, 0xf8, 0xef, 0x91, 0x9b, 0xb8, 0x61, 0x53, 0xcb,
            0xcb, 0x16, 0x63, 0x0f, 0xb6, 0x8a, 0xed, 0x0a, 0x76, 0x6a, 0x3e, 0xc6, 0x93, 0xd6,
            0x8e, 0x6a, 0xfa, 0x41,
        ];
        let mut r = FieldElement::zero();
        r.set_b32_mod(&b);
        r
    }

    /// Inverse of xswiftec: find t such that xswiftec_var(u, t) = x.
    /// c in 0..8 selects branch. Returns true if successful.
    fn xswiftec_inv_var(
        t: &mut FieldElement,
        x_in: &FieldElement,
        u_in: &FieldElement,
        c: i32,
    ) -> bool {
        debug_assert!((0..8).contains(&c));
        let mut x = *x_in;
        let mut u = *u_in;
        x.normalize_weak();
        u.normalize_weak();

        if !Ge::x_on_curve_var(&x) {
            return false;
        }

        #[allow(unused_assignments)]
        let (mut g, mut v, mut s, mut m, mut r, mut q) = (
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
            FieldElement::zero(),
        );

        if (c & 2) == 0 {
            // c in {0, 1, 4, 5}: x1 or x2 formula
            m = x;
            m.add_assign(&u);
            let m_val = m;
            m.negate(&m_val, 2);
            if Ge::x_on_curve_var(&m) {
                return false;
            }
            s = m;
            let m_sq = m;
            s.sqr(&m_sq);
            let s_val = s;
            s.negate(&s_val, 1);
            m.mul(&u, &x);
            s.add_assign(&m);
            if s.normalizes_to_zero_var() {
                return false;
            }
            g = u;
            let u_sq = u;
            g.sqr(&u_sq);
            let g_val = g;
            g.mul(&g_val, &u);
            g.add_int(SECP256K1_B);
            m.mul(&s, &g);
            if !FieldElement::is_square_var(&m) {
                return false;
            }
            let mut s_inv = FieldElement::zero();
            s_inv.inv(&s);
            s.mul(&s_inv, &g);
            v = x;
        } else {
            // c in {2, 3, 6, 7}: x3 formula
            m.negate(&u, 1);
            s = m;
            s.add_assign(&x);
            if !FieldElement::is_square_var(&s) {
                return false;
            }
            g = u;
            let u_sq = u;
            g.sqr(&u_sq);
            q.mul(&s, &g);
            q.mul_int(3);
            let g_val = g;
            g.mul(&g_val, &u);
            g.mul_int(4);
            g.add_int(4 * SECP256K1_B);
            q.add_assign(&g);
            let q_val = q;
            let s_val = s;
            q.mul(&q_val, &s_val);
            let q_val2 = q;
            q.negate(&q_val2, 1);
            if !FieldElement::is_square_var(&q) {
                return false;
            }
            if !r.sqrt(&q) {
                return false;
            }
            if (c & 1) != 0 && r.normalizes_to_zero_var() {
                return false;
            }
            if s.normalizes_to_zero_var() {
                return false;
            }
            let mut s_inv = FieldElement::zero();
            s_inv.inv(&s);
            v.mul(&s_inv, &r);
            v.add_assign(&m);
            v.half();
        }

        if !m.sqrt(&s) {
            return false;
        }
        if (c & 5) == 0 || (c & 5) == 5 {
            let m_val = m;
            m.negate(&m_val, 1);
        }
        let c3 = ellswift_c3();
        let c4 = ellswift_c4();
        let u_val = u;
        u.mul(&u_val, if (c & 1) != 0 { &c4 } else { &c3 });
        u.add_assign(&v);
        t.mul(&m, &u);
        true
    }

    /// PRNG: SHA256(prefix || cnt_le) for deterministic randomness.
    fn ellswift_prng(prefix: &[u8], cnt: u32) -> [u8; 32] {
        let mut hasher = Sha256::new();
        hasher.update(prefix);
        hasher.update(cnt.to_le_bytes());
        hasher.finalize().into()
    }

    /// Find ElligatorSwift encoding (u, t) for point P using prefix as PRNG seed.
    fn elligatorswift_var(u32_out: &mut [u8; 32], t: &mut FieldElement, p: &Ge, prefix: &[u8]) {
        let mut cnt: u32 = 0;
        let mut branch_hash = [0u8; 32];
        let mut branches_left = 0usize;

        loop {
            if branches_left == 0 {
                branch_hash = ellswift_prng(prefix, cnt);
                cnt += 1;
                branches_left = 64;
            }
            branches_left -= 1;
            let branch = (branch_hash[branches_left >> 1] >> ((branches_left & 1) << 2)) & 7;

            *u32_out = ellswift_prng(prefix, cnt);
            cnt += 1;

            let mut u = FieldElement::zero();
            u.set_b32_mod(u32_out);

            if xswiftec_inv_var(t, &p.x, &u, branch as i32) {
                break;
            }
        }

        t.normalize_weak();
        let mut t_norm = *t;
        t_norm.normalize();
        let mut p_y_norm = p.y;
        p_y_norm.normalize();
        if t_norm.is_odd() != p_y_norm.is_odd() {
            let t_val = *t;
            t.negate(&t_val, 1);
            t.normalize_weak();
        }
    }

    /// Decode (u, t) to fraction (xn, xd) representing X coordinate.
    fn xswiftec_frac_var(
        xn: &mut FieldElement,
        xd: &mut FieldElement,
        u: &FieldElement,
        t: &FieldElement,
    ) {
        let mut u1 = *u;
        if u1.normalizes_to_zero_var() {
            u1.set_int(1);
        }
        let mut s = FieldElement::zero();
        s.sqr(t);
        if t.normalizes_to_zero_var() {
            s.set_int(1);
        }
        let mut l = FieldElement::zero();
        l.sqr(&u1);
        let mut g = FieldElement::zero();
        g.mul(&l, &u1);
        g.add_int(SECP256K1_B);
        let mut p = g;
        p.add_assign(&s);
        if p.normalizes_to_zero_var() {
            s.mul_int(4);
            p = g;
            p.add_assign(&s);
        }
        let mut d = FieldElement::zero();
        d.mul(&s, &l);
        d.mul_int(3);
        l.sqr(&p);
        let l_val = l;
        l.negate(&l_val, 1);
        let mut n = FieldElement::zero();
        n.mul(&d, &u1);
        n.add_assign(&l);
        if Ge::x_frac_on_curve_var(&n, &d) {
            *xn = n;
            *xd = d;
            return;
        }
        *xd = p;
        let c1 = ellswift_c1();
        let c2 = ellswift_c2();
        l.mul(&c1, &s);
        n.mul(&c2, &g);
        n.add_assign(&l);
        let n_val = n;
        n.mul(&n_val, &u1);
        if Ge::x_frac_on_curve_var(&n, &p) {
            *xn = n;
            return;
        }
        l.mul(&p, &u1);
        n.add_assign(&l);
        xn.negate(&n, 2);
    }

    /// Decode (u, t) to X coordinate.
    fn xswiftec_var(x: &mut FieldElement, u: &FieldElement, t: &FieldElement) {
        let mut xn = FieldElement::zero();
        let mut xd = FieldElement::zero();
        xswiftec_frac_var(&mut xn, &mut xd, u, t);
        let mut xd_inv = FieldElement::zero();
        xd_inv.inv(&xd);
        x.mul(&xn, &xd_inv);
    }

    /// Decode (u, t) to curve point. Y parity from t.
    fn swiftec_var(p: &mut Ge, u: &FieldElement, t: &FieldElement) {
        let mut x = FieldElement::zero();
        xswiftec_var(&mut x, u, t);
        let mut t_norm = *t;
        t_norm.normalize_weak();
        p.set_xo_var(&x, t_norm.is_odd());
    }

    /// Encode pubkey to 64-byte ElligatorSwift with random rnd32.
    /// pubkey33: 33-byte compressed pubkey. rnd32: 32-byte randomness.
    pub fn ellswift_encode(pubkey33: &[u8; 33], rnd32: &[u8; 32]) -> Option<[u8; 64]> {
        let p = ge_from_compressed(pubkey33)?;
        let mut prefix = Vec::with_capacity(160);
        let tag_hash = Sha256::digest(b"secp256k1_ellswift_encode");
        prefix.extend_from_slice(&tag_hash);
        prefix.extend_from_slice(&tag_hash);
        let mut p64 = [0u8; 64];
        p64[0..33].copy_from_slice(pubkey33);
        prefix.extend_from_slice(&p64);
        prefix.extend_from_slice(rnd32);
        let mut u32_out = [0u8; 32];
        let mut t = FieldElement::zero();
        elligatorswift_var(&mut u32_out, &mut t, &p, &prefix);
        let mut ell64 = [0u8; 64];
        ell64[0..32].copy_from_slice(&u32_out);
        t.normalize_weak();
        let mut t_bytes = [0u8; 32];
        t.get_b32(&mut t_bytes);
        ell64[32..64].copy_from_slice(&t_bytes);
        Some(ell64)
    }

    /// Create 64-byte ElligatorSwift from seckey. auxrnd32 optional.
    pub fn ellswift_create(seckey32: &[u8; 32], auxrnd32: Option<&[u8; 32]>) -> Option<[u8; 64]> {
        let mut d = Scalar::zero();
        if d.set_b32(seckey32) {
            return None;
        }
        if d.is_zero() {
            return None;
        }
        let p = pubkey_from_secret(&d);
        let mut prefix = Vec::with_capacity(160);
        let tag_hash = Sha256::digest(b"secp256k1_ellswift_create");
        prefix.extend_from_slice(&tag_hash);
        prefix.extend_from_slice(&tag_hash);
        prefix.extend_from_slice(seckey32);
        prefix.extend_from_slice(&[0u8; 32]);
        if let Some(aux) = auxrnd32 {
            prefix.extend_from_slice(aux);
        }
        let mut u32_out = [0u8; 32];
        let mut t = FieldElement::zero();
        elligatorswift_var(&mut u32_out, &mut t, &p, &prefix);
        let mut ell64 = [0u8; 64];
        ell64[0..32].copy_from_slice(&u32_out);
        t.normalize_weak();
        let mut t_bytes = [0u8; 32];
        t.get_b32(&mut t_bytes);
        ell64[32..64].copy_from_slice(&t_bytes);
        Some(ell64)
    }

    /// Decode 64-byte ElligatorSwift encoding to curve point.
    pub fn ellswift_decode(ell64: &[u8; 64]) -> Ge {
        let mut u = FieldElement::zero();
        let mut t = FieldElement::zero();
        let u_bytes: [u8; 32] = ell64[0..32].try_into().unwrap();
        let t_bytes: [u8; 32] = ell64[32..64].try_into().unwrap();
        u.set_b32_mod(&u_bytes);
        t.set_b32_mod(&t_bytes);
        t.normalize_weak();
        let mut p = Ge::default();
        swiftec_var(&mut p, &u, &t);
        p
    }

    /// BIP340 tagged hash for "bip324_ellswift_xonly_ecdh".
    fn tagged_hash_bip324_ellswift_xdh(data: &[u8]) -> [u8; 32] {
        let tag = b"bip324_ellswift_xonly_ecdh";
        let tag_hash = Sha256::digest(tag);
        let mut hasher = Sha256::new();
        hasher.update(tag_hash);
        hasher.update(tag_hash);
        hasher.update(data);
        hasher.finalize().into()
    }

    /// BIP324 x-only ECDH hash: H_tag(ell_a64 || ell_b64 || x32).
    pub fn ellswift_xdh_hash_bip324(
        output: &mut [u8; 32],
        x32: &[u8; 32],
        ell_a64: &[u8; 64],
        ell_b64: &[u8; 64],
    ) {
        let mut data = Vec::with_capacity(160);
        data.extend_from_slice(ell_a64);
        data.extend_from_slice(ell_b64);
        data.extend_from_slice(x32);
        *output = tagged_hash_bip324_ellswift_xdh(&data);
    }

    /// X-only ECDH using ElligatorSwift. party: false = A, true = B.
    /// seckey must correspond to our ell encoding (ell_a64 if party false, ell_b64 if true).
    /// Multiplication by `seckey` uses constant-time `ecmult_const`.
    pub fn ellswift_xdh(
        ell_a64: &[u8; 64],
        ell_b64: &[u8; 64],
        seckey: &[u8; 32],
        party: bool,
    ) -> Option<[u8; 32]> {
        let theirs = if party { ell_a64 } else { ell_b64 };
        let mut u = FieldElement::zero();
        let mut t = FieldElement::zero();
        u.set_b32_mod(&theirs[0..32].try_into().unwrap());
        t.set_b32_mod(&theirs[32..64].try_into().unwrap());
        let mut xn = FieldElement::zero();
        let mut xd = FieldElement::zero();
        xswiftec_frac_var(&mut xn, &mut xd, &u, &t);
        let mut xd_inv = FieldElement::zero();
        xd_inv.inv(&xd);
        let mut x = FieldElement::zero();
        x.mul(&xn, &xd_inv);
        let mut pt = Ge::default();
        if !pt.set_xo_var(&x, false) {
            return None;
        }
        let mut seckey_scalar = Scalar::zero();
        if seckey_scalar.set_b32(seckey) {
            return None;
        }
        if seckey_scalar.is_zero() {
            return None;
        }
        let mut res = Gej::default();
        ecmult_const(&mut res, &pt, &seckey_scalar);
        let mut res_ge = Ge::default();
        res_ge.set_gej_var(&res);
        if res_ge.infinity {
            return None;
        }
        res_ge.x.normalize();
        let mut x32 = [0u8; 32];
        res_ge.x.get_b32(&mut x32);
        let mut output = [0u8; 32];
        ellswift_xdh_hash_bip324(&mut output, &x32, ell_a64, ell_b64);
        Some(output)
    }
}

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
pub use imp::{
    ellswift_create, ellswift_decode, ellswift_encode, ellswift_xdh, ellswift_xdh_hash_bip324,
};
