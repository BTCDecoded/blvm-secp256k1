//! Group operations for secp256k1.
//!
//! Affine (Ge) and Jacobian (Gej) point representation.
//! Curve: y² = x³ + 7

use std::sync::OnceLock;
use subtle::Choice;

use crate::field::{FeStorage, FieldElement};

/// secp256k1 curve constant: y² = x³ + B
const SECP256K1_B: u32 = 7;

/// Beta: nontrivial cube root of 1 in Fp. (x,y) -> (beta*x, y) is an endomorphism.
pub(crate) fn const_beta() -> FieldElement {
    let b: [u8; 32] = [
        0x7a, 0xe9, 0x6a, 0x2b, 0x65, 0x7c, 0x07, 0x10, 0x6e, 0x64, 0x47, 0x9e, 0xac, 0x34, 0x34,
        0xe9, 0x9c, 0xf0, 0x49, 0x75, 0x12, 0xf5, 0x89, 0x95, 0xc1, 0x39, 0x6c, 0x28, 0x71, 0x95,
        0x01, 0xee,
    ];
    let mut r = FieldElement::zero();
    r.set_b32_mod(&b);
    r
}

/// Affine point (x, y) or infinity.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Ge {
    pub x: FieldElement,
    pub y: FieldElement,
    pub infinity: bool,
}

/// Compact storage for affine point. Matches libsecp256k1 ge_storage.
#[derive(Clone, Copy, Debug)]
pub struct GeStorage {
    pub x: FeStorage,
    pub y: FeStorage,
}

impl GeStorage {
    pub fn cmov(&mut self, a: &GeStorage, flag: Choice) {
        self.x.cmov(&a.x, flag);
        self.y.cmov(&a.y, flag);
    }
}

/// Jacobian point (X, Y, Z) representing (X/Z², Y/Z³) or infinity.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct Gej {
    pub x: FieldElement,
    pub y: FieldElement,
    pub z: FieldElement,
    pub infinity: bool,
}

impl Ge {
    pub fn set_xy(&mut self, x: &FieldElement, y: &FieldElement) {
        self.infinity = false;
        self.x = *x;
        self.y = *y;
    }

    pub fn set_infinity(&mut self) {
        self.infinity = true;
        self.x.set_int(0);
        self.y.set_int(0);
    }

    pub fn is_infinity(&self) -> bool {
        self.infinity
    }

    pub fn neg(&mut self, a: &Ge) {
        *self = *a;
        if !a.infinity {
            self.y.normalize_weak();
            let y = self.y;
            self.y.negate(&y, 1);
        }
    }

    /// Set from Jacobian. Modifies a (inverts z in place).
    pub fn set_gej(&mut self, a: &mut Gej) {
        if a.infinity {
            self.set_infinity();
            return;
        }
        self.infinity = false;
        let z = a.z;
        a.z.inv(&z);
        let mut z2 = FieldElement::zero();
        let mut z3 = FieldElement::zero();
        z2.sqr(&a.z);
        z3.mul(&z2, &a.z);
        let ax = a.x;
        let ay = a.y;
        a.x.mul(&ax, &z2);
        a.y.mul(&ay, &z3);
        a.z.set_int(1);
        self.x = a.x;
        self.y = a.y;
    }

    /// Set from Jacobian (variable-time, doesn't modify a).
    #[inline(always)]
    pub fn set_gej_var(&mut self, a: &Gej) {
        if a.is_infinity() {
            self.set_infinity();
            return;
        }
        self.infinity = false;
        let mut zi = FieldElement::zero();
        zi.inv(&a.z);
        let mut z2 = FieldElement::zero();
        let mut z3 = FieldElement::zero();
        z2.sqr(&zi);
        z3.mul(&z2, &zi);
        let mut x = FieldElement::zero();
        let mut y = FieldElement::zero();
        x.mul(&a.x, &z2);
        y.mul(&a.y, &z3);
        self.set_xy(&x, &y);
    }

    /// Set r = (a.x*zi², a.y*zi³). a must not be infinity.
    #[inline(always)]
    pub fn set_gej_zinv(&mut self, a: &Gej, zi: &FieldElement) {
        debug_assert!(!a.infinity);
        self.infinity = false;
        let mut zi2 = FieldElement::zero();
        let mut zi3 = FieldElement::zero();
        zi2.sqr(zi);
        zi3.mul(&zi2, zi);
        self.x.mul(&a.x, &zi2);
        self.y.mul(&a.y, &zi3);
    }

    /// Set r = (a.x*zi², a.y*zi³). a must not be infinity.
    pub fn set_ge_zinv(&mut self, a: &Ge, zi: &FieldElement) {
        debug_assert!(!a.infinity);
        self.infinity = false;
        let mut zi2 = FieldElement::zero();
        let mut zi3 = FieldElement::zero();
        zi2.sqr(zi);
        zi3.mul(&zi2, zi);
        self.x.mul(&a.x, &zi2);
        self.y.mul(&a.y, &zi3);
    }

    /// Endomorphism: r = (beta*a.x, a.y). lambda*(x,y) = (beta*x, y).
    pub fn mul_lambda(&mut self, a: &Ge) {
        *self = *a;
        if !a.infinity {
            let beta = const_beta();
            let x = self.x;
            self.x.mul(&x, &beta);
        }
    }

    /// Set from compact storage.
    pub fn from_storage(&mut self, a: &GeStorage) {
        self.infinity = false;
        self.x.from_storage(&a.x);
        self.y.from_storage(&a.y);
    }

    /// Check if x is a valid X coordinate on the curve (y² = x³ + 7).
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    pub fn x_on_curve_var(x: &FieldElement) -> bool {
        let mut c = FieldElement::zero();
        c.sqr(x);
        let c_val = c;
        c.mul(&c_val, x);
        c.add_int(SECP256K1_B);
        FieldElement::is_square_var(&c)
    }

    /// Check if fraction xn/xd is a valid X coordinate on the curve.
    /// xd must not be zero.
    #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
    pub fn x_frac_on_curve_var(xn: &FieldElement, xd: &FieldElement) -> bool {
        let mut r = FieldElement::zero();
        let mut t = FieldElement::zero();
        r.mul(xd, xn);
        t.sqr(xn);
        let r_val = r;
        r.mul(&r_val, &t);
        t.sqr(xd);
        let t_val = t;
        t.sqr(&t_val);
        t.mul_int(SECP256K1_B);
        r.add_assign(&t);
        FieldElement::is_square_var(&r)
    }

    /// Set from x coordinate and y oddness. Returns true if x is on curve.
    #[inline(always)]
    pub fn set_xo_var(&mut self, x: &FieldElement, odd: bool) -> bool {
        let mut x2 = FieldElement::zero();
        let mut x3 = FieldElement::zero();
        x2.sqr(x);
        x3.mul(x, &x2);
        self.x = *x;
        self.infinity = false;
        x3.add_int(SECP256K1_B);
        let ok = self.y.sqrt(&x3);
        self.y.normalize();
        if self.y.is_odd() != odd {
            let y = self.y;
            self.y.negate(&y, 1);
        }
        ok
    }
}

impl Gej {
    pub fn set_infinity(&mut self) {
        self.infinity = true;
        self.x.set_int(0);
        self.y.set_int(0);
        self.z.set_int(0);
    }

    pub fn is_infinity(&self) -> bool {
        self.infinity
    }

    #[inline(always)]
    pub fn set_ge(&mut self, a: &Ge) {
        self.infinity = a.infinity;
        self.x = a.x;
        self.y = a.y;
        self.z.set_int(1);
    }

    pub fn neg(&mut self, a: &Gej) {
        self.infinity = a.infinity;
        self.x = a.x;
        self.y = a.y;
        self.z = a.z;
        if !a.infinity {
            self.y.normalize_weak();
            let y = self.y;
            self.y.negate(&y, 1);
        }
    }

    #[inline(always)]
    pub fn double(&mut self, a: &Gej) {
        let mut l = FieldElement::zero();
        let mut s = FieldElement::zero();
        let mut t = FieldElement::zero();

        self.infinity = a.infinity;
        if a.infinity {
            return;
        }

        self.z.mul(&a.z, &a.y);
        s.sqr(&a.y);
        l.sqr(&a.x);
        l.mul_int(3);
        l.half();
        t.negate(&s, 1);
        let t_in = t;
        t.mul(&t_in, &a.x);
        self.x.sqr(&l);
        self.x.add_assign(&t);
        self.x.add_assign(&t);
        let s_in = s;
        s.sqr(&s_in);
        t.add_assign(&self.x);
        self.y.mul(&t, &l);
        self.y.add_assign(&s);
        let y_in = self.y;
        self.y.negate(&y_in, 2);
    }

    #[inline(always)]
    pub fn double_var(&mut self, a: &Gej) {
        if a.infinity {
            self.set_infinity();
            return;
        }
        self.double(a);
    }

    #[inline(always)]
    pub fn add_var(&mut self, a: &Gej, b: &Gej) {
        if a.infinity {
            *self = *b;
            return;
        }
        if b.infinity {
            *self = *a;
            return;
        }

        let mut z22 = FieldElement::zero();
        let mut z12 = FieldElement::zero();
        let mut u1 = FieldElement::zero();
        let mut u2 = FieldElement::zero();
        let mut s1 = FieldElement::zero();
        let mut s2 = FieldElement::zero();
        let mut h = FieldElement::zero();
        let mut i = FieldElement::zero();
        let mut h2 = FieldElement::zero();
        let mut h3 = FieldElement::zero();
        let mut t = FieldElement::zero();

        z22.sqr(&b.z);
        z12.sqr(&a.z);
        u1.mul(&a.x, &z22);
        u2.mul(&b.x, &z12);
        s1.mul(&a.y, &z22);
        let s1_in = s1;
        s1.mul(&s1_in, &b.z);
        s2.mul(&b.y, &z12);
        let s2_in = s2;
        s2.mul(&s2_in, &a.z);
        h.negate(&u1, 1);
        h.add_assign(&u2);
        i.negate(&s2, 1);
        i.add_assign(&s1);

        if h.normalizes_to_zero_var() {
            if i.normalizes_to_zero_var() {
                self.double_var(a);
            } else {
                self.set_infinity();
            }
            return;
        }

        self.infinity = false;
        t.mul(&h, &b.z);
        self.z.mul(&a.z, &t);
        h2.sqr(&h);
        let h2_in = h2;
        h2.negate(&h2_in, 1);
        h3.mul(&h2, &h);
        t.mul(&u1, &h2);
        self.x.sqr(&i);
        self.x.add_assign(&h3);
        self.x.add_assign(&t);
        self.x.add_assign(&t);
        t.add_assign(&self.x);
        self.y.mul(&t, &i);
        let h3_in = h3;
        h3.mul(&h3_in, &s1);
        self.y.add_assign(&h3);
    }

    /// Rescale: r = (r.x*s², r.y*s³, r.z*s). s must not be zero.
    pub fn rescale(&mut self, s: &FieldElement) {
        debug_assert!(!s.normalizes_to_zero_var());
        let mut zz = FieldElement::zero();
        zz.sqr(s);
        let x = self.x;
        let y = self.y;
        let z = self.z;
        self.x.mul(&x, &zz);
        self.y.mul(&y, &zz);
        let y_zz = self.y;
        self.y.mul(&y_zz, s);
        self.z.mul(&z, s);
    }

    /// r = a + b with b's z known as bzinv (1/z). For adding affine to Jacobian.
    #[inline(always)]
    pub fn add_zinv_var(&mut self, a: &Gej, b: &Ge, bzinv: &FieldElement) {
        if a.infinity {
            self.infinity = b.infinity;
            if !b.infinity {
                let mut bzinv2 = FieldElement::zero();
                let mut bzinv3 = FieldElement::zero();
                bzinv2.sqr(bzinv);
                bzinv3.mul(&bzinv2, bzinv);
                self.x.mul(&b.x, &bzinv2);
                self.y.mul(&b.y, &bzinv3);
                self.z.set_int(1);
            }
            return;
        }
        if b.infinity {
            *self = *a;
            return;
        }
        let mut az = FieldElement::zero();
        az.mul(&a.z, bzinv);
        let mut z12 = FieldElement::zero();
        z12.sqr(&az);
        let u1 = a.x;
        let mut u2 = FieldElement::zero();
        u2.mul(&b.x, &z12);
        let s1 = a.y;
        let mut s2 = FieldElement::zero();
        s2.mul(&b.y, &z12);
        let s2_in = s2;
        s2.mul(&s2_in, &az);
        let mut h = FieldElement::zero();
        h.negate(&u1, 4);
        h.add_assign(&u2);
        let mut i = FieldElement::zero();
        i.negate(&s2, 1);
        i.add_assign(&s1);
        if h.normalizes_to_zero_var() {
            if i.normalizes_to_zero_var() {
                self.double_var(a);
            } else {
                self.set_infinity();
            }
            return;
        }
        self.infinity = false;
        self.z.mul(&a.z, &h);
        let mut h2 = FieldElement::zero();
        h2.sqr(&h);
        let h2_in = h2;
        h2.negate(&h2_in, 1);
        let mut h3 = FieldElement::zero();
        h3.mul(&h2, &h);
        let mut t = FieldElement::zero();
        t.mul(&u1, &h2);
        self.x.sqr(&i);
        self.x.add_assign(&h3);
        self.x.add_assign(&t);
        self.x.add_assign(&t);
        t.add_assign(&self.x);
        self.y.mul(&t, &i);
        let h3_in = h3;
        h3.mul(&h3_in, &s1);
        self.y.add_assign(&h3);
    }

    /// add_ge_var with optional z-ratio output (for odd_multiples_table).
    #[inline(always)]
    pub fn add_ge_var_rzr(&mut self, a: &Gej, b: &Ge, rzr: Option<&mut FieldElement>) {
        if a.infinity {
            if let Some(rzr) = rzr {
                rzr.set_int(1);
            }
            self.set_ge(b);
            return;
        }
        if b.infinity {
            if let Some(rzr) = rzr {
                rzr.set_int(0);
            }
            *self = *a;
            return;
        }

        let mut z12 = FieldElement::zero();
        let u1 = a.x;
        let mut u2 = FieldElement::zero();
        let s1 = a.y;
        let mut s2 = FieldElement::zero();
        let mut h = FieldElement::zero();
        let mut i = FieldElement::zero();
        let mut h2 = FieldElement::zero();
        let mut h3 = FieldElement::zero();
        let mut t = FieldElement::zero();

        z12.sqr(&a.z);
        u2.mul(&b.x, &z12);
        s2.mul(&b.y, &z12);
        let s2_in = s2;
        s2.mul(&s2_in, &a.z);
        h.negate(&u1, 4); // GEJ_X_MAGNITUDE_MAX
        h.add_assign(&u2);
        i.negate(&s2, 1);
        i.add_assign(&s1);

        if h.normalizes_to_zero_var() {
            if let Some(rzr) = rzr {
                rzr.set_int(0);
            }
            if i.normalizes_to_zero_var() {
                self.double_var(a);
            } else {
                self.set_infinity();
            }
            return;
        }

        self.infinity = false;
        if let Some(rzr) = rzr {
            *rzr = h;
        }
        self.z.mul(&a.z, &h);
        h2.sqr(&h);
        let h2_in = h2;
        h2.negate(&h2_in, 1);
        h3.mul(&h2, &h);
        t.mul(&u1, &h2);
        self.x.sqr(&i);
        self.x.add_assign(&h3);
        self.x.add_assign(&t);
        self.x.add_assign(&t);
        t.add_assign(&self.x);
        self.y.mul(&t, &i);
        let h3_in = h3;
        h3.mul(&h3_in, &s1);
        self.y.add_assign(&h3);
    }

    #[inline(always)]
    pub fn add_ge_var(&mut self, a: &Gej, b: &Ge) {
        self.add_ge_var_rzr(a, b, None);
    }

    /// Check whether affine x of this Jacobian point equals the given x.
    /// Returns x * z² == self.x (avoids inversion in ECDSA verify).
    pub fn eq_x_var(&self, x: &FieldElement) -> bool {
        debug_assert!(!self.infinity);
        let mut z2 = FieldElement::zero();
        z2.sqr(&self.z);
        let mut r = FieldElement::zero();
        r.mul(&z2, x);
        FieldElement::fe_equal(&r, &self.x)
    }

    /// Jacobian double (libsecp256k1 `gej_double`): same formula for all inputs; no early return.
    /// Used on secret-dependent paths with `ecmult_const`.
    pub fn double_ct(&mut self, a: &Gej) {
        let mut l = FieldElement::zero();
        let mut s = FieldElement::zero();
        let mut t = FieldElement::zero();
        self.infinity = a.infinity;
        self.z.mul(&a.z, &a.y);
        s.sqr(&a.y);
        l.sqr(&a.x);
        l.mul_int(3);
        l.half();
        t.negate(&s, 1);
        let t_in = t;
        t.mul(&t_in, &a.x);
        self.x.sqr(&l);
        self.x.add_assign(&t);
        self.x.add_assign(&t);
        let s_in = s;
        s.sqr(&s_in);
        t.add_assign(&self.x);
        self.y.mul(&t, &l);
        self.y.add_assign(&s);
        let y_in = self.y;
        self.y.negate(&y_in, 2);
        self.infinity = a.infinity;
    }

    /// Mixed add (libsecp256k1 `gej_add_ge`): Brier–Joye unified formula, constant-time.
    /// `b` must not be infinity.
    pub fn add_ge(&mut self, a: &Gej, b: &Ge) {
        debug_assert!(!b.infinity);
        let mut zz = FieldElement::zero();
        zz.sqr(&a.z);
        let u1 = a.x;
        let mut u2 = FieldElement::zero();
        u2.mul(&b.x, &zz);
        let s1 = a.y;
        let mut s2 = FieldElement::zero();
        let mut rr = FieldElement::zero();
        let mut m_alt = FieldElement::zero();
        let mut tt = FieldElement::zero();
        s2.mul(&b.y, &zz);
        let s2_in = s2;
        s2.mul(&s2_in, &a.z);
        let mut t = u1;
        t.add_assign(&u2);
        let mut m = s1;
        m.add_assign(&s2);
        rr.sqr(&t);
        m_alt.negate(&u2, 1);
        tt.mul(&u1, &m_alt);
        rr.add_assign(&tt);
        let degenerate = m.normalizes_to_zero();
        let mut rr_alt = s1;
        rr_alt.mul_int(2);
        m_alt.add_assign(&u1);
        let not_deg = !degenerate as u8;
        rr_alt.cmov(&rr, Choice::from(not_deg));
        m_alt.cmov(&m, Choice::from(not_deg));
        let mut n = FieldElement::zero();
        n.sqr(&m_alt);
        let mut q = FieldElement::zero();
        q.negate(&t, 5);
        let q0 = q;
        q.mul(&q0, &n);
        n.sqr_assign();
        n.cmov(&m, Choice::from(degenerate as u8));
        t.sqr(&rr_alt);
        self.z.mul(&a.z, &m_alt);
        t.add_assign(&q);
        self.x = t;
        t.mul_int(2);
        t.add_assign(&q);
        let t0 = t;
        t.mul(&t0, &rr_alt);
        t.add_assign(&n);
        self.y.negate(&t, 5);
        self.y.half();
        self.x.cmov(&b.x, Choice::from(a.infinity as u8));
        self.y.cmov(&b.y, Choice::from(a.infinity as u8));
        let mut one = FieldElement::zero();
        one.set_int(1);
        self.z.cmov(&one, Choice::from(a.infinity as u8));
        self.infinity = self.z.normalizes_to_zero();
    }
}

/// Batch convert Gej -> Ge using one field inversion instead of N.
/// Uses r[i].x as scratch for z-products; after: r[i] = affine(a[i]).
pub fn ge_set_all_gej_var(r: &mut [Ge], a: &[Gej]) {
    let len = r.len().min(a.len());
    let mut last_i: Option<usize> = None;

    for i in 0..len {
        if a[i].infinity {
            r[i].set_infinity();
        } else {
            if let Some(li) = last_i {
                let r_li_x = r[li].x;
                let a_i_z = a[i].z;
                r[i].x.mul(&r_li_x, &a_i_z);
            } else {
                r[i].x = a[i].z;
            }
            r[i].infinity = false; // will be overwritten by set_gej_zinv
            last_i = Some(i);
        }
    }

    let Some(mut last_i) = last_i else {
        return;
    };

    let mut u = FieldElement::zero();
    u.inv(&r[last_i].x);

    let mut i = last_i;
    while i > 0 {
        i -= 1;
        if !a[i].infinity {
            let r_i_x = r[i].x;
            r[last_i].x.mul(&r_i_x, &u);
            let u_in = u;
            let a_last_z = a[last_i].z;
            u.mul(&u_in, &a_last_z);
            last_i = i;
        }
    }
    r[last_i].x = u;

    for i in 0..len {
        if !a[i].infinity {
            let zi = r[i].x;
            r[i].set_gej_zinv(&a[i], &zi);
        }
    }
}

/// Bring pre_a[0..len] to same global z. zr[i] = z(pre_a[i])/z(pre_a[i-1]), zr[0] unused.
/// Uses r[i].x as scratch. After: z(pre_a[i]) = z(pre_a[len-1]) for all i.
#[inline(always)]
pub fn ge_table_set_globalz(len: usize, pre_a: &mut [Ge], zr: &[FieldElement]) {
    if len == 0 {
        return;
    }
    let mut i = len - 1;
    pre_a[i].y.normalize_weak();
    let mut zs = zr[i];
    while i > 0 {
        if i != len - 1 {
            let zs_in = zs;
            zs.mul(&zs_in, &zr[i]);
        }
        i -= 1;
        let ai = pre_a[i];
        pre_a[i].set_ge_zinv(&ai, &zs);
    }
}

/// secp256k1 generator G (SEC2 2.7.1). Cached for performance.
pub fn generator_g() -> Ge {
    *GENERATOR_G.get_or_init(|| {
        let gx = [
            0x79, 0xbe, 0x66, 0x7e, 0xf9, 0xdc, 0xbb, 0xac, 0x55, 0xa0, 0x62, 0x95, 0xce, 0x87,
            0x0b, 0x07, 0x02, 0x9b, 0xfc, 0xdb, 0x2d, 0xce, 0x28, 0xd9, 0x59, 0xf2, 0x81, 0x5b,
            0x16, 0xf8, 0x17, 0x98,
        ];
        let gy = [
            0x48, 0x3a, 0xda, 0x77, 0x26, 0xa3, 0xc4, 0x65, 0x5d, 0xa4, 0xfb, 0xfc, 0x0e, 0x11,
            0x08, 0xa8, 0xfd, 0x17, 0xb4, 0x48, 0xa6, 0x85, 0x54, 0x19, 0x9c, 0x47, 0xd0, 0x8f,
            0xfb, 0x10, 0xd4, 0xb8,
        ];
        let mut x = FieldElement::zero();
        let mut y = FieldElement::zero();
        x.set_b32_mod(&gx);
        y.set_b32_mod(&gy);
        let mut g = Ge {
            x: FieldElement::zero(),
            y: FieldElement::zero(),
            infinity: false,
        };
        g.set_xy(&x, &y);
        g
    })
}

static GENERATOR_G: OnceLock<Ge> = OnceLock::new();
