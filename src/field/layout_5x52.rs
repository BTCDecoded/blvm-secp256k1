//! 5x52 field element layout for x86_64 and aarch64.
//!
//! Pure Rust implementation using u128 for wide multiplication.
//! Field modulus: p = 2^256 - 2^32 - 977 (SEC2 secp256k1)

use crate::modinv64::{self, SECP256K1_FE_MODINV_MODINFO, Signed62};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// Field element in 5x52 limb representation.
/// Each limb holds 52 bits (except limb 4 which holds 48).
/// Layout matches libsecp256k1 field_5x52.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default)]
pub struct FieldElement {
    pub n: [u64; 5],
}

const M: u64 = 0xFFFFFFFFFFFFF; // 52 bits
const M4: u64 = 0x0FFFFFFFFFFFF; // 48 bits for limb 4
const R: u64 = 0x1000003D10;
const R_LOW: u64 = 0x1000003D1; // p = 2^256 - 2^32 - 977
const P0: u64 = 0xFFFFEFFFFFC2F; // high limb of p for comparison

impl FieldElement {
    pub fn zero() -> Self {
        Self { n: [0, 0, 0, 0, 0] }
    }

    pub fn one() -> Self {
        Self { n: [1, 0, 0, 0, 0] }
    }

    /// Set to a small integer a in [0, 0x7FFF].
    pub fn set_int(&mut self, a: u32) {
        debug_assert!(a <= 0x7FFF);
        self.n[0] = a as u64;
        self.n[1] = 0;
        self.n[2] = 0;
        self.n[3] = 0;
        self.n[4] = 0;
    }

    /// True if the normalized value is odd. Caller must ensure magnitude <= 1.
    pub fn is_odd(&self) -> bool {
        self.n[0] & 1 == 1
    }

    /// Set from 32-byte big-endian encoding (mod p).
    /// Returns false if the value is >= p.
    pub fn set_b32_limit(&mut self, bytes: &[u8; 32]) -> bool {
        self.set_b32_mod(bytes);
        // Check if value >= p
        !((self.n[4] == M4) && (self.n[3] & self.n[2] & self.n[1]) == M && (self.n[0] >= P0))
    }

    /// Set from 32-byte big-endian encoding, reducing mod p.
    pub fn set_b32_mod(&mut self, a: &[u8; 32]) {
        self.n[0] = u64::from(a[31])
            | (u64::from(a[30]) << 8)
            | (u64::from(a[29]) << 16)
            | (u64::from(a[28]) << 24)
            | (u64::from(a[27]) << 32)
            | (u64::from(a[26]) << 40)
            | ((u64::from(a[25]) & 0xF) << 48);
        self.n[1] = ((a[25] >> 4) & 0xF) as u64
            | (u64::from(a[24]) << 4)
            | (u64::from(a[23]) << 12)
            | (u64::from(a[22]) << 20)
            | (u64::from(a[21]) << 28)
            | (u64::from(a[20]) << 36)
            | (u64::from(a[19]) << 44);
        self.n[2] = u64::from(a[18])
            | (u64::from(a[17]) << 8)
            | (u64::from(a[16]) << 16)
            | (u64::from(a[15]) << 24)
            | (u64::from(a[14]) << 32)
            | (u64::from(a[13]) << 40)
            | ((u64::from(a[12]) & 0xF) << 48);
        self.n[3] = ((a[12] >> 4) & 0xF) as u64
            | (u64::from(a[11]) << 4)
            | (u64::from(a[10]) << 12)
            | (u64::from(a[9]) << 20)
            | (u64::from(a[8]) << 28)
            | (u64::from(a[7]) << 36)
            | (u64::from(a[6]) << 44);
        self.n[4] = u64::from(a[5])
            | (u64::from(a[4]) << 8)
            | (u64::from(a[3]) << 16)
            | (u64::from(a[2]) << 24)
            | (u64::from(a[1]) << 32)
            | (u64::from(a[0]) << 40);
    }

    /// Write normalized value to 32-byte big-endian.
    pub fn get_b32(&self, r: &mut [u8; 32]) {
        let a = &self.n;
        r[0] = (a[4] >> 40) as u8;
        r[1] = (a[4] >> 32) as u8;
        r[2] = (a[4] >> 24) as u8;
        r[3] = (a[4] >> 16) as u8;
        r[4] = (a[4] >> 8) as u8;
        r[5] = a[4] as u8;
        r[6] = (a[3] >> 44) as u8;
        r[7] = (a[3] >> 36) as u8;
        r[8] = (a[3] >> 28) as u8;
        r[9] = (a[3] >> 20) as u8;
        r[10] = (a[3] >> 12) as u8;
        r[11] = (a[3] >> 4) as u8;
        r[12] = (((a[2] >> 48) & 0xF) | ((a[3] & 0xF) << 4)) as u8;
        r[13] = (a[2] >> 40) as u8;
        r[14] = (a[2] >> 32) as u8;
        r[15] = (a[2] >> 24) as u8;
        r[16] = (a[2] >> 16) as u8;
        r[17] = (a[2] >> 8) as u8;
        r[18] = a[2] as u8;
        r[19] = (a[1] >> 44) as u8;
        r[20] = (a[1] >> 36) as u8;
        r[21] = (a[1] >> 28) as u8;
        r[22] = (a[1] >> 20) as u8;
        r[23] = (a[1] >> 12) as u8;
        r[24] = (a[1] >> 4) as u8;
        r[25] = (((a[0] >> 48) & 0xF) | ((a[1] & 0xF) << 4)) as u8;
        r[26] = (a[0] >> 40) as u8;
        r[27] = (a[0] >> 32) as u8;
        r[28] = (a[0] >> 24) as u8;
        r[29] = (a[0] >> 16) as u8;
        r[30] = (a[0] >> 8) as u8;
        r[31] = a[0] as u8;
    }

    /// Normalize to canonical form (constant-time).
    #[inline(always)]
    pub fn normalize(&mut self) {
        let mut t0 = self.n[0];
        let mut t1 = self.n[1];
        let mut t2 = self.n[2];
        let mut t3 = self.n[3];
        let mut t4 = self.n[4];

        let mut x = t4 >> 48;
        t4 &= M4;

        t0 += x * R_LOW;
        t1 += t0 >> 52;
        t0 &= M;
        let mut m = t1;
        t2 += t1 >> 52;
        t1 &= M;
        m &= t2;
        t3 += t2 >> 52;
        t2 &= M;
        m &= t3;
        t4 += t3 >> 52;
        t3 &= M;
        m &= t4;

        // For value p we have t4=M4, m=M4 (not M). Check (m == M) | (m == M4) for p.
        x = (t4 >> 48) | (((t4 == M4) as u64) & ((m == M || m == M4) as u64) & ((t0 >= P0) as u64));

        t0 += x * R_LOW;
        t1 += t0 >> 52;
        t0 &= M;
        t2 += t1 >> 52;
        t1 &= M;
        t3 += t2 >> 52;
        t2 &= M;
        t4 += t3 >> 52;
        t3 &= M;
        t4 &= M4;

        self.n[0] = t0;
        self.n[1] = t1;
        self.n[2] = t2;
        self.n[3] = t3;
        self.n[4] = t4;
    }

    pub fn is_zero(&self) -> bool {
        (self.n[0] | self.n[1] | self.n[2] | self.n[3] | self.n[4]) == 0
    }

    /// Add a small integer to r.n[0]. `a` must be in [0, 0x7FFF].
    pub fn add_int(&mut self, a: u32) {
        debug_assert!(a <= 0x7FFF);
        self.n[0] += a as u64;
    }

    /// Multiply all limbs by a small integer.
    #[inline(always)]
    pub fn mul_int(&mut self, a: u32) {
        let a = a as u64;
        self.n[0] *= a;
        self.n[1] *= a;
        self.n[2] *= a;
        self.n[3] *= a;
        self.n[4] *= a;
    }

    /// Half: r = r/2 mod p. If r is odd, add (p+1)/2 first.
    pub fn half(&mut self) {
        let mut t0 = self.n[0];
        let mut t1 = self.n[1];
        let mut t2 = self.n[2];
        let mut t3 = self.n[3];
        let mut t4 = self.n[4];
        let one = 1u64;
        let mask = if t0 & one == 1 { u64::MAX >> 12 } else { 0 };
        t0 += P0 & mask;
        t1 += mask;
        t2 += mask;
        t3 += mask;
        t4 += mask >> 4;
        self.n[0] = (t0 >> 1) + ((t1 & one) << 51);
        self.n[1] = (t1 >> 1) + ((t2 & one) << 51);
        self.n[2] = (t2 >> 1) + ((t3 & one) << 51);
        self.n[3] = (t3 >> 1) + ((t4 & one) << 51);
        self.n[4] = t4 >> 1;
    }

    /// Normalize weakly: propagate carries but skip final reduction check.
    /// Result may be in [0, p) or exactly p (represented as overflow in limb 4).
    pub fn normalize_weak(&mut self) {
        let mut t0 = self.n[0];
        let mut t1 = self.n[1];
        let mut t2 = self.n[2];
        let mut t3 = self.n[3];
        let mut t4 = self.n[4];
        let x = t4 >> 48;
        t4 &= M4;
        t0 += x * R_LOW;
        t1 += t0 >> 52;
        t0 &= M;
        t2 += t1 >> 52;
        t1 &= M;
        t3 += t2 >> 52;
        t2 &= M;
        t4 += t3 >> 52;
        t3 &= M;
        self.n[0] = t0;
        self.n[1] = t1;
        self.n[2] = t2;
        self.n[3] = t3;
        self.n[4] = t4;
    }

    /// Returns true if this value normalizes to 0 (either raw 0 or p).
    pub fn normalizes_to_zero(&self) -> bool {
        let mut t0 = self.n[0];
        let mut t1 = self.n[1];
        let mut t2 = self.n[2];
        let mut t3 = self.n[3];
        let mut t4 = self.n[4];
        let x = t4 >> 48;
        t4 &= M4;
        t0 += x * R_LOW;
        t1 += t0 >> 52;
        t0 &= M;
        let mut z0 = t0;
        let mut z1 = t0 ^ 0x1000003D0;
        t2 += t1 >> 52;
        t1 &= M;
        z0 |= t1;
        z1 &= t1;
        t3 += t2 >> 52;
        t2 &= M;
        z0 |= t2;
        z1 &= t2;
        t4 += t3 >> 52;
        t3 &= M;
        z0 |= t3;
        z1 &= t3;
        z0 |= t4;
        z1 &= t4 ^ 0xF000000000000;
        (z0 == 0) || (z1 == M)
    }

    /// Variable-time version of normalizes_to_zero (early exit).
    pub fn normalizes_to_zero_var(&self) -> bool {
        let t0 = self.n[0];
        let t4 = self.n[4];
        let x = t4 >> 48;
        let t0_new = t0 + x * R_LOW;
        let mut z0 = t0_new & M;
        let mut z1 = z0 ^ 0x1000003D0;
        if z0 != 0 && z1 != M {
            return false;
        }
        let mut t1 = self.n[1];
        let mut t2 = self.n[2];
        let mut t3 = self.n[3];
        let mut t4 = t4 & M4;
        t1 += t0_new >> 52;
        t2 += t1 >> 52;
        t1 &= M;
        z0 |= t1;
        z1 &= t1;
        t3 += t2 >> 52;
        t2 &= M;
        z0 |= t2;
        z1 &= t2;
        t4 += t3 >> 52;
        t3 &= M;
        z0 |= t3;
        z1 &= t3;
        z0 |= t4;
        z1 &= t4 ^ 0xF000000000000;
        (z0 == 0) || (z1 == M)
    }

    /// Returns true if a is a quadratic residue (square) mod p.
    /// Variable-time. Uses sqrt internally.
    pub fn is_square_var(a: &FieldElement) -> bool {
        let mut t = FieldElement::zero();
        t.sqrt(a)
    }

    /// Square root: r = sqrt(a) mod p if it exists.
    /// Returns true if a square root exists. Uses a^((p+1)/4) since p ≡ 3 (mod 4).
    #[inline]
    pub fn sqrt(&mut self, a: &FieldElement) -> bool {
        // Addition chain for (p+1)/4 with blocks {2, 22, 223}
        let mut x2 = FieldElement::zero();
        let mut x3 = FieldElement::zero();
        let mut x6;
        let mut x9;
        let mut x11;
        let mut x22;
        let mut x44;
        let mut x88;
        let mut x176;
        let mut x220;
        let mut x223;
        let mut t1;

        x2.sqr(a);
        let mut t = FieldElement::zero();
        t.mul(&x2, a);
        x2 = t;
        x3.sqr(&x2);
        t.mul(&x3, a);
        x3 = t;
        // Same addition chain as libsecp256k1; use in-place sqr to match their C.
        x6 = x3;
        for _ in 0..3 {
            x6.sqr_assign();
        }
        t.mul(&x6, &x3);
        x6 = t;
        x9 = x6;
        for _ in 0..3 {
            x9.sqr_assign();
        }
        t.mul(&x9, &x3);
        x9 = t;
        x11 = x9;
        for _ in 0..2 {
            x11.sqr_assign();
        }
        t.mul(&x11, &x2);
        x11 = t;
        x22 = x11;
        for _ in 0..11 {
            x22.sqr_assign();
        }
        t.mul(&x22, &x11);
        x22 = t;
        x44 = x22;
        for _ in 0..22 {
            x44.sqr_assign();
        }
        t.mul(&x44, &x22);
        x44 = t;
        x88 = x44;
        for _ in 0..44 {
            x88.sqr_assign();
        }
        t.mul(&x88, &x44);
        x88 = t;
        x176 = x88;
        for _ in 0..88 {
            x176.sqr_assign();
        }
        t.mul(&x176, &x88);
        x176 = t;
        x220 = x176;
        for _ in 0..44 {
            x220.sqr_assign();
        }
        t.mul(&x220, &x44);
        x220 = t;
        x223 = x220;
        for _ in 0..3 {
            x223.sqr_assign();
        }
        t.mul(&x223, &x3);
        t1 = t;
        for _ in 0..23 {
            t1.sqr_assign();
        }
        t.mul(&t1, &x22);
        t1 = t;
        for _ in 0..6 {
            t1.sqr_assign();
        }
        t.mul(&t1, &x2);
        t1 = t;
        t1.sqr_assign();
        self.sqr(&t1);

        let mut check = FieldElement::zero();
        check.sqr(self);
        FieldElement::fe_equal(&check, a)
    }

    /// Field equality (works with unnormalized): a == b iff (-a + b) normalizes to zero.
    pub fn fe_equal(a: &FieldElement, b: &FieldElement) -> bool {
        let mut na = *a;
        let a_copy = *a;
        na.negate(&a_copy, 1);
        na.add_assign(b);
        na.normalizes_to_zero_var()
    }

    /// Add another field element (in place).
    #[inline(always)]
    pub fn add_assign(&mut self, a: &FieldElement) {
        self.n[0] += a.n[0];
        self.n[1] += a.n[1];
        self.n[2] += a.n[2];
        self.n[3] += a.n[3];
        self.n[4] += a.n[4];
    }

    /// Add: r = a + b.
    #[inline(always)]
    pub fn add(&mut self, a: &FieldElement, b: &FieldElement) {
        *self = *a;
        self.add_assign(b);
    }

    /// Compare normalized field elements. Returns -1 if a < b, 0 if a == b, 1 if a > b.
    pub fn cmp_var(a: &FieldElement, b: &FieldElement) -> i32 {
        let mut na = *a;
        let mut nb = *b;
        na.normalize();
        nb.normalize();
        for i in (0..5).rev() {
            if na.n[i] < nb.n[i] {
                return -1;
            }
            if na.n[i] > nb.n[i] {
                return 1;
            }
        }
        0
    }

    /// Negate (r = -a mod p). Assumes magnitude <= 32.
    #[inline(always)]
    pub fn negate(&mut self, a: &FieldElement, m: u32) {
        let m = m as u64;
        self.n[0] = P0 * 2 * (m + 1) - a.n[0];
        self.n[1] = M * 2 * (m + 1) - a.n[1];
        self.n[2] = M * 2 * (m + 1) - a.n[2];
        self.n[3] = M * 2 * (m + 1) - a.n[3];
        self.n[4] = M4 * 2 * (m + 1) - a.n[4];
    }

    /// Multiply: r = a * b (mod p)
    #[inline(always)]
    pub fn mul(&mut self, a: &FieldElement, b: &FieldElement) {
        fe_mul_inner(&mut self.n, &a.n, &b.n);
    }

    /// Square: r = a^2 (mod p)
    #[inline(always)]
    pub fn sqr(&mut self, a: &FieldElement) {
        fe_sqr_inner(&mut self.n, &a.n);
    }

    /// Square in place: self = self^2. Matches libsecp256k1's in-place sqr in sqrt loops.
    #[inline(always)]
    pub fn sqr_assign(&mut self) {
        let a = *self;
        self.sqr(&a);
    }

    /// Modular inverse via modinv64 (Bernstein-Yang safegcd). Constant-time.
    #[inline(always)]
    pub fn inv(&mut self, a: &FieldElement) {
        let mut tmp = *a;
        tmp.normalize();
        let mut s = fe_to_signed62(&tmp);
        modinv64::modinv64(&mut s, &SECP256K1_FE_MODINV_MODINFO);
        fe_from_signed62(self, &s);
    }

    /// Constant-time conditional move: if flag, set self = a.
    pub fn cmov(&mut self, a: &FieldElement, flag: Choice) {
        for i in 0..5 {
            self.n[i] = u64::conditional_select(&self.n[i], &a.n[i], flag);
        }
    }

    /// Convert to compact storage (4x64). Input must be normalized.
    pub fn to_storage(&self, r: &mut FeStorage) {
        let a = &self.n;
        r.n[0] = a[0] | a[1] << 52;
        r.n[1] = a[1] >> 12 | a[2] << 40;
        r.n[2] = a[2] >> 24 | a[3] << 28;
        r.n[3] = a[3] >> 36 | a[4] << 16;
    }

    /// Convert from compact storage. Output limbs are already in canonical form
    /// (each fits its bit-width), matching libsecp256k1 which marks this as normalized=1, magnitude=1.
    pub fn from_storage(&mut self, a: &FeStorage) {
        self.n[0] = a.n[0] & M;
        self.n[1] = a.n[0] >> 52 | ((a.n[1] << 12) & M);
        self.n[2] = a.n[1] >> 40 | ((a.n[2] << 24) & M);
        self.n[3] = a.n[2] >> 28 | ((a.n[3] << 36) & M);
        self.n[4] = a.n[3] >> 16;
    }
}

/// Compact field storage (4x64). Matches libsecp256k1 fe_storage.
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct FeStorage {
    pub n: [u64; 4],
}

impl FeStorage {
    pub fn cmov(&mut self, a: &FeStorage, flag: Choice) {
        for i in 0..4 {
            self.n[i] = u64::conditional_select(&self.n[i], &a.n[i], flag);
        }
    }
}

impl ConstantTimeEq for FieldElement {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.n[0].ct_eq(&other.n[0])
            & self.n[1].ct_eq(&other.n[1])
            & self.n[2].ct_eq(&other.n[2])
            & self.n[3].ct_eq(&other.n[3])
            & self.n[4].ct_eq(&other.n[4])
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl Eq for FieldElement {}

const M62: u64 = u64::MAX >> 2;

/// Convert 5x52 field element to signed62 (for modinv64).
#[inline(always)]
fn fe_to_signed62(a: &FieldElement) -> Signed62 {
    let a0 = a.n[0];
    let a1 = a.n[1];
    let a2 = a.n[2];
    let a3 = a.n[3];
    let a4 = a.n[4];
    Signed62 {
        v: [
            ((a0 | a1 << 52) & M62) as i64,
            ((a1 >> 10 | a2 << 42) & M62) as i64,
            ((a2 >> 20 | a3 << 32) & M62) as i64,
            ((a3 >> 30 | a4 << 22) & M62) as i64,
            (a4 >> 40) as i64,
        ],
    }
}

/// Convert signed62 to 5x52 field element.
#[inline(always)]
fn fe_from_signed62(r: &mut FieldElement, a: &Signed62) {
    let a0 = a.v[0] as u64;
    let a1 = a.v[1] as u64;
    let a2 = a.v[2] as u64;
    let a3 = a.v[3] as u64;
    let a4 = a.v[4] as u64;
    r.n[0] = a0 & M;
    r.n[1] = (a0 >> 52 | a1 << 10) & M;
    r.n[2] = (a1 >> 42 | a2 << 20) & M;
    r.n[3] = (a2 >> 32 | a3 << 30) & M;
    r.n[4] = a3 >> 22 | a4 << 40;
}

/// 5x52 field multiplication using u128 (from field_5x52_int128_impl.h)
#[inline(always)]
fn fe_mul_inner(r: &mut [u64; 5], a: &[u64; 5], b: &[u64; 5]) {
    let a0 = a[0];
    let a1 = a[1];
    let a2 = a[2];
    let a3 = a[3];
    let a4 = a[4];

    let mut d: u128 = (a0 as u128) * (b[3] as u128);
    d += (a1 as u128) * (b[2] as u128);
    d += (a2 as u128) * (b[1] as u128);
    d += (a3 as u128) * (b[0] as u128);

    let mut c: u128 = (a4 as u128) * (b[4] as u128);
    d += (R as u128) * ((c as u64) as u128);
    c >>= 64;

    let t3 = (d as u64) & M;
    d >>= 52;

    d += (a0 as u128) * (b[4] as u128);
    d += (a1 as u128) * (b[3] as u128);
    d += (a2 as u128) * (b[2] as u128);
    d += (a3 as u128) * (b[1] as u128);
    d += (a4 as u128) * (b[0] as u128);
    d += ((R as u128) << 12) * ((c as u64) as u128);

    let mut t4 = (d as u64) & M;
    d >>= 52;
    let tx = t4 >> 48;
    t4 &= M >> 4;

    c = (a0 as u128) * (b[0] as u128);
    d += (a1 as u128) * (b[4] as u128);
    d += (a2 as u128) * (b[3] as u128);
    d += (a3 as u128) * (b[2] as u128);
    d += (a4 as u128) * (b[1] as u128);

    let mut u0 = (d as u64) & M;
    d >>= 52;
    u0 = (u0 << 4) | tx;
    c += (u0 as u128) * ((R >> 4) as u128);

    r[0] = (c as u64) & M;
    c >>= 52;

    c += (a0 as u128) * (b[1] as u128);
    c += (a1 as u128) * (b[0] as u128);
    d += (a2 as u128) * (b[4] as u128);
    d += (a3 as u128) * (b[3] as u128);
    d += (a4 as u128) * (b[2] as u128);
    c += (R as u128) * (((d as u64) & M) as u128);
    d >>= 52;

    r[1] = (c as u64) & M;
    c >>= 52;

    c += (a0 as u128) * (b[2] as u128);
    c += (a1 as u128) * (b[1] as u128);
    c += (a2 as u128) * (b[0] as u128);
    d += (a3 as u128) * (b[4] as u128);
    d += (a4 as u128) * (b[3] as u128);
    c += (R as u128) * ((d as u64) as u128);
    d >>= 64;

    r[2] = (c as u64) & M;
    c >>= 52;
    c += ((R as u128) << 12) * ((d as u64) as u128);
    c += t3 as u128;

    r[3] = (c as u64) & M;
    c >>= 52;
    r[4] = (c as u64) + t4;
}

/// 5x52 field squaring using u128
#[inline(always)]
fn fe_sqr_inner(r: &mut [u64; 5], a: &[u64; 5]) {
    let a0 = a[0];
    let a1 = a[1];
    let a2 = a[2];
    let a3 = a[3];
    let mut a4 = a[4];

    let mut d: u128 = (a0 as u128) * 2 * (a3 as u128);
    d += (a1 as u128) * 2 * (a2 as u128);

    let mut c: u128 = (a4 as u128) * (a4 as u128);
    d += (R as u128) * ((c as u64) as u128);
    c >>= 64;

    let t3 = (d as u64) & M;
    d >>= 52;

    a4 *= 2;
    d += (a0 as u128) * (a4 as u128);
    d += (a1 as u128) * 2 * (a3 as u128);
    d += (a2 as u128) * (a2 as u128);
    d += ((R as u128) << 12) * ((c as u64) as u128);

    let mut t4 = (d as u64) & M;
    d >>= 52;
    let tx = t4 >> 48;
    t4 &= M >> 4;

    c = (a0 as u128) * (a0 as u128);
    d += (a1 as u128) * (a4 as u128);
    d += (a2 as u128) * 2 * (a3 as u128);

    let mut u0 = (d as u64) & M;
    d >>= 52;
    u0 = (u0 << 4) | tx;
    c += (u0 as u128) * ((R >> 4) as u128);

    r[0] = (c as u64) & M;
    c >>= 52;

    let a0_2 = a0 * 2;
    c += (a0_2 as u128) * (a1 as u128);
    d += (a2 as u128) * (a4 as u128);
    d += (a3 as u128) * (a3 as u128);
    c += (R as u128) * (((d as u64) & M) as u128);
    d >>= 52;

    r[1] = (c as u64) & M;
    c >>= 52;

    c += (a0_2 as u128) * (a2 as u128);
    c += (a1 as u128) * (a1 as u128);
    d += (a3 as u128) * (a4 as u128);
    c += (R as u128) * ((d as u64) as u128);
    d >>= 64;

    r[2] = (c as u64) & M;
    c >>= 52;
    c += ((R as u128) << 12) * ((d as u64) as u128);
    c += t3 as u128;

    r[3] = (c as u64) & M;
    c >>= 52;
    r[4] = (c as u64) + t4;
}
