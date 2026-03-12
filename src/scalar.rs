//! Scalar arithmetic modulo the secp256k1 group order n.
//!
//! n = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

#[cfg(target_arch = "x86_64")]
mod scalar_asm {
    use super::Scalar;

    extern "C" {
        /// libsecp256k1 scalar_mul_512: l8 = a * b (512-bit product).
        /// SysV: rdi=l8, rsi=a, rdx=b.
        fn blvm_secp256k1_scalar_mul_512(l8: *mut u64, a: *const Scalar, b: *const Scalar);

        /// libsecp256k1 scalar_reduce_512: reduce 512-bit l mod n into r.
        /// Returns overflow for final reduction. SysV: rdi=r, rsi=l.
        fn blvm_secp256k1_scalar_reduce_512(r: *mut Scalar, l: *const u64) -> u64;
    }

    #[inline(always)]
    pub(super) unsafe fn scalar_mul_512_asm(l: *mut u64, a: *const Scalar, b: *const Scalar) {
        blvm_secp256k1_scalar_mul_512(l, a, b);
    }

    #[inline(always)]
    pub(super) unsafe fn scalar_reduce_512_asm(r: *mut Scalar, l: *const u64) -> u64 {
        blvm_secp256k1_scalar_reduce_512(r, l)
    }
}

use num_bigint::BigUint;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

/// Scalar modulo group order n. 4x64 limb layout (x86_64, aarch64).
#[repr(C)]
#[derive(Clone, Copy, Debug)]
pub struct Scalar {
    pub d: [u64; 4],
}

// secp256k1 group order n
const N_0: u64 = 0xBFD25E8CD0364141;
const N_1: u64 = 0xBAAEDCE6AF48A03B;
const N_2: u64 = 0xFFFFFFFFFFFFFFFE;
const N_3: u64 = 0xFFFFFFFFFFFFFFFF;

// 2^256 - n (for reduction)
const N_C_0: u64 = 0x402DA1732FC9BEBF;
const N_C_1: u64 = 0x4551231950B75FC4;
const N_C_2: u64 = 1;

// n/2 (for is_high)
const N_H_0: u64 = 0xDFE92F46681B20A0;
const N_H_1: u64 = 0x5D576E7357A4501D;
const N_H_2: u64 = 0xFFFFFFFFFFFFFFFF;
const N_H_3: u64 = 0x7FFFFFFFFFFFFFFF;

// modulus n (for safegcd inv)
#[allow(dead_code)]
const N: Scalar = Scalar {
    d: [N_0, N_1, N_2, N_3],
};

const LAMBDA: Scalar = Scalar {
    d: [
        0xDF02967C1B23BD72,
        0x122E22EA20816678,
        0xA5261C028812645A,
        0x5363AD4CC05C30E0,
    ],
};

impl Scalar {
    pub fn zero() -> Self {
        Self { d: [0, 0, 0, 0] }
    }

    pub fn one() -> Self {
        Self { d: [1, 0, 0, 0] }
    }

    pub fn set_int(&mut self, v: u32) {
        self.d[0] = v as u64;
        self.d[1] = 0;
        self.d[2] = 0;
        self.d[3] = 0;
    }

    /// Set from 32-byte big-endian. Reduces mod n.
    pub fn set_b32(&mut self, bin: &[u8; 32]) -> bool {
        self.d[0] = read_be64(&bin[24..32]);
        self.d[1] = read_be64(&bin[16..24]);
        self.d[2] = read_be64(&bin[8..16]);
        self.d[3] = read_be64(&bin[0..8]);
        let overflow = self.check_overflow();
        self.reduce(overflow as u64);
        overflow
    }

    pub fn get_b32(&self, bin: &mut [u8; 32]) {
        write_be64(&mut bin[0..8], self.d[3]);
        write_be64(&mut bin[8..16], self.d[2]);
        write_be64(&mut bin[16..24], self.d[1]);
        write_be64(&mut bin[24..32], self.d[0]);
    }

    fn check_overflow(&self) -> bool {
        let mut yes = 0u64;
        let mut no = 0u64;
        no |= (self.d[3] < N_3) as u64;
        no |= (self.d[2] < N_2) as u64;
        yes |= (self.d[2] > N_2) as u64 & !no;
        no |= (self.d[1] < N_1) as u64;
        yes |= (self.d[1] > N_1) as u64 & !no;
        yes |= (self.d[0] >= N_0) as u64 & !no;
        yes != 0
    }

    fn reduce(&mut self, overflow: u64) {
        let mut t: u128 = self.d[0] as u128 + (overflow as u128 * N_C_0 as u128);
        self.d[0] = t as u64;
        t >>= 64;
        t += self.d[1] as u128 + (overflow as u128 * N_C_1 as u128);
        self.d[1] = t as u64;
        t >>= 64;
        t += self.d[2] as u128 + (overflow as u128 * N_C_2 as u128);
        self.d[2] = t as u64;
        t >>= 64;
        t += self.d[3] as u128;
        self.d[3] = t as u64;
    }

    pub fn is_zero(&self) -> bool {
        (self.d[0] | self.d[1] | self.d[2] | self.d[3]) == 0
    }

    pub fn is_one(&self) -> bool {
        (self.d[0] ^ 1) | self.d[1] | self.d[2] | self.d[3] == 0
    }

    /// True if scalar is odd (d[0] & 1).
    #[allow(dead_code)]
    fn is_odd(&self) -> bool {
        self.d[0] & 1 != 0
    }

    /// True if scalar is even. Used by wnaf_fixed (Pippenger).
    pub(crate) fn is_even(&self) -> bool {
        self.d[0] & 1 == 0
    }

    /// self = a - b (mod n). Result in [0, n-1].
    #[allow(dead_code)]
    fn sub(&mut self, a: &Scalar, b: &Scalar) {
        let mut neg_b = Scalar::zero();
        neg_b.negate(b);
        self.add(a, &neg_b);
    }

    /// self = self / 2. Only valid when self is even.
    #[allow(dead_code)]
    fn half(&mut self) {
        self.d[0] = (self.d[0] >> 1) | (self.d[1] << 63);
        self.d[1] = (self.d[1] >> 1) | (self.d[2] << 63);
        self.d[2] = (self.d[2] >> 1) | (self.d[3] << 63);
        self.d[3] >>= 1;
    }

    /// self = (self + n) / 2. Only valid when self is odd. Result in [0, n-1].
    #[allow(dead_code)]
    fn half_add_n(&mut self) {
        let mut t: u128 = self.d[0] as u128 + N_0 as u128;
        let c0 = t as u64;
        let mut c1 = (t >> 64) as u64;
        t = self.d[1] as u128 + N_1 as u128 + c1 as u128;
        c1 = t as u64;
        let mut c2 = (t >> 64) as u64;
        t = self.d[2] as u128 + N_2 as u128 + c2 as u128;
        c2 = t as u64;
        let mut c3 = (t >> 64) as u64;
        t = self.d[3] as u128 + N_3 as u128 + c3 as u128;
        c3 = t as u64;
        let c4 = (t >> 64) as u64;
        self.d[0] = (c0 >> 1) | (c1 << 63);
        self.d[1] = (c1 >> 1) | (c2 << 63);
        self.d[2] = (c2 >> 1) | (c3 << 63);
        self.d[3] = (c3 >> 1) | (c4 << 63);
        self.reduce(self.check_overflow() as u64);
    }

    /// div2(M, x): x/2 mod n when x even, (x+n)/2 mod n when x odd.
    #[allow(dead_code)]
    fn div2(&mut self) {
        if self.is_odd() {
            self.half_add_n();
        } else {
            self.half();
        }
    }

    /// tmp = a + b (full 257-bit add, no reduction). Used when both are odd and we need (a+b)/2.
    #[allow(dead_code)]
    fn add_no_reduce(a: &Scalar, b: &Scalar) -> [u64; 5] {
        let mut t: u128 = a.d[0] as u128 + b.d[0] as u128;
        let c0 = t as u64;
        let mut c1 = (t >> 64) as u64;
        t = a.d[1] as u128 + b.d[1] as u128 + c1 as u128;
        c1 = t as u64;
        let mut c2 = (t >> 64) as u64;
        t = a.d[2] as u128 + b.d[2] as u128 + c2 as u128;
        c2 = t as u64;
        let mut c3 = (t >> 64) as u64;
        t = a.d[3] as u128 + b.d[3] as u128 + c3 as u128;
        c3 = t as u64;
        let c4 = (t >> 64) as u64;
        [c0, c1, c2, c3, c4]
    }

    /// self = (c0..c4) >> 1, then reduce mod n.
    #[allow(dead_code)]
    fn set_from_5limb_half(&mut self, c: &[u64; 5]) {
        self.d[0] = (c[0] >> 1) | (c[1] << 63);
        self.d[1] = (c[1] >> 1) | (c[2] << 63);
        self.d[2] = (c[2] >> 1) | (c[3] << 63);
        self.d[3] = (c[3] >> 1) | (c[4] << 63);
        self.reduce(self.check_overflow() as u64);
    }

    /// self = (a - b) >> 1 mod n. a and b odd. When a>=b, a-b is even; when a<b, a-b+n is odd, div2 adds n.
    #[allow(dead_code)]
    fn sub_half(&mut self, a: &Scalar, b: &Scalar) {
        self.sub(a, b);
        self.div2();
    }

    pub fn add(&mut self, a: &Scalar, b: &Scalar) -> bool {
        let mut t: u128 = a.d[0] as u128 + b.d[0] as u128;
        self.d[0] = t as u64;
        t >>= 64;
        t += a.d[1] as u128 + b.d[1] as u128;
        self.d[1] = t as u64;
        t >>= 64;
        t += a.d[2] as u128 + b.d[2] as u128;
        self.d[2] = t as u64;
        t >>= 64;
        t += a.d[3] as u128 + b.d[3] as u128;
        self.d[3] = t as u64;
        t >>= 64;
        let overflow = t as u64 + self.check_overflow() as u64;
        debug_assert!(overflow <= 1);
        self.reduce(overflow);
        overflow != 0
    }

    pub fn negate(&mut self, a: &Scalar) {
        let nonzero = if a.is_zero() { 0u64 } else { u64::MAX };
        let mut t: u128 = (!a.d[0]) as u128 + (N_0 + 1) as u128;
        self.d[0] = (t as u64) & nonzero;
        t >>= 64;
        t += (!a.d[1]) as u128 + N_1 as u128;
        self.d[1] = (t as u64) & nonzero;
        t >>= 64;
        t += (!a.d[2]) as u128 + N_2 as u128;
        self.d[2] = (t as u64) & nonzero;
        t >>= 64;
        t += (!a.d[3]) as u128 + N_3 as u128;
        self.d[3] = (t as u64) & nonzero;
    }

    pub fn mul(&mut self, a: &Scalar, b: &Scalar) {
        let mut l = [0u64; 8];
        scalar_mul_512(&mut l, a, b);
        scalar_reduce_512(self, &l);
    }

    /// split_lambda: find r1, r2 such that r1 + r2*lambda == k (mod n)
    pub fn split_lambda(r1: &mut Scalar, r2: &mut Scalar, k: &Scalar) {
        const MINUS_B1: Scalar = Scalar {
            d: [
                (0x6F547FA9u64 << 32) | 0x0ABFE4C3,
                (0xE4437ED6u64 << 32) | 0x010E8828,
                0,
                0,
            ],
        };
        const MINUS_B2: Scalar = Scalar {
            d: [
                (0xD765CDA8u64 << 32) | 0x3DB1562C,
                (0x8A280AC5u64 << 32) | 0x0774346D,
                (0xFFFFFFFFu64 << 32) | 0xFFFFFFFE,
                (0xFFFFFFFFu64 << 32) | 0xFFFFFFFF,
            ],
        };
        const G1: Scalar = Scalar {
            d: [
                (0xE893209Au64 << 32) | 0x45DBB031,
                (0x3DAA8A14u64 << 32) | 0x71E8CA7F,
                (0xE86C90E4u64 << 32) | 0x9284EB15,
                (0x3086D221u64 << 32) | 0xA7D46BCD,
            ],
        };
        const G2: Scalar = Scalar {
            d: [
                (0x1571B4AEu64 << 32) | 0x8AC47F71,
                (0x221208ACu64 << 32) | 0x9DF506C6,
                (0x6F547FA9u64 << 32) | 0x0ABFE4C4,
                (0xE4437ED6u64 << 32) | 0x010E8828,
            ],
        };

        let mut c1 = Scalar::zero();
        let mut c2 = Scalar::zero();
        scalar_mul_shift_var(&mut c1, k, &G1, 384);
        scalar_mul_shift_var(&mut c2, k, &G2, 384);
        let mut t = Scalar::zero();
        t.mul(&c1, &MINUS_B1);
        c1 = t;
        t.mul(&c2, &MINUS_B2);
        c2 = t;
        r2.add(&c1, &c2);
        r1.mul(r2, &LAMBDA);
        let mut neg = Scalar::zero();
        neg.negate(r1);
        r1.add(&neg, k);
    }

    /// Extract `count` bits at `offset` (0..256). Count in [1,32].
    /// For single-limb (offset and offset+count-1 in same u64) uses fast path.
    pub fn get_bits_limb32(&self, offset: u32, count: u32) -> u32 {
        debug_assert!(count > 0 && count <= 32);
        debug_assert!((offset + count - 1) >> 6 == offset >> 6);
        let limb = offset >> 6;
        let shift = offset & 0x3F;
        let mask = if count == 32 {
            u32::MAX
        } else {
            (1u32 << count) - 1
        };
        ((self.d[limb as usize] >> shift) as u32) & mask
    }

    /// Extract `count` bits at `offset`. Count in [1,32], offset+count <= 256.
    pub fn get_bits_var(&self, offset: u32, count: u32) -> u32 {
        debug_assert!(count > 0 && count <= 32);
        debug_assert!(offset + count <= 256);
        if (offset + count - 1) >> 6 == offset >> 6 {
            self.get_bits_limb32(offset, count)
        } else {
            let limb = (offset >> 6) as usize;
            let shift = offset & 0x3F;
            let mask = if count == 32 {
                u32::MAX
            } else {
                (1u32 << count) - 1
            };
            let lo = self.d[limb] >> shift;
            let hi = self.d[limb + 1].wrapping_shl(64u32 - shift);
            ((lo | hi) as u32) & mask
        }
    }

    pub fn split_128(r1: &mut Scalar, r2: &mut Scalar, k: &Scalar) {
        r1.d[0] = k.d[0];
        r1.d[1] = k.d[1];
        r1.d[2] = 0;
        r1.d[3] = 0;
        r2.d[0] = k.d[2];
        r2.d[1] = k.d[3];
        r2.d[2] = 0;
        r2.d[3] = 0;
    }

    /// Modular inverse. Variable-time. r = a^(-1) mod n. If a is zero, r is zero.
    /// On x86_64/aarch64: modinv64 (safegcd). Else: Fermat via num-bigint modpow.
    pub fn inv_var(&mut self, a: &Scalar) {
        if a.is_zero() {
            *self = Scalar::zero();
            return;
        }
        #[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
        {
            use crate::modinv64::{modinv64, SECP256K1_SCALAR_MODINV_MODINFO};
            let mut x = scalar_to_signed62(a);
            modinv64(&mut x, &SECP256K1_SCALAR_MODINV_MODINFO);
            scalar_from_signed62(self, &x);
        }
        #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
        {
            let a_big = scalar_to_biguint(a);
            let n_big = scalar_to_biguint(&N);
            let exp = &n_big - 2u32;
            let inv_big = a_big.modpow(&exp, &n_big);
            biguint_to_scalar(self, &inv_big);
        }
    }

    /// True if scalar is in the upper half [n/2, n).
    pub fn is_high(&self) -> bool {
        let mut yes = 0u64;
        let mut no = 0u64;
        no |= (self.d[3] < N_H_3) as u64;
        yes |= (self.d[3] > N_H_3) as u64 & !no;
        no |= (self.d[2] < N_H_2) as u64 & !yes;
        no |= (self.d[1] < N_H_1) as u64 & !yes;
        yes |= (self.d[1] > N_H_1) as u64 & !no;
        yes |= (self.d[0] > N_H_0) as u64 & !no;
        yes != 0
    }

    /// Conditionally negate: if flag != 0, negate in place. Returns 1 if negated, -1 if not.
    pub fn cond_negate(&mut self, flag: i32) -> i32 {
        let mask = if flag != 0 { u64::MAX } else { 0 };
        let nonzero = if self.is_zero() { 0 } else { u64::MAX };
        let mut t: u128 = (self.d[0] ^ mask) as u128;
        t += ((N_0 + 1) & mask) as u128;
        self.d[0] = (t as u64) & nonzero;
        t >>= 64;
        t += (self.d[1] ^ mask) as u128;
        t += (N_1 & mask) as u128;
        self.d[1] = (t as u64) & nonzero;
        t >>= 64;
        t += (self.d[2] ^ mask) as u128;
        t += (N_2 & mask) as u128;
        self.d[2] = (t as u64) & nonzero;
        t >>= 64;
        t += (self.d[3] ^ mask) as u128;
        t += (N_3 & mask) as u128;
        self.d[3] = (t as u64) & nonzero;
        if mask == 0 {
            -1
        } else {
            1
        }
    }
}

impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.d[0].ct_eq(&other.d[0])
            & self.d[1].ct_eq(&other.d[1])
            & self.d[2].ct_eq(&other.d[2])
            & self.d[3].ct_eq(&other.d[3])
    }
}

impl ConditionallySelectable for Scalar {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Self {
            d: [
                u64::conditional_select(&a.d[0], &b.d[0], choice),
                u64::conditional_select(&a.d[1], &b.d[1], choice),
                u64::conditional_select(&a.d[2], &b.d[2], choice),
                u64::conditional_select(&a.d[3], &b.d[3], choice),
            ],
        }
    }
}

/// Pack 4×64 scalar limbs into 5×62 signed limbs for modinv64.
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
fn scalar_to_signed62(a: &Scalar) -> crate::modinv64::Signed62 {
    const M62: u64 = u64::MAX >> 2;
    let d = &a.d;
    crate::modinv64::Signed62 {
        v: [
            (d[0] & M62) as i64,
            ((d[0] >> 62 | d[1] << 2) & M62) as i64,
            ((d[1] >> 60 | d[2] << 4) & M62) as i64,
            ((d[2] >> 58 | d[3] << 6) & M62) as i64,
            (d[3] >> 56) as i64,
        ],
    }
}

/// Unpack 5×62 signed limbs back to 4×64 scalar limbs.
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
fn scalar_from_signed62(r: &mut Scalar, a: &crate::modinv64::Signed62) {
    let v = &a.v;
    r.d[0] = (v[0] as u64) | ((v[1] as u64) << 62);
    r.d[1] = ((v[1] as u64) >> 2) | ((v[2] as u64) << 60);
    r.d[2] = ((v[2] as u64) >> 4) | ((v[3] as u64) << 58);
    r.d[3] = ((v[3] as u64) >> 6) | ((v[4] as u64) << 56);
}

#[allow(dead_code)]
fn scalar_to_biguint(s: &Scalar) -> BigUint {
    let mut bytes = [0u8; 32];
    s.get_b32(&mut bytes);
    BigUint::from_bytes_be(&bytes)
}

#[allow(dead_code)]
fn biguint_to_scalar(r: &mut Scalar, b: &BigUint) {
    let bytes = b.to_bytes_be();
    let mut buf = [0u8; 32];
    let len = bytes.len().min(32);
    let start = 32 - len;
    buf[start..].copy_from_slice(&bytes[..len]);
    r.set_b32(&buf);
}

fn read_be64(b: &[u8]) -> u64 {
    u64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
}

fn write_be64(b: &mut [u8], v: u64) {
    b[0..8].copy_from_slice(&v.to_be_bytes());
}

fn scalar_mul_512(l: &mut [u64; 8], a: &Scalar, b: &Scalar) {
    #[cfg(target_arch = "x86_64")]
    {
        unsafe {
            scalar_asm::scalar_mul_512_asm(l.as_mut_ptr(), a, b);
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        scalar_mul_512_rust(l, a, b);
    }
}

#[cfg(not(target_arch = "x86_64"))]
fn scalar_mul_512_rust(l: &mut [u64; 8], a: &Scalar, b: &Scalar) {
    let mut c0: u64 = 0;
    let mut c1: u64 = 0;
    let mut c2: u32 = 0;

    macro_rules! muladd_fast {
        ($a:expr, $b:expr) => {{
            let prod = ($a as u128) * ($b as u128);
            let prod_lo = prod as u64;
            let prod_hi = (prod >> 64) as u64;
            let (lo, o) = c0.overflowing_add(prod_lo);
            c0 = lo;
            c1 += prod_hi + o as u64; // C: th = prod_hi + (c0 < tl)
        }};
    }
    macro_rules! muladd {
        ($a:expr, $b:expr) => {{
            let prod = ($a as u128) * ($b as u128);
            let hi = (prod >> 64) as u64;
            let (lo, o1) = c0.overflowing_add(prod as u64);
            c0 = lo;
            let th = hi + o1 as u64;
            let (mid, o2) = c1.overflowing_add(th);
            c1 = mid;
            c2 += o2 as u32;
        }};
    }
    macro_rules! sumadd {
        ($a:expr) => {{
            let (lo, o) = c0.overflowing_add($a);
            c0 = lo;
            c1 += o as u64;
            c2 += (c1 == 0 && o) as u32;
        }};
    }
    macro_rules! extract {
        () => {{
            let n = c0;
            c0 = c1;
            c1 = c2 as u64;
            c2 = 0;
            n
        }};
    }
    macro_rules! extract_fast {
        () => {{
            let n = c0;
            c0 = c1;
            c1 = 0;
            n
        }};
    }

    muladd_fast!(a.d[0], b.d[0]);
    l[0] = extract_fast!();
    muladd!(a.d[0], b.d[1]);
    muladd!(a.d[1], b.d[0]);
    l[1] = extract!();
    muladd!(a.d[0], b.d[2]);
    muladd!(a.d[1], b.d[1]);
    muladd!(a.d[2], b.d[0]);
    l[2] = extract!();
    muladd!(a.d[0], b.d[3]);
    muladd!(a.d[1], b.d[2]);
    muladd!(a.d[2], b.d[1]);
    muladd!(a.d[3], b.d[0]);
    l[3] = extract!();
    muladd!(a.d[1], b.d[3]);
    muladd!(a.d[2], b.d[2]);
    muladd!(a.d[3], b.d[1]);
    l[4] = extract!();
    muladd!(a.d[2], b.d[3]);
    muladd!(a.d[3], b.d[2]);
    l[5] = extract!();
    muladd_fast!(a.d[3], b.d[3]);
    l[6] = extract_fast!();
    l[7] = c0;
}

#[allow(dead_code)]
fn limbs_512_to_biguint(l: &[u64; 8]) -> BigUint {
    let mut acc = BigUint::from(0u64);
    for (i, &limb) in l.iter().enumerate() {
        acc += BigUint::from(limb) << (64 * i);
    }
    acc
}

/// Limb-based 512→256 reduction mod n. Replaces BigUint for hot path.
/// Port of libsecp256k1 scalar_reduce_512 C fallback (muladd/extract).
#[cfg(not(target_arch = "x86_64"))]
fn scalar_reduce_512_limbs(r: &mut Scalar, l: &[u64; 8]) {
    let n0 = l[4];
    let n1 = l[5];
    let n2 = l[6];
    let n3 = l[7];

    let mut c0: u64 = l[0];
    let mut c1: u64 = 0;
    let mut c2: u32 = 0;

    macro_rules! muladd_fast {
        ($a:expr, $b:expr) => {{
            let prod = ($a as u128) * ($b as u128);
            let (lo, o) = c0.overflowing_add(prod as u64);
            c0 = lo;
            c1 += (prod >> 64) as u64 + o as u64;
        }};
    }
    macro_rules! muladd {
        ($a:expr, $b:expr) => {{
            let prod = ($a as u128) * ($b as u128);
            let (lo, o1) = c0.overflowing_add(prod as u64);
            c0 = lo;
            let th = (prod >> 64) as u64 + o1 as u64;
            let (mid, o2) = c1.overflowing_add(th);
            c1 = mid;
            c2 += o2 as u32;
        }};
    }
    macro_rules! sumadd_fast {
        ($a:expr) => {{
            let (lo, o) = c0.overflowing_add($a);
            c0 = lo;
            c1 += o as u64;
        }};
    }
    macro_rules! sumadd {
        ($a:expr) => {{
            let (lo, o) = c0.overflowing_add($a);
            c0 = lo;
            let (mid, o2) = c1.overflowing_add(o as u64);
            c1 = mid;
            c2 += o2 as u32;
        }};
    }
    macro_rules! extract {
        () => {{
            let n = c0;
            c0 = c1;
            c1 = c2 as u64;
            c2 = 0;
            n
        }};
    }
    macro_rules! extract_fast {
        () => {{
            let n = c0;
            c0 = c1;
            c1 = 0;
            n
        }};
    }

    // Reduce 512 bits into 385: m[0..6] = l[0..3] + n[0..3] * N_C
    muladd_fast!(n0, N_C_0);
    let m0 = extract_fast!();
    sumadd_fast!(l[1]);
    muladd!(n1, N_C_0);
    muladd!(n0, N_C_1);
    let m1 = extract!();
    sumadd!(l[2]);
    muladd!(n2, N_C_0);
    muladd!(n1, N_C_1);
    sumadd!(n0);
    let m2 = extract!();
    sumadd!(l[3]);
    muladd!(n3, N_C_0);
    muladd!(n2, N_C_1);
    sumadd!(n1);
    let m3 = extract!();
    muladd!(n3, N_C_1);
    sumadd!(n2);
    let m4 = extract!();
    sumadd_fast!(n3);
    let m5 = extract_fast!();
    let m6 = c0 as u32;

    // Reduce 385 into 258: p[0..4] = m[0..3] + m[4..6] * N_C
    c0 = m0;
    c1 = 0;
    c2 = 0;
    muladd_fast!(m4, N_C_0);
    let p0 = extract_fast!();
    sumadd_fast!(m1);
    muladd!(m5, N_C_0);
    muladd!(m4, N_C_1);
    let p1 = extract!();
    sumadd!(m2);
    muladd!(m6 as u64, N_C_0);
    muladd!(m5, N_C_1);
    sumadd!(m4);
    let p2 = extract!();
    sumadd_fast!(m3);
    muladd_fast!(m6 as u64, N_C_1);
    sumadd_fast!(m5);
    let p3 = extract_fast!();
    let p4 = (c0 + m6 as u64) as u32;

    // Reduce 258 into 256: r = p[0..3] + p4 * N_C
    let mut t: u128 = p0 as u128;
    t += (N_C_0 as u128) * (p4 as u128);
    r.d[0] = t as u64;
    t >>= 64;
    t += p1 as u128;
    t += (N_C_1 as u128) * (p4 as u128);
    r.d[1] = t as u64;
    t >>= 64;
    t += p2 as u128;
    t += p4 as u128;
    r.d[2] = t as u64;
    t >>= 64;
    t += p3 as u128;
    r.d[3] = t as u64;
    let c = (t >> 64) as u64;

    // Final reduction
    scalar_reduce(r, c + scalar_check_overflow(r));
}

/// Add overflow*N_C to r. overflow is 0 or 1.
fn scalar_reduce(r: &mut Scalar, overflow: u64) {
    let of = overflow as u128;
    let mut t: u128 = r.d[0] as u128;
    t += of * (N_C_0 as u128);
    r.d[0] = t as u64;
    t >>= 64;
    t += r.d[1] as u128;
    t += of * (N_C_1 as u128);
    r.d[1] = t as u64;
    t >>= 64;
    t += r.d[2] as u128;
    t += of * (N_C_2 as u128);
    r.d[2] = t as u64;
    t >>= 64;
    r.d[3] = (t as u64).wrapping_add(r.d[3]);
}

/// Returns 1 if r >= n, else 0.
fn scalar_check_overflow(r: &Scalar) -> u64 {
    let mut yes = 0u64;
    let mut no = 0u64;
    no |= (r.d[3] < N_3) as u64;
    no |= (r.d[2] < N_2) as u64;
    yes |= (r.d[2] > N_2) as u64 & !no;
    no |= (r.d[1] < N_1) as u64;
    yes |= (r.d[1] > N_1) as u64 & !no;
    yes |= (r.d[0] >= N_0) as u64 & !no;
    yes
}

fn scalar_reduce_512(r: &mut Scalar, l: &[u64; 8]) {
    #[cfg(target_arch = "x86_64")]
    {
        let c = unsafe { scalar_asm::scalar_reduce_512_asm(r, l.as_ptr()) };
        scalar_reduce(r, c + scalar_check_overflow(r));
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        scalar_reduce_512_limbs(r, l);
    }
}

#[cfg(test)]
#[test]
fn test_scalar_reduce_n_plus_1() {
    let l = [N_0 + 1, N_1, N_2, N_3, 0, 0, 0, 0];
    let mut r = Scalar::zero();
    scalar_reduce_512(&mut r, &l);
    assert!(r.is_one(), "(n+1) mod n = 1, got r.d = {:?}", r.d);
}

#[cfg(test)]
#[test]
fn test_scalar_mul_inv2_times_2() {
    let inv2_hex = "7fffffffffffffffffffffffffffffff5d576e7357a4501ddfe92f46681b20a1";
    let inv2_bytes = hex::decode(inv2_hex).unwrap();
    let mut buf = [0u8; 32];
    buf.copy_from_slice(&inv2_bytes);
    let mut inv2 = Scalar::zero();
    inv2.set_b32(&buf);
    let mut two = Scalar::zero();
    two.set_int(2);
    let mut l = [0u64; 8];
    scalar_mul_512(&mut l, &inv2, &two);
    let mut r = Scalar::zero();
    scalar_reduce_512(&mut r, &l);
    assert!(r.is_one(), "inv2*2 mod n = 1");
}

fn scalar_mul_shift_var(r: &mut Scalar, a: &Scalar, b: &Scalar, shift: u32) {
    assert!(shift >= 256);
    let mut l = [0u64; 8];
    scalar_mul_512(&mut l, a, b);
    let shiftlimbs = (shift >> 6) as usize;
    let shiftlow = shift & 0x3F;
    let shifthigh = 64 - shiftlow;
    r.d[0] = if shift < 512 {
        (l[shiftlimbs] >> shiftlow)
            | (if shift < 448 && shiftlow != 0 {
                l[1 + shiftlimbs] << shifthigh
            } else {
                0
            })
    } else {
        0
    };
    r.d[1] = if shift < 448 {
        (l[1 + shiftlimbs] >> shiftlow)
            | (if shift < 384 && shiftlow != 0 {
                l[2 + shiftlimbs] << shifthigh
            } else {
                0
            })
    } else {
        0
    };
    r.d[2] = if shift < 384 {
        (l[2 + shiftlimbs] >> shiftlow)
            | (if shift < 320 && shiftlow != 0 {
                l[3 + shiftlimbs] << shifthigh
            } else {
                0
            })
    } else {
        0
    };
    r.d[3] = if shift < 320 {
        l[3 + shiftlimbs] >> shiftlow
    } else {
        0
    };
    let bit = (l[(shift - 1) as usize >> 6] >> ((shift - 1) & 0x3F)) & 1;
    scalar_cadd_bit(r, 0, bit != 0);
}

fn scalar_cadd_bit(r: &mut Scalar, bit: u32, flag: bool) {
    let bit = if flag { bit } else { bit + 256 };
    if bit >= 256 {
        return;
    }
    let mut t: u128 = r.d[0] as u128
        + if (bit >> 6) == 0 {
            1u128 << (bit & 0x3F)
        } else {
            0
        };
    r.d[0] = t as u64;
    t >>= 64;
    t += r.d[1] as u128
        + if (bit >> 6) == 1 {
            1u128 << (bit & 0x3F)
        } else {
            0
        };
    r.d[1] = t as u64;
    t >>= 64;
    t += r.d[2] as u128
        + if (bit >> 6) == 2 {
            1u128 << (bit & 0x3F)
        } else {
            0
        };
    r.d[2] = t as u64;
    t >>= 64;
    t += r.d[3] as u128
        + if (bit >> 6) == 3 {
            1u128 << (bit & 0x3F)
        } else {
            0
        };
    r.d[3] = t as u64;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_lambda_identity() {
        // r1 + lambda*r2 == k (mod n)
        let mut k = Scalar::zero();
        k.set_int(42);

        let mut r1 = Scalar::zero();
        let mut r2 = Scalar::zero();
        Scalar::split_lambda(&mut r1, &mut r2, &k);

        let mut lambda_r2 = Scalar::zero();
        lambda_r2.mul(&r2, &LAMBDA);
        let mut check = Scalar::zero();
        check.add(&r1, &lambda_r2);
        assert!(bool::from(check.ct_eq(&k)), "r1 + lambda*r2 should equal k");
    }

    #[test]
    fn test_split_lambda_neg_three() {
        let mut three = Scalar::zero();
        three.set_int(3);
        let mut k = Scalar::zero();
        k.negate(&three); // k = -3 mod n

        let mut r1 = Scalar::zero();
        let mut r2 = Scalar::zero();
        Scalar::split_lambda(&mut r1, &mut r2, &k);

        let mut lambda_r2 = Scalar::zero();
        lambda_r2.mul(&r2, &LAMBDA);
        let mut check = Scalar::zero();
        check.add(&r1, &lambda_r2);
        assert!(
            bool::from(check.ct_eq(&k)),
            "r1 + lambda*r2 should equal k for k=-3"
        );
    }

    #[test]
    fn test_split_lambda_ecdsa_scalar() {
        let mut k = Scalar::zero();
        k.d = [
            11125243483441707226,
            2149109665766520832,
            14302025600096445326,
            4162584031737161978,
        ];

        let n_big = scalar_to_biguint(&N);

        let mut r1 = Scalar::zero();
        let mut r2 = Scalar::zero();
        Scalar::split_lambda(&mut r1, &mut r2, &k);

        let r1_big = scalar_to_biguint(&r1);
        let r2_big = scalar_to_biguint(&r2);
        let n_half = &n_big / BigUint::from(2u64);
        let r1_abs = if r1_big > n_half {
            &n_big - &r1_big
        } else {
            r1_big.clone()
        };
        let r2_abs = if r2_big > n_half {
            &n_big - &r2_big
        } else {
            r2_big.clone()
        };
        assert!(
            r1_abs.bits() <= 128,
            "|r1| exceeds 128 bits: {}",
            r1_abs.bits()
        );
        assert!(
            r2_abs.bits() <= 128,
            "|r2| exceeds 128 bits: {}",
            r2_abs.bits()
        );

        let mut lambda_r2 = Scalar::zero();
        lambda_r2.mul(&r2, &LAMBDA);
        let mut check = Scalar::zero();
        check.add(&r1, &lambda_r2);
        assert!(bool::from(check.ct_eq(&k)), "r1 + lambda*r2 should equal k");
    }

    #[test]
    fn test_split_128_identity() {
        // r1 + 2^128*r2 == k
        let mut k = Scalar::zero();
        k.set_int(0x1234_5678);

        let mut r1 = Scalar::zero();
        let mut r2 = Scalar::zero();
        Scalar::split_128(&mut r1, &mut r2, &k);

        let mut two_128 = Scalar::zero();
        two_128.d[2] = 1;
        let mut r2_shifted = Scalar::zero();
        r2_shifted.mul(&r2, &two_128);
        let mut check = Scalar::zero();
        check.add(&r1, &r2_shifted);
        assert!(bool::from(check.ct_eq(&k)), "r1 + 2^128*r2 should equal k");
    }
}
