//! Elliptic curve multiplication: R = na*A + ng*G
//!
//! Strauss WNAF with endomorphism optimization.
//! Precomputed tables for G and 2^128*G.
//!
//! **Timing:** [`ecmult`] and [`ecmult_gen`] are optimized for public-data workloads (variable
//! WNAF, table lookup). For secret scalars, use the crate’s [`ecmult_const`](crate::ecmult_const)
//! and [`ecmult_gen_const`](crate::ecmult_gen_const) entry points instead.

use crate::field::FieldElement;
use crate::group::{ge_table_set_globalz, generator_g, Ge, GeStorage, Gej};
use crate::scalar::Scalar;

pub(crate) const WINDOW_A: i32 = 5;
const WINDOW_G: i32 = 15;

fn ecmult_table_size(w: i32) -> usize {
    1 << (w - 2)
}

const TABLE_A: usize = 8; // ECMULT_TABLE_SIZE(WINDOW_A)
#[allow(dead_code)]
const TABLE_G: usize = 8192; // ECMULT_TABLE_SIZE(WINDOW_G)

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
mod precomputed {
    include!("ecmult_precomputed.rs");
}

#[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
mod precomputed {
    use crate::field::FeStorage;
    use crate::group::GeStorage;
    pub const PRE_G: [GeStorage; 1] = [GeStorage {
        x: FeStorage { n: [0, 0, 0, 0] },
        y: FeStorage { n: [0, 0, 0, 0] },
    }];
    pub const PRE_G_128: [GeStorage; 1] = [GeStorage {
        x: FeStorage { n: [0, 0, 0, 0] },
        y: FeStorage { n: [0, 0, 0, 0] },
    }];
}

use precomputed::{PRE_G, PRE_G_128};

/// Fill pre_a with odd multiples [1*a, 3*a, ..., (2*n-1)*a]. zr gets z-ratios, z gets final z.
#[inline(always)]
pub(crate) fn ecmult_odd_multiples_table(
    n: usize,
    pre_a: &mut [Ge],
    zr: &mut [FieldElement],
    z: &mut FieldElement,
    a: &Gej,
) {
    debug_assert!(!a.infinity);
    let mut d = Gej::default();
    d.double_var(a);
    let mut d_ge = Ge::default();
    d_ge.set_xy(&d.x, &d.y);
    pre_a[0].set_gej_zinv(a, &d.z);
    let mut ai = Gej::default();
    ai.set_ge(&pre_a[0]);
    ai.z = a.z;
    zr[0] = d.z;
    for i in 1..n {
        let mut sum = Gej::default();
        sum.add_ge_var_rzr(&ai, &d_ge, Some(&mut zr[i]));
        ai = sum;
        pre_a[i].set_xy(&ai.x, &ai.y);
    }
    z.mul(&ai.z, &d.z);
}

/// WNAF: represent scalar as sum of 2^i * wnaf[i]. Returns last_set_bit + 1.
pub(crate) fn ecmult_wnaf(wnaf: &mut [i32], len: usize, a: &Scalar, w: i32) -> usize {
    let mut s = *a;
    let mut sign: i32 = 1;
    if s.get_bits_limb32(255, 1) != 0 {
        let s_in = s;
        s.negate(&s_in);
        sign = -1;
    }
    for w in wnaf.iter_mut().take(len) {
        *w = 0;
    }
    let mut bit: usize = 0;
    let mut last_set_bit: i32 = -1;
    let mut carry: i32 = 0;
    while bit < len {
        if s.get_bits_limb32(bit as u32, 1) as i32 == carry {
            bit += 1;
            continue;
        }
        let now = (w as usize).min(len - bit);
        let word = s.get_bits_var(bit as u32, now as u32) as i32 + carry;
        carry = (word >> (w - 1)) & 1;
        let wval = word - (carry << w);
        wnaf[bit] = sign * wval;
        last_set_bit = bit as i32;
        bit += now;
    }
    (last_set_bit + 1) as usize
}

#[inline(always)]
fn table_get_ge(r: &mut Ge, pre: &[Ge], n: i32, _w: i32) {
    debug_assert!(n != 0 && (n & 1) == 1);
    if n > 0 {
        *r = pre[(n - 1) as usize / 2];
    } else {
        *r = pre[(-n - 1) as usize / 2];
        r.y.normalize_weak();
        let y = r.y;
        r.y.negate(&y, 1);
    }
}

#[inline(always)]
fn table_get_ge_lambda(r: &mut Ge, pre: &[Ge], x_beta: &[FieldElement], n: i32, _w: i32) {
    debug_assert!(n != 0 && (n & 1) == 1);
    if n > 0 {
        let idx = (n - 1) as usize / 2;
        r.set_xy(&x_beta[idx], &pre[idx].y);
    } else {
        let idx = (-n - 1) as usize / 2;
        r.set_xy(&x_beta[idx], &pre[idx].y);
        r.y.normalize_weak();
        let y = r.y;
        r.y.negate(&y, 1);
    }
}

#[inline(always)]
fn table_get_ge_storage(r: &mut Ge, pre: &[GeStorage], n: i32, _w: i32) {
    debug_assert!(n != 0 && (n & 1) == 1);
    if n > 0 {
        r.from_storage(&pre[(n - 1) as usize / 2]);
    } else {
        r.from_storage(&pre[(-n - 1) as usize / 2]);
        r.y.normalize_weak();
        let y = r.y;
        r.y.negate(&y, 1);
    }
}

/// Inner Strauss WNAF. No workarounds. For testing/debugging.
#[cfg(test)]
pub(crate) fn ecmult_inner(r: &mut Gej, a: &Gej, na: &Scalar, ng: Option<&Scalar>) {
    ecmult_strauss(r, a, na, ng);
}

/// R = na*A + ng*G. Double multiply with endomorphism optimization.
#[inline(always)]
pub fn ecmult(r: &mut Gej, a: &Gej, na: &Scalar, ng: Option<&Scalar>) {
    ecmult_strauss(r, a, na, ng);
}

#[inline(always)]
fn ecmult_strauss(r: &mut Gej, a: &Gej, na: &Scalar, ng: Option<&Scalar>) {
    let mut aux = [FieldElement::zero(); TABLE_A];
    let mut pre_a = [
        Ge::default(),
        Ge::default(),
        Ge::default(),
        Ge::default(),
        Ge::default(),
        Ge::default(),
        Ge::default(),
        Ge::default(),
    ];
    let mut wnaf_na_1 = [0i32; 129];
    let mut wnaf_na_lam = [0i32; 129];
    let mut wnaf_ng_1 = [0i32; 129];
    let mut wnaf_ng_128 = [0i32; 129];

    let mut na_1 = Scalar::zero();
    let mut na_lam = Scalar::zero();
    Scalar::split_lambda(&mut na_1, &mut na_lam, na);
    let bits_na_1 = ecmult_wnaf(&mut wnaf_na_1, 129, &na_1, WINDOW_A);
    let bits_na_lam = ecmult_wnaf(&mut wnaf_na_lam, 129, &na_lam, WINDOW_A);
    let mut bits = bits_na_1.max(bits_na_lam);

    let mut ng_1 = Scalar::zero();
    let mut ng_128 = Scalar::zero();
    let mut bits_ng_1 = 0;
    let mut bits_ng_128 = 0;
    if let Some(ng) = ng {
        Scalar::split_128(&mut ng_1, &mut ng_128, ng);
        bits_ng_1 = ecmult_wnaf(&mut wnaf_ng_1, 129, &ng_1, WINDOW_G);
        bits_ng_128 = ecmult_wnaf(&mut wnaf_ng_128, 129, &ng_128, WINDOW_G);
        bits = bits.max(bits_ng_1).max(bits_ng_128);
    }

    let mut z = FieldElement::one();
    if !na.is_zero() && !a.is_infinity() {
        ecmult_odd_multiples_table(TABLE_A, &mut pre_a, &mut aux, &mut z, a);
        ge_table_set_globalz(TABLE_A, &mut pre_a, &aux);
        let beta = crate::group::const_beta();
        for i in 0..TABLE_A {
            aux[i].mul(&pre_a[i].x, &beta);
        }
    }

    r.set_infinity();
    let mut tmpa = Ge::default();
    for i in (0..bits).rev() {
        let r_in = *r;
        r.double_var(&r_in);
        if !na.is_zero() && !a.is_infinity() {
            let n1 = wnaf_na_1[i];
            if n1 != 0 {
                table_get_ge(&mut tmpa, &pre_a, n1, WINDOW_A);
                let r_in = *r;
                r.add_ge_var(&r_in, &tmpa);
            }
            let nlam = wnaf_na_lam[i];
            if nlam != 0 {
                table_get_ge_lambda(&mut tmpa, &pre_a, &aux, nlam, WINDOW_A);
                let r_in = *r;
                r.add_ge_var(&r_in, &tmpa);
            }
        }
        if ng.is_some() {
            let n1 = if i < bits_ng_1 { wnaf_ng_1[i] } else { 0 };
            if n1 != 0 {
                table_get_ge_storage(&mut tmpa, &PRE_G, n1, WINDOW_G);
                let r_in = *r;
                r.add_zinv_var(&r_in, &tmpa, &z);
            }
            let n128 = if i < bits_ng_128 { wnaf_ng_128[i] } else { 0 };
            if n128 != 0 {
                table_get_ge_storage(&mut tmpa, &PRE_G_128, n128, WINDOW_G);
                let r_in = *r;
                r.add_zinv_var(&r_in, &tmpa, &z);
            }
        }
    }
    if !r.is_infinity() {
        let r_z = r.z;
        r.z.mul(&r_z, &z);
    }
}

/// R = ng*G (scalar multiplication by generator only).
/// Dedicated path: skips split_lambda and wnaf for na when na=0.
#[inline(always)]
pub fn ecmult_gen(r: &mut Gej, ng: &Scalar) {
    ecmult_gen_strauss(r, ng);
}

/// Strauss WNAF for ng*G only. No na-path overhead.
fn ecmult_gen_strauss(r: &mut Gej, ng: &Scalar) {
    let mut wnaf_ng_1 = [0i32; 129];
    let mut wnaf_ng_128 = [0i32; 129];

    let mut ng_1 = Scalar::zero();
    let mut ng_128 = Scalar::zero();
    Scalar::split_128(&mut ng_1, &mut ng_128, ng);
    let bits_ng_1 = ecmult_wnaf(&mut wnaf_ng_1, 129, &ng_1, WINDOW_G);
    let bits_ng_128 = ecmult_wnaf(&mut wnaf_ng_128, 129, &ng_128, WINDOW_G);
    let bits = bits_ng_1.max(bits_ng_128);

    let z = FieldElement::one();
    r.set_infinity();
    let mut tmpa = Ge::default();
    for i in (0..bits).rev() {
        let r_in = *r;
        r.double_var(&r_in);
        let n1 = if i < bits_ng_1 { wnaf_ng_1[i] } else { 0 };
        if n1 != 0 {
            table_get_ge_storage(&mut tmpa, &PRE_G, n1, WINDOW_G);
            let r_in = *r;
            r.add_zinv_var(&r_in, &tmpa, &z);
        }
        let n128 = if i < bits_ng_128 { wnaf_ng_128[i] } else { 0 };
        if n128 != 0 {
            table_get_ge_storage(&mut tmpa, &PRE_G_128, n128, WINDOW_G);
            let r_in = *r;
            r.add_zinv_var(&r_in, &tmpa, &z);
        }
    }
}

const WNAF_BITS: usize = 128;
const PIPPENGER_MAX_BUCKET_WINDOW: i32 = 12;
const ECMULT_PIPPENGER_THRESHOLD: usize = 88;

fn wnaf_size(w: i32) -> usize {
    WNAF_BITS.div_ceil(w as usize).max(1)
}

/// Fixed-size WNAF for Pippenger. Returns skew (0 or 1).
fn wnaf_fixed(wnaf: &mut [i32], s: &Scalar, w: i32) -> i32 {
    let n_wnaf = wnaf_size(w);
    if s.is_zero() {
        for wv in wnaf.iter_mut().take(n_wnaf) {
            *wv = 0;
        }
        return 0;
    }
    let skew = if s.is_even() { 1 } else { 0 };
    let last_w = WNAF_BITS - (n_wnaf - 1) * (w as usize);

    let mut max_pos = 0;
    for pos in (1..n_wnaf).rev() {
        let bits = if pos == n_wnaf - 1 {
            last_w as u32
        } else {
            w as u32
        };
        let val = s.get_bits_var((pos * w as usize) as u32, bits) as i32;
        if val != 0 {
            max_pos = pos;
            break;
        }
        wnaf[pos] = 0;
    }

    wnaf[0] = s.get_bits_var(0, w as u32) as i32 + skew;

    let mut pos = 1;
    while pos <= max_pos {
        let bits = if pos == n_wnaf - 1 {
            last_w as u32
        } else {
            w as u32
        };
        let val = s.get_bits_var((pos * w as usize) as u32, bits) as i32;
        if (val & 1) == 0 {
            wnaf[pos - 1] -= 1 << w;
            wnaf[pos] = val + 1;
        } else {
            wnaf[pos] = val;
        }
        if pos >= 2
            && ((wnaf[pos - 1] == 1 && wnaf[pos - 2] < 0)
                || (wnaf[pos - 1] == -1 && wnaf[pos - 2] > 0))
        {
            if wnaf[pos - 1] == 1 {
                wnaf[pos - 2] += 1 << w;
            } else {
                wnaf[pos - 2] -= 1 << w;
            }
            wnaf[pos - 1] = 0;
        }
        pos += 1;
    }
    skew
}

fn pippenger_bucket_window(n: usize) -> i32 {
    match n {
        0..=1 => 1,
        2..=4 => 2,
        5..=20 => 3,
        21..=57 => 4,
        58..=136 => 5,
        137..=235 => 6,
        236..=1260 => 7,
        1261..=4420 => 9,
        4421..=7880 => 10,
        7881..=16050 => 11,
        _ => PIPPENGER_MAX_BUCKET_WINDOW,
    }
}

/// Endomorphism split: s1,s2 from split_lambda; p2 = lambda*p1; negate if is_high.
fn ecmult_endo_split(s1: &mut Scalar, s2: &mut Scalar, p1: &mut Ge, p2: &mut Ge) {
    let tmp = *s1;
    Scalar::split_lambda(s1, s2, &tmp);
    p2.mul_lambda(p1);
    if s1.is_high() {
        let s1_val = *s1;
        s1.negate(&s1_val);
        let p1_val = *p1;
        p1.neg(&p1_val);
    }
    if s2.is_high() {
        let s2_val = *s2;
        s2.negate(&s2_val);
        let p2_val = *p2;
        p2.neg(&p2_val);
    }
}

struct PippengerState {
    wnaf_na: Vec<i32>,
    skew_na: Vec<i32>,
    input_pos: Vec<usize>,
}

fn pippenger_wnaf(
    buckets: &mut [Gej],
    bucket_window: i32,
    state: &PippengerState,
    r: &mut Gej,
    points: &[Ge],
) {
    let w = bucket_window + 1;
    let n_wnaf = wnaf_size(w);
    let n_buckets = ecmult_table_size(bucket_window + 2);
    let no = state.input_pos.len();

    r.set_infinity();
    if no == 0 {
        return;
    }

    let mut tmp = Ge::default();
    for i in (0..n_wnaf).rev() {
        for b in buckets.iter_mut().take(n_buckets) {
            b.set_infinity();
        }

        for np in 0..no {
            let idx = state.input_pos[np];
            let n = state.wnaf_na[np * n_wnaf + i];
            let skew = state.skew_na[np];

            if i == 0 && skew != 0 {
                tmp.neg(&points[idx]);
                let b0 = buckets[0];
                buckets[0].add_ge_var(&b0, &tmp);
            }
            if n > 0 {
                let bi = (n - 1) / 2;
                let b = buckets[bi as usize];
                buckets[bi as usize].add_ge_var(&b, &points[idx]);
            } else if n < 0 {
                let bi = (-(n + 1)) / 2;
                tmp.neg(&points[idx]);
                let b = buckets[bi as usize];
                buckets[bi as usize].add_ge_var(&b, &tmp);
            }
        }

        for _ in 0..bucket_window {
            let r_in = *r;
            r.double_var(&r_in);
        }

        let mut running_sum = Gej::default();
        running_sum.set_infinity();
        for j in (1..n_buckets).rev() {
            let rs_old = running_sum;
            running_sum.add_var(&rs_old, &buckets[j]);
            let r_old = *r;
            r.add_var(&r_old, &running_sum);
        }
        let rs_old = running_sum;
        running_sum.add_var(&rs_old, &buckets[0]);
        let r_in = *r;
        r.double_var(&r_in);
        let r_old = *r;
        r.add_var(&r_old, &running_sum);
    }
}

/// Pippenger multi-multiply: R = g_scalar*G + sum_i(scalar_i * point_i).
/// Faster than ecmult_multi_simple for n >= 88. Uses heap allocation.
/// Uses endomorphism split; split outputs are < 2^128 so 128-bit WNAF suffices.
pub fn ecmult_multi_pippenger(r: &mut Gej, g_scalar: &Scalar, scalars: &[Scalar], points: &[Ge]) {
    let n = scalars.len().min(points.len());
    if n == 0 {
        if g_scalar.is_zero() {
            r.set_infinity();
        } else {
            let mut inf = Gej::default();
            inf.set_infinity();
            let zero = Scalar::zero();
            ecmult(r, &inf, &zero, Some(g_scalar));
        }
        return;
    }

    let bucket_window = pippenger_bucket_window(n);
    let w = bucket_window + 1;
    let n_wnaf = wnaf_size(w);
    let n_buckets = ecmult_table_size(bucket_window + 2);
    let entries = 2 * n + 2;

    let mut pipp_scalars = vec![Scalar::zero(); entries];
    let mut pipp_points = vec![Ge::default(); entries];
    let mut buckets = vec![Gej::default(); n_buckets];
    let mut wnaf_na = vec![0i32; entries * n_wnaf];
    let mut skew_na = Vec::with_capacity(entries);
    let mut input_pos = Vec::with_capacity(entries);

    let mut idx = 0;
    if !g_scalar.is_zero() {
        pipp_scalars[0] = *g_scalar;
        pipp_points[0] = generator_g();
        let (s01, s1) = pipp_scalars.split_at_mut(1);
        let (p01, p1) = pipp_points.split_at_mut(1);
        ecmult_endo_split(&mut s01[0], &mut s1[0], &mut p01[0], &mut p1[0]);
        idx = 2;
    }

    let mut point_idx = 0;
    while point_idx < n {
        pipp_scalars[idx] = scalars[point_idx];
        pipp_points[idx] = points[point_idx];
        let (s_lo, s_hi) = pipp_scalars.split_at_mut(idx + 1);
        let (p_lo, p_hi) = pipp_points.split_at_mut(idx + 1);
        ecmult_endo_split(&mut s_lo[idx], &mut s_hi[0], &mut p_lo[idx], &mut p_hi[0]);
        idx += 2;
        point_idx += 1;
    }

    let mut no = 0;
    for np in 0..idx {
        if pipp_scalars[np].is_zero() || pipp_points[np].is_infinity() {
            continue;
        }
        let skew = wnaf_fixed(
            &mut wnaf_na[no * n_wnaf..(no + 1) * n_wnaf],
            &pipp_scalars[np],
            w,
        );
        skew_na.push(skew);
        input_pos.push(np);
        no += 1;
    }

    let state = PippengerState {
        wnaf_na,
        skew_na,
        input_pos,
    };
    pippenger_wnaf(&mut buckets, bucket_window, &state, r, &pipp_points);
}

/// Pippenger without endomorphism: for scalars that may not fit 128 bits after split.
/// Uses 256-bit WNAF. Fallback when endo path fails (e.g. Schnorr batch).
fn ecmult_multi_pippenger_no_endo(
    r: &mut Gej,
    g_scalar: &Scalar,
    scalars: &[Scalar],
    points: &[Ge],
) {
    const WNAF_BITS_256: usize = 256;
    let n = scalars.len().min(points.len());
    if n == 0 {
        if g_scalar.is_zero() {
            r.set_infinity();
        } else {
            ecmult_gen(r, g_scalar);
        }
        return;
    }

    fn wnaf_size_256(w: i32) -> usize {
        WNAF_BITS_256.div_ceil(w as usize).max(1)
    }

    let bucket_window = pippenger_bucket_window(n);
    let w = bucket_window + 1;
    let n_wnaf = wnaf_size_256(w);
    let n_buckets = ecmult_table_size(bucket_window + 2);
    let entries = n + 1; // G + n points, no split

    let mut pipp_scalars = vec![Scalar::zero(); entries];
    let mut pipp_points = vec![Ge::default(); entries];
    let mut buckets = vec![Gej::default(); n_buckets];
    let mut wnaf_na = vec![0i32; entries * n_wnaf];
    let mut skew_na = Vec::with_capacity(entries);
    let mut input_pos = Vec::with_capacity(entries);

    let mut idx = 0;
    if !g_scalar.is_zero() {
        pipp_scalars[0] = *g_scalar;
        pipp_points[0] = generator_g();
        idx = 1;
    }

    for i in 0..n {
        pipp_scalars[idx] = scalars[i];
        pipp_points[idx] = points[i];
        idx += 1;
    }

    let last_w = WNAF_BITS_256 - (n_wnaf - 1) * (w as usize);
    let mut no = 0;
    for np in 0..idx {
        if pipp_scalars[np].is_zero() || pipp_points[np].is_infinity() {
            continue;
        }
        let skew = wnaf_fixed_256(
            &mut wnaf_na[no * n_wnaf..(no + 1) * n_wnaf],
            &pipp_scalars[np],
            w,
            n_wnaf,
            last_w,
        );
        skew_na.push(skew);
        input_pos.push(np);
        no += 1;
    }

    let state = PippengerState {
        wnaf_na,
        skew_na,
        input_pos,
    };
    pippenger_wnaf_256(
        &mut buckets,
        bucket_window,
        w,
        n_wnaf,
        &state,
        r,
        &pipp_points,
    );
}

fn wnaf_fixed_256(wnaf: &mut [i32], s: &Scalar, w: i32, n_wnaf: usize, last_w: usize) -> i32 {
    if s.is_zero() {
        for wv in wnaf.iter_mut().take(n_wnaf) {
            *wv = 0;
        }
        return 0;
    }
    let skew = if s.is_even() { 1 } else { 0 };

    let mut max_pos = 0;
    for pos in (1..n_wnaf).rev() {
        let bits = if pos == n_wnaf - 1 {
            last_w as u32
        } else {
            w as u32
        };
        let val = s.get_bits_var((pos * w as usize) as u32, bits) as i32;
        if val != 0 {
            max_pos = pos;
            break;
        }
        wnaf[pos] = 0;
    }

    wnaf[0] = s.get_bits_var(0, w as u32) as i32 + skew;

    let mut pos = 1;
    while pos <= max_pos {
        let bits = if pos == n_wnaf - 1 {
            last_w as u32
        } else {
            w as u32
        };
        let val = s.get_bits_var((pos * w as usize) as u32, bits) as i32;
        if (val & 1) == 0 {
            wnaf[pos - 1] -= 1 << w;
            wnaf[pos] = val + 1;
        } else {
            wnaf[pos] = val;
        }
        if pos >= 2
            && ((wnaf[pos - 1] == 1 && wnaf[pos - 2] < 0)
                || (wnaf[pos - 1] == -1 && wnaf[pos - 2] > 0))
        {
            if wnaf[pos - 1] == 1 {
                wnaf[pos - 2] += 1 << w;
            } else {
                wnaf[pos - 2] -= 1 << w;
            }
            wnaf[pos - 1] = 0;
        }
        pos += 1;
    }
    skew
}

fn pippenger_wnaf_256(
    buckets: &mut [Gej],
    bucket_window: i32,
    _w: i32,
    n_wnaf: usize,
    state: &PippengerState,
    r: &mut Gej,
    points: &[Ge],
) {
    let n_buckets = ecmult_table_size(bucket_window + 2);
    let no = state.input_pos.len();

    r.set_infinity();
    if no == 0 {
        return;
    }

    let mut tmp = Ge::default();
    for i in (0..n_wnaf).rev() {
        for b in buckets.iter_mut().take(n_buckets) {
            b.set_infinity();
        }

        for np in 0..no {
            let idx = state.input_pos[np];
            let n = state.wnaf_na[np * n_wnaf + i];
            let skew = state.skew_na[np];

            if i == 0 && skew != 0 {
                tmp.neg(&points[idx]);
                let b0 = buckets[0];
                buckets[0].add_ge_var(&b0, &tmp);
            }
            if n > 0 {
                let bi = (n - 1) / 2;
                let b = buckets[bi as usize];
                buckets[bi as usize].add_ge_var(&b, &points[idx]);
            } else if n < 0 {
                let bi = (-(n + 1)) / 2;
                tmp.neg(&points[idx]);
                let b = buckets[bi as usize];
                buckets[bi as usize].add_ge_var(&b, &tmp);
            }
        }

        for _ in 0..bucket_window {
            let r_in = *r;
            r.double_var(&r_in);
        }

        let mut running_sum = Gej::default();
        running_sum.set_infinity();
        for j in (1..n_buckets).rev() {
            let rs_old = running_sum;
            running_sum.add_var(&rs_old, &buckets[j]);
            let r_old = *r;
            r.add_var(&r_old, &running_sum);
        }
        let rs_old = running_sum;
        running_sum.add_var(&rs_old, &buckets[0]);
        let r_in = *r;
        r.double_var(&r_in);
        let r_old = *r;
        r.add_var(&r_old, &running_sum);
    }
}

/// Multi-multiply: R = g_scalar*G + sum_i(scalar_i * point_i).
/// Uses Pippenger for n >= 88, ecmult_multi_simple otherwise.
/// Uses no-endo Pippenger (256-bit WNAF) for correctness with full scalars.
pub fn ecmult_multi(r: &mut Gej, g_scalar: &Scalar, scalars: &[Scalar], points: &[Ge]) {
    let n = scalars.len().min(points.len());
    if n >= ECMULT_PIPPENGER_THRESHOLD {
        ecmult_multi_pippenger_no_endo(r, g_scalar, scalars, points);
    } else {
        ecmult_multi_simple(r, g_scalar, scalars, points);
    }
}

/// Multi-multiply (simple algorithm, no scratch): R = g_scalar*G + sum_i(scalar_i * point_i).
/// Used when scratch is NULL in libsecp256k1. Correct for any n; for n >= 88 consider Pippenger.
pub fn ecmult_multi_simple(r: &mut Gej, g_scalar: &Scalar, scalars: &[Scalar], points: &[Ge]) {
    debug_assert_eq!(scalars.len(), points.len());
    let n = scalars.len().min(points.len());

    let mut inf = Gej::default();
    inf.set_infinity();
    let zero = Scalar::zero();

    if n == 0 {
        if g_scalar.is_zero() {
            r.set_infinity();
        } else {
            ecmult(r, &inf, &zero, Some(g_scalar));
        }
        return;
    }

    // r = g_scalar * G
    ecmult(r, &inf, &zero, Some(g_scalar));

    let mut tmp = Gej::default();
    for i in 0..n {
        if scalars[i].is_zero() || points[i].is_infinity() {
            continue;
        }
        let mut pj = Gej::default();
        pj.set_ge(&points[i]);
        ecmult(&mut tmp, &pj, &scalars[i], None);
        let r_old = *r;
        r.add_var(&r_old, &tmp);
    }
}

#[cfg(test)]
fn test_ecmult_multi_simple_two_g_impl() {
    use crate::group::generator_g;
    let g = generator_g();
    let mut one = Scalar::zero();
    one.set_int(1);
    let scalars = [one];
    let points = [g];
    let mut r = Gej::default();
    ecmult_multi_simple(&mut r, &one, &scalars, &points);
    let mut two_g = Gej::default();
    two_g.set_ge(&g);
    let mut sum = Gej::default();
    sum.add_ge_var(&two_g, &g);
    let mut r_aff = Ge::default();
    r_aff.set_gej_var(&r);
    let mut two_g_aff = Ge::default();
    two_g_aff.set_gej_var(&sum);
    r_aff.x.normalize();
    r_aff.y.normalize();
    two_g_aff.x.normalize();
    two_g_aff.y.normalize();
    assert!(
        FieldElement::fe_equal(&r_aff.x, &two_g_aff.x)
            && FieldElement::fe_equal(&r_aff.y, &two_g_aff.y),
        "ecmult_multi_simple(1*G + 1*G) should equal 2*G"
    );
}

/// Reference: R = na*A using double-and-add (no endomorphism). Correct fallback when Strauss fails.
#[allow(dead_code)]
pub(crate) fn ecmult_simple(r: &mut Gej, a: &Gej, na: &Scalar) {
    if a.is_infinity() || na.is_zero() {
        r.set_infinity();
        return;
    }
    // Handle negative: na*a = -(|na|*a)
    if na.get_bits_limb32(255, 1) != 0 {
        let mut na_abs = Scalar::zero();
        na_abs.negate(na);
        ecmult_simple(r, a, &na_abs);
        if !r.is_infinity() {
            let r_in = *r;
            r.neg(&r_in);
        }
        return;
    }
    r.set_infinity();
    let mut tmp = Gej::default();
    for i in (0..256).rev() {
        let r_in = *r;
        r.double_var(&r_in);
        if na.get_bits_limb32(i, 1) != 0 {
            let r_in = *r;
            tmp.add_var(&r_in, a);
            *r = tmp;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group::{generator_g, Ge, Gej};
    use crate::scalar::Scalar;

    /// Debug: print WNAF for na=-3
    #[test]
    fn test_ecmult_debug_neg3() {
        let mut three = Scalar::zero();
        three.set_int(3);
        let mut na = Scalar::zero();
        na.negate(&three); // -3
        let mut na_1 = Scalar::zero();
        let mut na_lam = Scalar::zero();
        Scalar::split_lambda(&mut na_1, &mut na_lam, &na);
        let mut wnaf_1 = [0i32; 129];
        let mut wnaf_lam = [0i32; 129];
        let bits_1 = super::ecmult_wnaf(&mut wnaf_1, 129, &na_1, super::WINDOW_A);
        let bits_lam = super::ecmult_wnaf(&mut wnaf_lam, 129, &na_lam, super::WINDOW_A);
        eprintln!("na=-3: na_1 bits={}, na_lam bits={}", bits_1, bits_lam);
        for i in 0..bits_1.max(bits_lam) {
            if wnaf_1[i] != 0 || wnaf_lam[i] != 0 {
                eprintln!("  i={}: wnaf_1={} wnaf_lam={}", i, wnaf_1[i], wnaf_lam[i]);
            }
        }
    }

    #[test]
    fn test_ecmult_multi_simple_two_g() {
        super::test_ecmult_multi_simple_two_g_impl();
    }

    /// ecmult_multi must match ecmult_multi_simple for small n.
    #[test]
    fn test_ecmult_multi_vs_simple() {
        let g = generator_g();
        let mut one = Scalar::zero();
        one.set_int(1);
        for n in [1, 2, 5, 10] {
            let scalars: Vec<Scalar> = (0..n).map(|_| one).collect();
            let points: Vec<Ge> = (0..n).map(|_| g).collect();
            let mut r_multi = Gej::default();
            super::ecmult_multi(&mut r_multi, &one, &scalars, &points);
            let mut r_simple = Gej::default();
            super::ecmult_multi_simple(&mut r_simple, &one, &scalars, &points);
            let mut aff_multi = Ge::default();
            aff_multi.set_gej_var(&r_multi);
            let mut aff_simple = Ge::default();
            aff_simple.set_gej_var(&r_simple);
            aff_multi.x.normalize();
            aff_multi.y.normalize();
            aff_simple.x.normalize();
            aff_simple.y.normalize();
            assert!(
                FieldElement::fe_equal(&aff_multi.x, &aff_simple.x)
                    && FieldElement::fe_equal(&aff_multi.y, &aff_simple.y),
                "ecmult_multi != ecmult_multi_simple for n={}",
                n
            );
        }
    }

    /// ecmult_multi (Pippenger) must match ecmult_multi_simple for n=88.
    #[test]
    fn test_ecmult_multi_pippenger_vs_simple() {
        let g = generator_g();
        let mut one = Scalar::zero();
        one.set_int(1);
        let n = 88;
        let scalars: Vec<Scalar> = (0..n).map(|_| one).collect();
        let points: Vec<Ge> = (0..n).map(|_| g).collect();
        let mut r_multi = Gej::default();
        super::ecmult_multi(&mut r_multi, &one, &scalars, &points);
        let mut r_simple = Gej::default();
        super::ecmult_multi_simple(&mut r_simple, &one, &scalars, &points);
        let mut aff_multi = Ge::default();
        aff_multi.set_gej_var(&r_multi);
        let mut aff_simple = Ge::default();
        aff_simple.set_gej_var(&r_simple);
        aff_multi.x.normalize();
        aff_multi.y.normalize();
        aff_simple.x.normalize();
        aff_simple.y.normalize();
        assert!(
            FieldElement::fe_equal(&aff_multi.x, &aff_simple.x)
                && FieldElement::fe_equal(&aff_multi.y, &aff_simple.y),
            "ecmult_multi (Pippenger) != ecmult_multi_simple for n=88"
        );
    }

    /// Pippenger vs simple with varied scalars and same point (Schnorr batch-like).
    #[test]
    fn test_ecmult_multi_pippenger_varied_scalars() {
        let g = generator_g();
        let mut g_scalar = Scalar::zero();
        g_scalar.set_int(12345);
        let n = 90;
        let mut scalars: Vec<Scalar> = Vec::with_capacity(n);
        for i in 0..n {
            let mut s = Scalar::zero();
            s.set_int((i as u32).wrapping_mul(0x9E3779B9)); // varied values
            scalars.push(s);
        }
        let points: Vec<Ge> = (0..n).map(|_| g).collect();
        let mut r_multi = Gej::default();
        super::ecmult_multi(&mut r_multi, &g_scalar, &scalars, &points);
        let mut r_simple = Gej::default();
        super::ecmult_multi_simple(&mut r_simple, &g_scalar, &scalars, &points);
        let mut aff_multi = Ge::default();
        aff_multi.set_gej_var(&r_multi);
        let mut aff_simple = Ge::default();
        aff_simple.set_gej_var(&r_simple);
        aff_multi.x.normalize();
        aff_multi.y.normalize();
        aff_simple.x.normalize();
        aff_simple.y.normalize();
        assert!(
            FieldElement::fe_equal(&aff_multi.x, &aff_simple.x)
                && FieldElement::fe_equal(&aff_multi.y, &aff_simple.y),
            "ecmult_multi (Pippenger) != ecmult_multi_simple for n=90 varied scalars"
        );
    }

    /// Pippenger (no-endo) vs simple with full 256-bit scalars.
    #[test]
    fn test_ecmult_multi_pippenger_full_scalars() {
        let g = generator_g();
        let mut g_scalar = Scalar::zero();
        g_scalar.set_b32(&[
            0x12, 0x34, 0x56, 0x78, 0x9a, 0xbc, 0xde, 0xf0, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66,
            0x77, 0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff, 0x00, 0x11, 0x22, 0x33, 0x44,
            0x55, 0x66, 0x77, 0x88,
        ]);
        let n = 90;
        let mut scalars: Vec<Scalar> = Vec::with_capacity(n);
        for i in 0..n {
            let mut s = Scalar::zero();
            let mut bytes = [0u8; 32];
            bytes[0] = (i as u8).wrapping_add(1);
            bytes[31] = (i as u8).wrapping_mul(0x9E);
            s.set_b32(&bytes);
            scalars.push(s);
        }
        let points: Vec<Ge> = (0..n).map(|_| g).collect();
        let mut r_multi = Gej::default();
        super::ecmult_multi(&mut r_multi, &g_scalar, &scalars, &points);
        let mut r_simple = Gej::default();
        super::ecmult_multi_simple(&mut r_simple, &g_scalar, &scalars, &points);
        let mut aff_multi = Ge::default();
        aff_multi.set_gej_var(&r_multi);
        let mut aff_simple = Ge::default();
        aff_simple.set_gej_var(&r_simple);
        aff_multi.x.normalize();
        aff_multi.y.normalize();
        aff_simple.x.normalize();
        aff_simple.y.normalize();
        assert!(
            FieldElement::fe_equal(&aff_multi.x, &aff_simple.x)
                && FieldElement::fe_equal(&aff_multi.y, &aff_simple.y),
            "ecmult_multi (Pippenger no-endo) != ecmult_multi_simple for n=90 full scalars"
        );
    }

    /// Pippenger vs simple with non-G point (pubkey from Schnorr batch).
    #[test]
    fn test_ecmult_multi_pippenger_non_g_point() {
        use crate::schnorr;
        let seckey = [0u8; 32];
        let mut seckey = seckey;
        seckey[31] = 1; // key 0x01
        let pk_bytes = schnorr::xonly_pubkey_from_secret(&seckey).unwrap();
        let pk_ge = schnorr::lift_x(&pk_bytes).unwrap();
        let mut g_scalar = Scalar::zero();
        g_scalar.set_int(1);
        let n = 90;
        let mut scalars: Vec<Scalar> = Vec::with_capacity(n);
        for i in 0..n {
            let mut s = Scalar::zero();
            s.set_int((i as u32).wrapping_add(1)); // 1..90
            scalars.push(s);
        }
        let points: Vec<Ge> = (0..n).map(|_| pk_ge).collect();
        let mut r_multi = Gej::default();
        super::ecmult_multi(&mut r_multi, &g_scalar, &scalars, &points);
        let mut r_simple = Gej::default();
        super::ecmult_multi_simple(&mut r_simple, &g_scalar, &scalars, &points);
        let mut aff_multi = Ge::default();
        aff_multi.set_gej_var(&r_multi);
        let mut aff_simple = Ge::default();
        aff_simple.set_gej_var(&r_simple);
        aff_multi.x.normalize();
        aff_multi.y.normalize();
        aff_simple.x.normalize();
        aff_simple.y.normalize();
        assert!(
            FieldElement::fe_equal(&aff_multi.x, &aff_simple.x)
                && FieldElement::fe_equal(&aff_multi.y, &aff_simple.y),
            "ecmult_multi (Pippenger) != ecmult_multi_simple for n=90 with pubkey point"
        );
    }

    /// Compare Strauss ecmult with double-and-add reference for na*a (no ng).
    #[test]
    fn test_ecmult_pre_a_vs_simple() {
        let g = generator_g();

        // Test cases: (a_mult, na, ng) where a = a_mult*G, we compute na*a + ng*G
        // ng=None means we only compute na*a
        let cases: Vec<(i32, i32, Option<i32>)> = vec![
            (1, 1, None),
            (1, 3, None),
            (1, -3, None), // known failing case (negative na)
            (2, 1, None),
            (2, 3, None),
            (2, -3, None),
            (1, 5, None),
            (1, -5, None),
            (2, 1, Some(1)),             // 1*2G + 1*G = 3G
            (2, 1, Some(2)),             // 1*2G + 2*G = 4G
            (17, 1, Some(1)),            // 1*17G + 1*G = 18G (pubkey=0x11*G case)
            (2, 1, Some(-1)),            // 1*2G + (-1)*G = G (ng with top bit set)
            (2, -1, Some(1)),            // (-1)*2G + 1*G = -G (na with top bit set)
            (2, -1, Some(-1)),           // (-1)*2G + (-1)*G = -3G (both)
            (17, 0x1234, Some(0x5678)),  // Larger scalars: 0x1234*17G + 0x5678*G
            (17, 0x1234, Some(-0x5678)), // ng with top bit set, larger scalar
        ];
        for (a_mult, na_val, ng_val) in cases {
            // a = a_mult * G
            let mut a = Gej::default();
            a.set_ge(&g);
            for _ in 1..a_mult {
                let mut s = Gej::default();
                s.add_ge_var(&a, &g);
                a = s;
            }

            let mut na = Scalar::zero();
            if na_val >= 0 {
                na.set_int(na_val as u32);
            } else {
                let mut abs = Scalar::zero();
                abs.set_int((-na_val) as u32);
                na.negate(&abs);
            }

            let ng_opt = ng_val.map(|v| {
                let mut s = Scalar::zero();
                if v >= 0 {
                    s.set_int(v as u32);
                } else {
                    let mut abs = Scalar::zero();
                    abs.set_int((-v) as u32);
                    s.negate(&abs);
                }
                s
            });

            // Reference: na*a + ng*G
            let mut ref_r = Gej::default();
            ecmult_simple(&mut ref_r, &a, &na);
            if let Some(ref ng) = ng_opt {
                let mut ng_g = Gej::default();
                super::ecmult_gen(&mut ng_g, ng);
                let mut sum = Gej::default();
                sum.add_var(&ref_r, &ng_g);
                ref_r = sum;
            }

            // Strauss
            let mut strauss_r = Gej::default();
            ecmult_inner(&mut strauss_r, &a, &na, ng_opt.as_ref());

            if ref_r.is_infinity() {
                assert!(
                    strauss_r.is_infinity(),
                    "a_mult={} na={} ng={:?}: ref=inf, strauss should be inf",
                    a_mult,
                    na_val,
                    ng_val
                );
                continue;
            }
            let mut ref_aff = Ge::default();
            ref_aff.set_gej_var(&ref_r);
            let mut strauss_aff = Ge::default();
            strauss_aff.set_gej_var(&strauss_r);
            ref_aff.x.normalize();
            ref_aff.y.normalize();
            strauss_aff.x.normalize();
            strauss_aff.y.normalize();
            assert!(
                FieldElement::fe_equal(&ref_aff.x, &strauss_aff.x)
                    && FieldElement::fe_equal(&ref_aff.y, &strauss_aff.y),
                "a_mult={} na={} ng={:?}: ref != strauss. ref=({:?},{:?}) strauss=({:?},{:?})",
                a_mult,
                na_val,
                ng_val,
                ref_aff.x,
                ref_aff.y,
                strauss_aff.x,
                strauss_aff.y
            );
        }
    }

    /// ecmult_simple(G, na) must equal ecmult_gen(na) for any na.
    #[test]
    fn test_ecmult_simple_vs_gen() {
        let g_aff = generator_g();
        let mut g = Gej::default();
        g.set_ge(&g_aff);
        for na_val in [1u32, 0x12345678, 0xFFFFFFFF] {
            let mut na = Scalar::zero();
            na.set_int(na_val);
            let mut simple_r = Gej::default();
            ecmult_simple(&mut simple_r, &g, &na);
            let mut gen_r = Gej::default();
            super::ecmult_gen(&mut gen_r, &na);
            let mut simple_aff = Ge::default();
            simple_aff.set_gej_var(&simple_r);
            let mut gen_aff = Ge::default();
            gen_aff.set_gej_var(&gen_r);
            simple_aff.x.normalize();
            simple_aff.y.normalize();
            gen_aff.x.normalize();
            gen_aff.y.normalize();
            assert!(
                FieldElement::fe_equal(&simple_aff.x, &gen_aff.x)
                    && FieldElement::fe_equal(&simple_aff.y, &gen_aff.y),
                "ecmult_simple(G,na) != ecmult_gen(na) for na={}",
                na_val
            );
        }
    }

    /// ecmult_simple(G, u2) with ECDSA u2 must equal ecmult_gen(u2).
    #[test]
    fn test_ecmult_simple_ecdsa_u2() {
        use crate::ecdsa::ecdsa_sig_sign;
        let seckey_bytes = [0x11u8; 32];
        let msg_bytes = [0x22u8; 32];
        let nonce_bytes = [0x33u8; 32];
        let seckey = {
            let mut s = Scalar::zero();
            s.set_b32(&seckey_bytes);
            s
        };
        let message = {
            let mut s = Scalar::zero();
            s.set_b32(&msg_bytes);
            s
        };
        let nonce = {
            let mut s = Scalar::zero();
            s.set_b32(&nonce_bytes);
            s
        };
        let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
        let mut sn = Scalar::zero();
        sn.inv_var(&sig.1);
        let mut u2 = Scalar::zero();
        u2.mul(&sn, &sig.0);
        let g_aff = generator_g();
        let mut g = Gej::default();
        g.set_ge(&g_aff);
        let mut simple_r = Gej::default();
        ecmult_simple(&mut simple_r, &g, &u2);
        let mut gen_r = Gej::default();
        super::ecmult_gen(&mut gen_r, &u2);
        let mut simple_aff = Ge::default();
        simple_aff.set_gej_var(&simple_r);
        let mut gen_aff = Ge::default();
        gen_aff.set_gej_var(&gen_r);
        simple_aff.x.normalize();
        simple_aff.y.normalize();
        gen_aff.x.normalize();
        gen_aff.y.normalize();
        assert!(
            FieldElement::fe_equal(&simple_aff.x, &gen_aff.x)
                && FieldElement::fe_equal(&simple_aff.y, &gen_aff.y),
            "ecmult_simple(G,u2) != ecmult_gen(u2) for ECDSA u2"
        );
    }

    /// ecmult_simple(17*G, u2) must equal ecmult_gen(u2*17) for ECDSA u2.
    /// Uses 17*G built by repeated add (avoids z=1).
    #[test]
    fn test_ecmult_simple_17g_ecdsa_u2() {
        use crate::ecdsa::ecdsa_sig_sign;
        let seckey_bytes = [0x11u8; 32];
        let msg_bytes = [0x22u8; 32];
        let nonce_bytes = [0x33u8; 32];
        let seckey = {
            let mut s = Scalar::zero();
            s.set_b32(&seckey_bytes);
            s
        };
        let message = {
            let mut s = Scalar::zero();
            s.set_b32(&msg_bytes);
            s
        };
        let nonce = {
            let mut s = Scalar::zero();
            s.set_b32(&nonce_bytes);
            s
        };
        let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
        let mut sn = Scalar::zero();
        sn.inv_var(&sig.1);
        let mut u2 = Scalar::zero();
        u2.mul(&sn, &sig.0);
        // Build 17*G by repeated addition (like test_ecmult_pre_a_vs_simple)
        let g_aff = generator_g();
        let mut pubkeyj = Gej::default();
        pubkeyj.set_ge(&g_aff);
        for _ in 1..17 {
            let mut s = Gej::default();
            s.add_ge_var(&pubkeyj, &g_aff);
            pubkeyj = s;
        }
        let mut simple_r = Gej::default();
        ecmult_simple(&mut simple_r, &pubkeyj, &u2);
        let mut seventeen = Scalar::zero();
        seventeen.set_int(17);
        let mut u2_17 = Scalar::zero();
        u2_17.mul(&u2, &seventeen);
        let mut gen_r = Gej::default();
        super::ecmult_gen(&mut gen_r, &u2_17);
        let mut simple_aff = Ge::default();
        simple_aff.set_gej_var(&simple_r);
        let mut gen_aff = Ge::default();
        gen_aff.set_gej_var(&gen_r);
        simple_aff.x.normalize();
        simple_aff.y.normalize();
        gen_aff.x.normalize();
        gen_aff.y.normalize();
        assert!(
            FieldElement::fe_equal(&simple_aff.x, &gen_aff.x)
                && FieldElement::fe_equal(&simple_aff.y, &gen_aff.y),
            "ecmult_simple(17*G,u2) != ecmult_gen(u2*17) for ECDSA u2"
        );
    }

    /// ecmult_gen(17) must equal 17*G from repeated add.
    #[test]
    fn test_ecmult_gen_17_vs_repeated_add() {
        let mut seventeen = Scalar::zero();
        seventeen.set_int(17);
        let mut gen_17g = Gej::default();
        super::ecmult_gen(&mut gen_17g, &seventeen);
        let g_aff = generator_g();
        let mut ref_17g = Gej::default();
        ref_17g.set_ge(&g_aff);
        for _ in 1..17 {
            let mut s = Gej::default();
            s.add_ge_var(&ref_17g, &g_aff);
            ref_17g = s;
        }
        let mut gen_aff = Ge::default();
        gen_aff.set_gej_var(&gen_17g);
        let mut ref_aff = Ge::default();
        ref_aff.set_gej_var(&ref_17g);
        gen_aff.x.normalize();
        gen_aff.y.normalize();
        ref_aff.x.normalize();
        ref_aff.y.normalize();
        assert!(
            FieldElement::fe_equal(&gen_aff.x, &ref_aff.x)
                && FieldElement::fe_equal(&gen_aff.y, &ref_aff.y),
            "ecmult_gen(17) must equal 17*G from repeated add"
        );
    }

    /// (ecmult_gen(17) -> set_ge -> rescale) must equal 17*G from repeated add.
    #[test]
    fn test_rescale_gen17_matches_repeated_add() {
        let mut seventeen = Scalar::zero();
        seventeen.set_int(17);
        let mut gen_17g = Gej::default();
        super::ecmult_gen(&mut gen_17g, &seventeen);
        let mut ge_aff = Ge::default();
        ge_aff.set_gej_var(&gen_17g);
        let mut pt = Gej::default();
        pt.set_ge(&ge_aff);
        let one = {
            let mut f = FieldElement::zero();
            f.set_int(1);
            f
        };
        assert!(
            FieldElement::fe_equal(&pt.z, &one),
            "pt should have z=1 after set_ge"
        );
        let two = {
            let mut f = FieldElement::zero();
            f.set_int(2);
            f
        };
        pt.rescale(&two);
        let g_aff = generator_g();
        let mut ref_17g = Gej::default();
        ref_17g.set_ge(&g_aff);
        for _ in 1..17 {
            let mut s = Gej::default();
            s.add_ge_var(&ref_17g, &g_aff);
            ref_17g = s;
        }
        let mut pt_aff = Ge::default();
        pt_aff.set_gej_var(&pt);
        let mut ref_aff = Ge::default();
        ref_aff.set_gej_var(&ref_17g);
        pt_aff.x.normalize();
        pt_aff.y.normalize();
        ref_aff.x.normalize();
        ref_aff.y.normalize();
        assert!(
            FieldElement::fe_equal(&pt_aff.x, &ref_aff.x)
                && FieldElement::fe_equal(&pt_aff.y, &ref_aff.y),
            "rescaled ecmult_gen(17) must equal 17*G from repeated add"
        );
    }

    /// ECDSA verify case: pubkey=17*G, u1=z/s, u2=r/s. Compare Strauss vs simple+gen.
    #[test]
    fn test_ecmult_ecdsa_case_vs_simple() {
        use crate::ecdsa::{ecdsa_sig_sign, pubkey_from_secret};
        let seckey_bytes = [0x11u8; 32];
        let msg_bytes = [0x22u8; 32];
        let nonce_bytes = [0x33u8; 32];
        let seckey = {
            let mut s = Scalar::zero();
            s.set_b32(&seckey_bytes);
            s
        };
        let message = {
            let mut s = Scalar::zero();
            s.set_b32(&msg_bytes);
            s
        };
        let nonce = {
            let mut s = Scalar::zero();
            s.set_b32(&nonce_bytes);
            s
        };
        let _pubkey = pubkey_from_secret(&seckey);
        let sig = ecdsa_sig_sign(&seckey, &message, &nonce).expect("sign");
        let mut sn = Scalar::zero();
        sn.inv_var(&sig.1);
        let mut u1 = Scalar::zero();
        u1.mul(&sn, &message);
        let mut u2 = Scalar::zero();
        u2.mul(&sn, &sig.0);

        // Use 17*G from repeated add (avoids z=1). pubkey_from_secret gives same point.
        let g_aff = generator_g();
        let mut pubkeyj = Gej::default();
        pubkeyj.set_ge(&g_aff);
        for _ in 1..17 {
            let mut s = Gej::default();
            s.add_ge_var(&pubkeyj, &g_aff);
            pubkeyj = s;
        }

        // Reference: ecmult_simple(17*G, u2) + ecmult_gen(u1)
        let mut na_part = Gej::default();
        ecmult_simple(&mut na_part, &pubkeyj, &u2);
        let mut ng_part = Gej::default();
        super::ecmult_gen(&mut ng_part, &u1);
        let mut ref_simple = Gej::default();
        ref_simple.add_var(&na_part, &ng_part);

        // Reference 2: (u2*17+u1)*G via ecmult_gen (mathematical truth)
        let mut seventeen = Scalar::zero();
        seventeen.set_int(17);
        let mut u2_17 = Scalar::zero();
        u2_17.mul(&u2, &seventeen);
        let mut expected_scalar = Scalar::zero();
        expected_scalar.add(&u2_17, &u1);
        let mut ref_gen = Gej::default();
        super::ecmult_gen(&mut ref_gen, &expected_scalar);

        // Strauss (ecmult_inner uses same rescale as public ecmult)
        let mut strauss_r = Gej::default();
        ecmult_inner(&mut strauss_r, &pubkeyj, &u2, Some(&u1));

        let mut ref_simple_aff = Ge::default();
        ref_simple_aff.set_gej_var(&ref_simple);
        let mut ref_gen_aff = Ge::default();
        ref_gen_aff.set_gej_var(&ref_gen);
        let mut strauss_aff = Ge::default();
        strauss_aff.set_gej_var(&strauss_r);
        ref_simple_aff.x.normalize();
        ref_simple_aff.y.normalize();
        ref_gen_aff.x.normalize();
        ref_gen_aff.y.normalize();
        strauss_aff.x.normalize();
        strauss_aff.y.normalize();

        // Diagnostic: test na-only and ng-only separately
        let mut na_only_strauss = Gej::default();
        ecmult_inner(&mut na_only_strauss, &pubkeyj, &u2, None);
        let mut na_only_simple = Gej::default();
        ecmult_simple(&mut na_only_simple, &pubkeyj, &u2);
        let mut na_str_aff = Ge::default();
        na_str_aff.set_gej_var(&na_only_strauss);
        let mut na_sim_aff = Ge::default();
        na_sim_aff.set_gej_var(&na_only_simple);
        na_str_aff.x.normalize();
        na_str_aff.y.normalize();
        na_sim_aff.x.normalize();
        na_sim_aff.y.normalize();
        assert!(
            FieldElement::fe_equal(&na_str_aff.x, &na_sim_aff.x)
                && FieldElement::fe_equal(&na_str_aff.y, &na_sim_aff.y),
            "ECDSA case: na-only Strauss != simple. u2.d={:?}",
            u2.d
        );

        // Strauss must match (u2*17+u1)*G
        assert!(
            FieldElement::fe_equal(&ref_gen_aff.x, &strauss_aff.x)
                && FieldElement::fe_equal(&ref_gen_aff.y, &strauss_aff.y),
            "ECDSA case: strauss != (u2*17+u1)*G"
        );
    }
}
