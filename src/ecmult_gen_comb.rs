//! Constant-time `k*G` via signed-digit multi-comb (libsecp256k1 `ecmult_gen` algorithm).
//!
//! This is a port of the algorithm described in
//! `libsecp256k1/src/ecmult_gen_impl.h` (Pieter Wuille, Peter Dettman).
//!
//! Trade-offs vs. the GLV-based `ecmult_const(g, k)` path:
//!   * one-time ~22 KiB precomputed table (built lazily on first call);
//!   * 44 mixed adds + 3 doublings per call, vs. ~52 adds + 125 doublings;
//!   * no per-call odd-multiples / global-z / lambda-table setup;
//!   * still constant-time: cmov scan over the table, unified `add_ge`, `double_ct`.
//!
//! No blinding is applied; the signing path is constant-time but does not randomize
//! intermediate `Z` coordinates between calls. This matches the security model of the
//! rest of the crate (we rely on CT formulas, not randomization).

use std::sync::OnceLock;
use subtle::Choice;

use crate::field::FieldElement;
use crate::group::{ge_set_all_gej_var, generator_g, Ge, GeStorage, Gej};
use crate::scalar::Scalar;

/// Constant-time conditional move on a `GeStorage`, parameterised by a bit-mask `mask`
/// that must be `0u64` (no copy) or `u64::MAX` (copy). Mirrors libsecp's
/// `secp256k1_ge_storage_cmov` formulation `(dst & ~mask) | (src & mask)`.
///
/// We construct `mask` from a `volatile`-loaded flag in the caller to keep the compiler
/// honest about the constant-time intent (matches libsecp's `volatile int vflag = flag`).
#[inline(always)]
fn ge_storage_cmov_mask(dst: &mut GeStorage, src: &GeStorage, mask: u64) {
    let nmask = !mask;
    dst.x.n[0] = (dst.x.n[0] & nmask) | (src.x.n[0] & mask);
    dst.x.n[1] = (dst.x.n[1] & nmask) | (src.x.n[1] & mask);
    dst.x.n[2] = (dst.x.n[2] & nmask) | (src.x.n[2] & mask);
    dst.x.n[3] = (dst.x.n[3] & nmask) | (src.x.n[3] & mask);
    dst.y.n[0] = (dst.y.n[0] & nmask) | (src.y.n[0] & mask);
    dst.y.n[1] = (dst.y.n[1] & nmask) | (src.y.n[1] & mask);
    dst.y.n[2] = (dst.y.n[2] & nmask) | (src.y.n[2] & mask);
    dst.y.n[3] = (dst.y.n[3] & nmask) | (src.y.n[3] & mask);
}

/// Number of independent comb blocks. Each has its own table indexed by `COMB_TEETH` bits.
pub(crate) const COMB_BLOCKS: usize = 11;
/// Number of bits read into a single table lookup.
pub(crate) const COMB_TEETH: usize = 6;
/// Distance (in scalar bit positions) between successive teeth of one block.
pub(crate) const COMB_SPACING: usize = 4;
/// Bits covered by the comb (≥ 256 to cover the scalar range).
pub(crate) const COMB_BITS: usize = COMB_BLOCKS * COMB_TEETH * COMB_SPACING; // 264
/// Stored entries per block: half of `2^COMB_TEETH` (negation symmetry).
pub(crate) const COMB_POINTS: usize = 1 << (COMB_TEETH - 1); // 32
/// Number of u32 limbs needed to hold the recoded scalar (`(COMB_BITS + 31) / 32`).
const RECODED_LIMBS: usize = (COMB_BITS + 31) >> 5; // 9

/// Lazy comb context.
struct CombContext {
    /// `prec_table[block][index]` = table entry for block `block`, packed mask digit `index`.
    /// `index` is the lower `(COMB_TEETH - 1)` bits of the digit; the top bit is encoded
    /// as a sign and resolved by negating `y` after the cmov scan.
    prec_table: [[GeStorage; COMB_POINTS]; COMB_BLOCKS],
    /// `scalar_offset = 1 + (2^COMB_BITS - 1) / 2  (mod n)`.
    scalar_offset: Scalar,
    /// `ge_offset = -G`. Added after the comb walk to cancel the `+1` baked into
    /// `scalar_offset` and the implicit `+(2^COMB_BITS - 1)/2 * G` contributed by the
    /// comb decomposition. See `libsecp256k1` derivation in `ecmult_gen_impl.h`.
    ge_offset: Ge,
}

static COMB_CTX: OnceLock<CombContext> = OnceLock::new();

#[inline]
fn comb_ctx() -> &'static CombContext {
    COMB_CTX.get_or_init(build_comb_context)
}

/// Build the comb context (precomputed table + scalar/ge offsets).
fn build_comb_context() -> CombContext {
    let g = generator_g();

    // ge_offset = -G (unblinded variant of libsecp's blinded ge_offset = b*G).
    let mut ge_offset = Ge::default();
    ge_offset.neg(&g);
    ge_offset.x.normalize();
    ge_offset.y.normalize();

    // Compute scalar_offset = 1 + (2^COMB_BITS - 1)/2  (mod n).
    //   diff = 2^(COMB_BITS-1) - 1/2  ==  (2^COMB_BITS - 1)/2  (mod n).
    // We mirror libsecp's `ecmult_gen_scalar_diff`: doubling-and-add to 2^(BITS-1), then
    // subtract 1/2 mod n.
    let mut diff = Scalar::one();
    for _ in 0..(COMB_BITS - 1) {
        let d = diff;
        diff.add(&d, &d);
    }
    // half = 1/2 mod n  ==  (n+1)/2 because n is odd.
    let mut half = Scalar::zero();
    half.half_modn(&Scalar::one());
    let mut neg_half = Scalar::zero();
    neg_half.negate(&half);
    let diff_in = diff;
    diff.add(&diff_in, &neg_half);
    let mut scalar_offset = Scalar::zero();
    scalar_offset.add(&Scalar::one(), &diff);

    let prec_table = compute_table(&g);

    CombContext {
        prec_table,
        scalar_offset,
        ge_offset,
    }
}

/// Precompute the comb table for generator `gen`.
///
/// Mirrors `secp256k1_ecmult_gen_compute_table` in
/// `libsecp256k1/src/ecmult_gen_compute_table_impl.h`.
fn compute_table(gen: &Ge) -> [[GeStorage; COMB_POINTS]; COMB_BLOCKS] {
    // u starts at gen/2 (multiplied via a simple ladder to avoid relying on ecmult).
    let mut u = Gej::default();
    u.set_infinity();
    let mut half_one = Scalar::zero();
    half_one.half_modn(&Scalar::one());
    for i in (0..256).rev() {
        let u_in = u;
        u.double_var(&u_in);
        if half_one.get_bits_var(i, 1) != 0 {
            let u_in = u;
            u.add_ge_var(&u_in, gen);
        }
    }

    // Build per-block tables in Jacobian, then batch-convert to affine.
    let mut vs: Vec<Gej> = vec![Gej::default(); COMB_BLOCKS * COMB_POINTS];
    let mut ds: [Gej; COMB_TEETH] = [Gej::default(); COMB_TEETH];
    let mut vs_pos = 0usize;

    for block in 0..COMB_BLOCKS {
        // u = 2^(block * COMB_TEETH * COMB_SPACING) * gen/2 at this point.
        let mut sum = Gej::default();
        sum.set_infinity();
        for tooth in 0..COMB_TEETH {
            // sum += 2^((block*COMB_TEETH + tooth) * COMB_SPACING) * gen/2
            let sum_in = sum;
            sum.add_var(&sum_in, &u);
            let u_in = u;
            u.double_var(&u_in);
            ds[tooth] = u;
            // Advance u by `COMB_SPACING - 1` more doublings, except for the very last
            // tooth of the very last block (libsecp condition).
            if block + tooth != COMB_BLOCKS + COMB_TEETH - 2 {
                for _ in 1..COMB_SPACING {
                    let u_in = u;
                    u.double_var(&u_in);
                }
            }
        }

        // First entry (i=0): all teeth contribute -1 → vs[..] = -sum.
        vs[vs_pos].neg(&sum);
        vs_pos += 1;
        // For each new tooth, double the table by adding ds[tooth] to existing entries.
        for tooth in 0..(COMB_TEETH - 1) {
            let stride = 1usize << tooth;
            for _index in 0..stride {
                let prev = vs[vs_pos - stride];
                let d = ds[tooth];
                vs[vs_pos].add_var(&prev, &d);
                vs_pos += 1;
            }
        }
    }
    debug_assert_eq!(vs_pos, COMB_BLOCKS * COMB_POINTS);

    // Batch convert Gej → Ge, then to GeStorage.
    let mut affine: Vec<Ge> = vec![Ge::default(); COMB_BLOCKS * COMB_POINTS];
    ge_set_all_gej_var(&mut affine, &vs);

    // Pack into [BLOCKS][POINTS] storage table.
    let zero_storage = GeStorage {
        x: crate::field::FeStorage { n: [0; 4] },
        y: crate::field::FeStorage { n: [0; 4] },
    };
    let mut table = [[zero_storage; COMB_POINTS]; COMB_BLOCKS];
    for block in 0..COMB_BLOCKS {
        for index in 0..COMB_POINTS {
            let pt = &mut affine[block * COMB_POINTS + index];
            debug_assert!(!pt.is_infinity());
            pt.x.normalize();
            pt.y.normalize();
            pt.x.to_storage(&mut table[block][index].x);
            pt.y.to_storage(&mut table[block][index].y);
        }
    }
    table
}

/// Constant-time `r = k*G` using the multi-comb. Public API.
///
/// On the secret-dependent path:
///   * the cmov scan over `COMB_POINTS` (=32) entries is independent of the secret digit;
///   * `add_ge` is the Brier–Joye unified formula (CT for any pair of valid points);
///   * `double_ct` is the unbranched doubling.
pub fn ecmult_gen_const_comb(r: &mut Gej, gn: &Scalar) {
    let ctx = comb_ctx();
    // Start from infinity so every block uses the same `add_ge` path (Brier–Joye handles
    // a.infinity with cmov). A separate `set_ge` on the first nonzero digit would branch
    // on how many leading comb digits are zero — secret-dependent timing.
    r.set_infinity();

    // d = gn + scalar_offset (mod n). Recoded into u32 limbs covering COMB_BITS bits.
    let mut d = Scalar::zero();
    d.add(gn, &ctx.scalar_offset);
    // Recoded buffer holds 9 u32 limbs (covers COMB_BITS=264 bits). The scalar itself is
    // only 256 bits, so the top limb stays zero — bits beyond 256 are pure padding for
    // the comb decomposition and do not contribute to the result. Mirrors libsecp's
    // `i < 8 && i < ((COMB_BITS + 31) >> 5)` bound in `ecmult_gen_impl.h`.
    let mut recoded = [0u32; RECODED_LIMBS];
    for (i, slot) in recoded.iter_mut().take(8).enumerate() {
        *slot = d.get_bits_limb32(32 * i as u32, 32);
    }

    // Outer loop over comb_off in (COMB_SPACING-1 .. 0). Inner loop over blocks.
    let mut comb_off: i32 = COMB_SPACING as i32 - 1;
    let mut adds = GeStorage {
        x: crate::field::FeStorage { n: [0; 4] },
        y: crate::field::FeStorage { n: [0; 4] },
    };
    let mut neg_y = FieldElement::zero();
    let mut add = Ge::default();

    loop {
        for block in 0..COMB_BLOCKS {
            // Build the COMB_TEETH-bit digit from the recoded scalar.
            let mut bits: u32 = 0;
            let mut bit_pos = (block * COMB_TEETH * COMB_SPACING) as u32 + comb_off as u32;
            for tooth in 0..COMB_TEETH {
                let bit = (recoded[(bit_pos >> 5) as usize] >> (bit_pos & 0x1f)) & 1;
                bits |= bit << tooth;
                bit_pos += COMB_SPACING as u32;
            }

            // Top bit is the sign; remaining (COMB_TEETH - 1) bits are the table index.
            let sign = (bits >> (COMB_TEETH - 1)) & 1;
            let abs = (bits ^ sign.wrapping_neg()) & (COMB_POINTS as u32 - 1);

            // Constant-time table read: scan all entries with cmov; secret-dependent
            // memory access pattern is uniform.
            //
            // Mask construction mirrors libsecp's `mask0 = vflag + ~0; mask1 = ~mask0`
            // pattern, but unrolled into one branchless u64 arithmetic per entry. The
            // `core::hint::black_box` on the equality flag prevents the optimizer from
            // turning the scan into a branch on `abs`.
            let abs_bx = core::hint::black_box(abs);
            for (index, entry) in ctx.prec_table[block].iter().enumerate() {
                let eq = ((index as u32) == abs_bx) as u64;
                let mask = 0u64.wrapping_sub(eq);
                ge_storage_cmov_mask(&mut adds, entry, mask);
            }

            // Promote to affine (storage is already canonical) and conditionally negate y.
            add.from_storage(&adds);
            neg_y.negate(&add.y, 1);
            add.y.cmov(&neg_y, Choice::from(sign as u8));

            // Accumulate into r (always mixed add from Jacobian; infinity + P handled in add_ge).
            let r_in = *r;
            r.add_ge(&r_in, &add);
        }

        if comb_off == 0 {
            break;
        }
        // Double `r` between consecutive comb iterations (CT, no early-out).
        let r_in = *r;
        r.double_ct(&r_in);
        comb_off -= 1;
    }

    // Final correction: + ge_offset (= -G) cancels the +1 baked into scalar_offset.
    let r_in = *r;
    r.add_ge(&r_in, &ctx.ge_offset);
}

/// Eagerly build the comb context. Useful for callers that don't want first-call latency.
pub fn precompute() {
    let _ = comb_ctx();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecmult::ecmult_gen;
    use sha2::{Digest, Sha256};

    fn assert_gej_eq_affine(a: &Gej, b: &Gej, label: &str) {
        let mut ag = Ge::default();
        let mut bg = Ge::default();
        ag.set_gej_var(a);
        bg.set_gej_var(b);
        ag.x.normalize();
        ag.y.normalize();
        bg.x.normalize();
        bg.y.normalize();
        assert_eq!(ag.x, bg.x, "{label}: x");
        assert_eq!(ag.y, bg.y, "{label}: y");
        assert_eq!(ag.infinity, bg.infinity, "{label}: infinity");
    }

    #[test]
    fn comb_matches_ecmult_gen_small() {
        for v in 1u32..50 {
            let mut k = Scalar::zero();
            k.set_int(v);
            let mut want = Gej::default();
            ecmult_gen(&mut want, &k);
            let mut got = Gej::default();
            ecmult_gen_const_comb(&mut got, &k);
            assert_gej_eq_affine(&got, &want, &format!("k={v}"));
        }
    }

    #[test]
    fn comb_matches_ecmult_gen_hashed() {
        for i in 0u32..256 {
            let h = Sha256::digest((i as u64).to_le_bytes());
            let b32: [u8; 32] = h.into();
            let mut k = Scalar::zero();
            let _ = k.set_b32(&b32);
            if k.is_zero() {
                continue;
            }
            let mut want = Gej::default();
            ecmult_gen(&mut want, &k);
            let mut got = Gej::default();
            ecmult_gen_const_comb(&mut got, &k);
            assert_gej_eq_affine(&got, &want, &format!("hash i={i}"));
        }
    }

    #[test]
    fn comb_matches_ecmult_const_g() {
        use crate::ecmult_const::ecmult_const;
        let g = generator_g();
        for i in 0u32..32 {
            let h = Sha256::digest((42u64 + i as u64).to_le_bytes());
            let b32: [u8; 32] = h.into();
            let mut k = Scalar::zero();
            let _ = k.set_b32(&b32);
            if k.is_zero() {
                continue;
            }
            let mut want = Gej::default();
            ecmult_const(&mut want, &g, &k);
            let mut got = Gej::default();
            ecmult_gen_const_comb(&mut got, &k);
            assert_gej_eq_affine(&got, &want, &format!("ecmult_const G k={i}"));
        }
    }
}
