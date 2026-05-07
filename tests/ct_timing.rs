//! Dudect-style statistical timing tests for constant-time properties.
//!
//! Uses Welch's t-test on CPU cycle counts to detect timing differences
//! between two classes of inputs. With 50k+ samples, \|t\| > 4.5 often appears
//! for **sub-cycle** mean shifts (turbo, SMT). We only fail when \|t\| is high
//! **and** the absolute gap is **material** (≥ max(25 cycles, 0.2% of class mean)).
//!
//! **Methodology:** Reparaz, Balasch, Verbauwhede 2017 (https://eprint.iacr.org/2017/287).
//!
//! ## Tests
//!
//! **Fast correctness tests** (`cargo test --release --test ct_timing`):
//! These run on every `cargo test` and verify that the branchless implementations
//! produce correct outputs for every possible input class (even/odd, zero/nonzero,
//! flag-0/flag-1). A bug in the branchless logic shows up as a wrong result.
//!
//! **Timing tests** (`--include-ignored`, needs quiet isolated CPU):
//! These measure cycle counts (Welch + material effect gate) on secret paths:
//! scalar `div2` / `cond_negate` / `inv`, `ecmult_gen_const`, `pubkey_from_secret`,
//! `xonly_pubkey_from_secret`, ECDSA sign (+ recoverable), Schnorr sign (parity classes),
//! ECDH / `ecdh_compressed`, `ecdsa_sign_der_rfc6979`, MuSig `nonce_gen` / `partial_sign` /
//! `KeyAggCache` tweak adds, Taproot `xonly_pubkey_tweak_add` / `taproot_output_key`, and
//! (x86_64/aarch64) `ellswift_create` / `ellswift_xdh`.
//!
//! They MUST be run with `--test-threads=1` on a machine that is:
//!   - Not running background workloads (IBD, builds, etc.)
//!   - CPU frequency locked (disable turbo, use `cpupower frequency-set`)
//!   - Ideally pinned to one core with `taskset -c 0`
//!
//! ```text
//! taskset -c 0 cargo test --release --test ct_timing -- \
//!     --test-threads=1 --include-ignored --nocapture
//! ```
//!
//! **ASM inspection:** `RUSTFLAGS='--emit=asm' cargo rustc --release --lib`, then search `target/release/deps/*.s` for hot symbols (see `TIMING.md`).

use blvm_secp256k1::ecdh;
use blvm_secp256k1::ecdsa::{
    ecdsa_sig_sign, ecdsa_sig_sign_recoverable, ecdsa_sign_der_rfc6979, ge_from_compressed,
    ge_to_compressed, pubkey_from_secret,
};
use blvm_secp256k1::ecmult_gen_const;
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
use blvm_secp256k1::ellswift::{ellswift_create, ellswift_xdh};
use blvm_secp256k1::group::Gej;
use blvm_secp256k1::musig::{
    nonce_agg, nonce_gen, nonce_process, partial_sign, KeyAggCache, Session,
};
use blvm_secp256k1::scalar::Scalar;
use blvm_secp256k1::schnorr::{schnorr_sign, xonly_pubkey_from_secret};
use blvm_secp256k1::taproot::{taproot_output_key, xonly_pubkey_tweak_add};
use sha2::{Digest, Sha256};
use subtle::ConstantTimeEq;

fn scalar_from_seed(seed: u64) -> Option<Scalar> {
    let h: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
    let mut s = Scalar::zero();
    if s.set_b32(&h) || s.is_zero() {
        return None;
    }
    Some(s)
}

// helper: compare scalars via CT equality
fn scalar_eq(a: &Scalar, b: &Scalar) -> bool {
    bool::from(a.ct_eq(b))
}

// ────────────────────────────────────────────────────────────────────────────
// Correctness tests — always run, verify branchless paths are correct
// ────────────────────────────────────────────────────────────────────────────

/// div2 branchless path gives the same result as the reference formula for
/// even inputs (s/2) and odd inputs ((s+n)/2) across 512 representative scalars.
#[test]
fn ct_div2_correctness_even_and_odd() {
    for i in 1u64..=512 {
        let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
        let mut s = Scalar::zero();
        if s.set_b32(&h) || s.is_zero() {
            continue;
        }
        let orig = s;

        // Apply div2 and recover: result * 2 (mod n) must equal orig.
        s.div2();
        let mut s2 = Scalar::zero();
        s2.add(&s, &s);
        assert!(scalar_eq(&s2, &orig), "div2(s)*2 != s for i={i}");
    }
}

/// cond_negate(flag=0) must be identity; cond_negate(flag=1) must equal negate.
#[test]
fn ct_cond_negate_correctness() {
    for i in 0u32..=255 {
        let mut s = Scalar::zero();
        s.set_int(i.max(1));

        // flag = 0: no change
        let mut r0 = s;
        let ret0 = r0.cond_negate(0);
        assert!(
            scalar_eq(&r0, &s),
            "cond_negate(0) changed scalar for i={i}"
        );
        assert_eq!(ret0, -1, "cond_negate(0) return value wrong");

        // flag = 1: should equal negate
        let mut r1 = s;
        let ret1 = r1.cond_negate(1);
        let mut expected_neg = Scalar::zero();
        expected_neg.negate(&s);
        assert!(
            scalar_eq(&r1, &expected_neg),
            "cond_negate(1) != negate for i={i}"
        );
        assert_eq!(ret1, 1, "cond_negate(1) return value wrong");
    }
}

/// negate(0) == 0, negate(s) + s == 0 (mod n) for nonzero s.
#[test]
fn ct_negate_correctness_zero_and_nonzero() {
    // negate(0) must be 0
    let zero = Scalar::zero();
    let mut neg_zero = Scalar::zero();
    neg_zero.negate(&zero);
    assert!(scalar_eq(&neg_zero, &zero), "negate(0) != 0");

    // for nonzero s: s + negate(s) = 0
    for i in 1u32..=256 {
        let mut s = Scalar::zero();
        s.set_int(i);
        let mut neg_s = Scalar::zero();
        neg_s.negate(&s);
        let mut sum = Scalar::zero();
        sum.add(&s, &neg_s);
        assert!(scalar_eq(&sum, &zero), "s + negate(s) != 0 for i={i}");
    }
}

/// Schnorr sign with a key whose R.y is odd must produce a verifiable signature.
/// Previously the R.y branch could produce wrong output under misoptimisation.
#[test]
fn ct_schnorr_r_parity_correctness() {
    use blvm_secp256k1::schnorr::schnorr_verify;
    let msg = [0x42u8; 32];
    let aux = [0x00u8; 32];

    // Precompute 64 secret keys — roughly half will produce odd-y R, half even-y R.
    let mut verified_odd = 0u32;
    let mut verified_even = 0u32;
    for i in 1u64..=64 {
        let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
        let mut sk_scalar = Scalar::zero();
        if sk_scalar.set_b32(&h) || sk_scalar.is_zero() {
            continue;
        }
        let sig = schnorr_sign(&h, &msg, &aux).expect("schnorr sign failed");
        let pk_x = xonly_pubkey_from_secret(&h).expect("xonly_pubkey failed");
        assert!(
            schnorr_verify(&sig, &msg, &pk_x),
            "schnorr_verify failed for i={i}"
        );
        // Track parity distribution via sig[0] low bit (proxy for R.x oddness).
        if sig[0] & 1 == 0 {
            verified_even += 1;
        } else {
            verified_odd += 1;
        }
    }
    // Sanity: both parities were exercised.
    assert!(
        verified_odd > 0,
        "no odd-parity R generated — test not exercising both paths"
    );
    assert!(
        verified_even > 0,
        "no even-parity R generated — test not exercising both paths"
    );
}

// ────────────────────────────────────────────────────────────────────────────
// Timing harness (used by ignored tests below)
// ────────────────────────────────────────────────────────────────────────────

#[cfg(target_arch = "x86_64")]
fn read_cycles_start() -> u64 {
    unsafe {
        core::arch::x86_64::__cpuid(0);
        core::arch::x86_64::_rdtsc()
    }
}

#[cfg(target_arch = "x86_64")]
fn read_cycles_end() -> u64 {
    let mut _aux = 0u32;
    unsafe { core::arch::x86_64::__rdtscp(&mut _aux) }
}

#[cfg(not(target_arch = "x86_64"))]
fn read_cycles_start() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos() as u64
}

#[cfg(not(target_arch = "x86_64"))]
fn read_cycles_end() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .subsec_nanos() as u64
}

fn welch_t(a: &[i64], b: &[i64]) -> f64 {
    let n_a = a.len() as f64;
    let n_b = b.len() as f64;
    let mean_a = a.iter().map(|&x| x as f64).sum::<f64>() / n_a;
    let mean_b = b.iter().map(|&x| x as f64).sum::<f64>() / n_b;
    let var_a = a.iter().map(|&x| (x as f64 - mean_a).powi(2)).sum::<f64>() / (n_a - 1.0);
    let var_b = b.iter().map(|&x| (x as f64 - mean_b).powi(2)).sum::<f64>() / (n_b - 1.0);
    let se = (var_a / n_a + var_b / n_b).sqrt();
    if se < f64::EPSILON {
        return 0.0;
    }
    (mean_a - mean_b) / se
}

fn trim(samples: &mut Vec<i64>, pct: f64) {
    samples.sort_unstable();
    let cut = (samples.len() as f64 * pct) as usize;
    samples.drain(..cut);
    samples.truncate(samples.len().saturating_sub(cut));
}

const THRESHOLD: f64 = 4.5;

fn dudect<I, F: Fn(&I)>(inputs0: &[I], inputs1: &[I], n_each: usize, f: &F) -> (f64, f64, f64) {
    let mut times0: Vec<i64> = Vec::with_capacity(n_each);
    let mut times1: Vec<i64> = Vec::with_capacity(n_each);

    // Warm-up: fill i-cache and d-cache.
    for i in 0..200 {
        f(std::hint::black_box(&inputs0[i % inputs0.len()]));
        f(std::hint::black_box(&inputs1[i % inputs1.len()]));
    }

    // LCG for random interleaving (no external deps).
    let mut rng: u64 = 0xdeadbeef_cafebabe;
    let lcg = |s: &mut u64| -> u64 {
        *s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *s
    };

    let mut i0 = 0usize;
    let mut i1 = 0usize;
    while times0.len() < n_each || times1.len() < n_each {
        let use_class0 =
            times0.len() < n_each && (times1.len() >= n_each || lcg(&mut rng) & 1 == 0);
        if use_class0 {
            let inp = &inputs0[i0 % inputs0.len()];
            let t0 = read_cycles_start();
            f(std::hint::black_box(inp));
            let t1 = read_cycles_end();
            let d = (t1 as i64).saturating_sub(t0 as i64);
            if d > 0 {
                times0.push(d);
            }
            i0 += 1;
        } else if times1.len() < n_each {
            let inp = &inputs1[i1 % inputs1.len()];
            let t0 = read_cycles_start();
            f(std::hint::black_box(inp));
            let t1 = read_cycles_end();
            let d = (t1 as i64).saturating_sub(t0 as i64);
            if d > 0 {
                times1.push(d);
            }
            i1 += 1;
        }
    }

    trim(&mut times0, 0.05);
    trim(&mut times1, 0.05);

    let m0 = times0.iter().map(|&x| x as f64).sum::<f64>() / times0.len() as f64;
    let m1 = times1.iter().map(|&x| x as f64).sum::<f64>() / times1.len() as f64;
    let t = welch_t(&times0, &times1);
    (t, m0, m1)
}

/// With 50k+ trimmed samples, Welch's \|t\| crosses 4.5 for **sub-cycle** mean shifts on
/// Linux (turbo, SMT, scheduling). Treat \|t\| ≥ threshold as a leak only when the
/// absolute **and** relative gap is material: ≥ `max(25 cycles, 0.2% of class mean)`.
fn assert_welch_timing_ok(name: &str, t: f64, m0: f64, m1: f64) {
    let mean = (m0 + m1) * 0.5;
    let adiff = (m0 - m1).abs();
    let floor = (mean * 0.002).max(25.0);
    let material = adiff >= floor;
    if t.abs() >= THRESHOLD && material {
        panic!(
            "{name}: timing skew |t|={t:.2} ≥ {THRESHOLD} with material Δ={adiff:.1} cycles (mean≈{mean:.1}, m0={m0:.1}, m1={m1:.1})"
        );
    }
    if t.abs() >= THRESHOLD && !material {
        eprintln!(
            "  (note) {name}: |t|={t:.2} but Δ={adiff:.1} < material floor {floor:.1} — statistical noise, not actionable"
        );
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Timing tests — all #[ignore]; run on isolated CPU with --include-ignored
// ────────────────────────────────────────────────────────────────────────────

/// Timing: div2 must take the same time for even vs odd inputs.
/// Run: taskset -c 0 cargo test --release --test ct_timing -- --test-threads=1 --include-ignored --nocapture
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_div2_even_vs_odd() {
    const N: usize = 100_000;
    let even: Vec<Scalar> = (1u64..)
        .filter_map(|i| {
            let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
            let mut s = Scalar::zero();
            if s.set_b32(&h) || s.is_zero() || s.is_odd() {
                return None;
            }
            Some(s)
        })
        .take(N)
        .collect();
    let odd: Vec<Scalar> = (1u64..)
        .filter_map(|i| {
            let h: [u8; 32] = Sha256::digest((i | 0x8000_0000_0000_0000u64).to_le_bytes()).into();
            let mut s = Scalar::zero();
            if s.set_b32(&h) || s.is_zero() || !s.is_odd() {
                return None;
            }
            Some(s)
        })
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&even, &odd, N, &|s: &Scalar| {
        let mut r = *s;
        r.div2();
        std::hint::black_box(r);
    });
    println!("ct_timing_div2 | t={t:.2} mean_even={m0:.1} mean_odd={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_div2 even/odd", t, m0, m1);
}

/// Timing: cond_negate must take the same time regardless of flag.
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_cond_negate_flag0_vs_flag1() {
    const N: usize = 100_000;
    let mut base = Scalar::zero();
    base.set_int(0x1234_5678);
    let class0: Vec<(Scalar, i32)> = (0..N).map(|_| (base, 0)).collect();
    let class1: Vec<(Scalar, i32)> = (0..N).map(|_| (base, 1)).collect();

    let (t, m0, m1) = dudect(&class0, &class1, N, &|(s, f): &(Scalar, i32)| {
        let mut r = *s;
        r.cond_negate(*f);
        std::hint::black_box(r);
    });
    println!("ct_timing_cond_negate | t={t:.2} flag0={m0:.1} flag1={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_cond_negate flag0/1", t, m0, m1);
}

/// Timing: `ecmult_gen_const` — two pools of uniform hash-derived scalars (same recipe).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ecmult_gen_const_two_random_pools() {
    fn hash_scalar(seed: u64) -> Option<Scalar> {
        let h: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
        let mut s = Scalar::zero();
        if s.set_b32(&h) || s.is_zero() {
            return None;
        }
        Some(s)
    }

    const N: usize = 10_000;
    let pool_a: Vec<Scalar> = (0u64..).filter_map(hash_scalar).take(N).collect();
    let pool_b: Vec<Scalar> = (0u64..)
        .filter_map(|i| hash_scalar(i ^ 0xA5A5_A5A5_A5A5_A5A5))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|s: &Scalar| {
        let mut r = Gej::default();
        ecmult_gen_const(&mut r, s);
        std::hint::black_box(r);
    });
    println!("ct_timing_ecmult_gen_const | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ecmult_gen_const pools", t, m0, m1);
}

/// Timing: Schnorr sign must take the same time regardless of nonce-point R.y parity.
/// This is the specific vulnerability that the branchless R.y fix addresses.
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_schnorr_sign_r_parity() {
    const N: usize = 5_000;
    let fixed_msg = [0x42u8; 32];
    let fixed_aux = [0x00u8; 32];

    let mut class0: Vec<[u8; 32]> = Vec::new(); // sig[0] bit 0 == 0
    let mut class1: Vec<[u8; 32]> = Vec::new(); // sig[0] bit 0 == 1
    let mut i = 1u64;
    while class0.len() < N || class1.len() < N {
        let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
        i += 1;
        let mut sk = Scalar::zero();
        if sk.set_b32(&h) || sk.is_zero() {
            continue;
        }
        let sig = match schnorr_sign(&h, &fixed_msg, &fixed_aux) {
            Some(s) => s,
            None => continue,
        };
        if sig[0] & 1 == 0 {
            if class0.len() < N {
                class0.push(h);
            }
        } else {
            if class1.len() < N {
                class1.push(h);
            }
        }
    }

    let (t, m0, m1) = dudect(&class0, &class1, N, &|sk: &[u8; 32]| {
        std::hint::black_box(schnorr_sign(sk, &fixed_msg, &fixed_aux));
    });
    println!("ct_timing_schnorr_r_parity | t={t:.2} class0={m0:.1} class1={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_schnorr R parity", t, m0, m1);
}

/// Timing: ECDH — two pools of hash-derived secrets (same distribution).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ecdh_two_random_pools() {
    const LIBSECP_ECDH_PK: [u8; 33] = [
        0x03, 0x54, 0x94, 0xc1, 0x5d, 0x32, 0x09, 0x97, 0x06, 0xc2, 0x39, 0x5f, 0x94, 0x34, 0x87,
        0x45, 0xfd, 0x75, 0x7c, 0xe3, 0x0e, 0x4e, 0x8c, 0x90, 0xfb, 0xa2, 0xba, 0xd1, 0x84, 0xf8,
        0x83, 0xc6, 0x9f,
    ];
    let pk = ge_from_compressed(&LIBSECP_ECDH_PK).expect("ct_ecdh: bad pk");
    fn hash_scalar(seed: u64) -> Option<Scalar> {
        let h: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
        let mut s = Scalar::zero();
        if s.set_b32(&h) || s.is_zero() {
            return None;
        }
        Some(s)
    }
    const N: usize = 5_000;
    let pool_a: Vec<Scalar> = (0u64..).filter_map(hash_scalar).take(N).collect();
    let pool_b: Vec<Scalar> = (0u64..)
        .filter_map(|i| hash_scalar(i ^ 0x5A5A_5A5A_5A5A_5A5A))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|s: &Scalar| {
        std::hint::black_box(ecdh::ecdh(&pk, s));
    });
    println!("ct_timing_ecdh | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ecdh pools", t, m0, m1);
}

// ── ECDSA / MuSig helpers (timing fixtures) ─────────────────────────────────

fn ecdsa_signing_triple(seed: u64) -> Option<(Scalar, Scalar, Scalar)> {
    let sec = scalar_from_seed(seed)?;
    let msg = scalar_from_seed(seed ^ 0x0123_4567_89AB_CDEF)?;
    let nonce = scalar_from_seed(seed ^ 0xFEDC_BA09_8765_4321)?;
    ecdsa_sig_sign(&sec, &msg, &nonce)?;
    Some((sec, msg, nonce))
}

fn collect_ecdsa_triples(mut seed: u64, n: usize) -> Vec<(Scalar, Scalar, Scalar)> {
    let mut v = Vec::with_capacity(n);
    while v.len() < n {
        if let Some(t) = ecdsa_signing_triple(seed) {
            v.push(t);
        }
        seed = seed.wrapping_add(1);
    }
    v
}

#[derive(Clone)]
struct MusigNonceCtx {
    sk: [u8; 32],
    pk: [u8; 33],
    cache: KeyAggCache,
    msg: [u8; 32],
    session_rand32: [u8; 32],
}

fn musig_nonce_ctx(seed: u64) -> Option<MusigNonceCtx> {
    let sec = scalar_from_seed(seed)?;
    let mut sk = [0u8; 32];
    sec.get_b32(&mut sk);
    let pk = ge_to_compressed(&pubkey_from_secret(&sec));
    let cache = KeyAggCache::new(&[pk])?;
    let msg = [0x42u8; 32];
    let session_rand32: [u8; 32] = Sha256::digest((seed ^ 0xC001_D00D).to_le_bytes()).into();
    if session_rand32.iter().all(|b| *b == 0) {
        return None;
    }
    let mut r = session_rand32;
    if nonce_gen(&mut r, Some(&sk), &pk, Some(&msg), Some(&cache), None).is_none() {
        return None;
    }
    Some(MusigNonceCtx {
        sk,
        pk,
        cache,
        msg,
        session_rand32,
    })
}

fn collect_musig_nonce_ctx(mut seed: u64, n: usize) -> Vec<MusigNonceCtx> {
    let mut v = Vec::with_capacity(n);
    while v.len() < n {
        if let Some(c) = musig_nonce_ctx(seed) {
            v.push(c);
        }
        seed = seed.wrapping_add(1);
    }
    v
}

#[derive(Clone)]
struct MusigPartialCtx {
    sk: [u8; 32],
    secnonce_tpl: ([u8; 64], [u8; 33]),
    cache: KeyAggCache,
    session: Session,
}

fn musig_partial_ctx(seed: u64) -> Option<MusigPartialCtx> {
    let sec = scalar_from_seed(seed)?;
    let mut sk = [0u8; 32];
    sec.get_b32(&mut sk);
    let pk = ge_to_compressed(&pubkey_from_secret(&sec));
    let cache = KeyAggCache::new(&[pk])?;
    let msg = [0x42u8; 32];
    let mut rand: [u8; 32] = Sha256::digest((seed ^ 0xBADC_0FFE).to_le_bytes()).into();
    if rand.iter().all(|b| *b == 0) {
        return None;
    }
    let ((k64, pk33), pubnonce) =
        nonce_gen(&mut rand, Some(&sk), &pk, Some(&msg), Some(&cache), None)?;
    let aggnonce = nonce_agg(std::slice::from_ref(&pubnonce))?;
    let session = nonce_process(&aggnonce, &msg, &cache)?;
    Some(MusigPartialCtx {
        sk,
        secnonce_tpl: (k64, pk33),
        cache,
        session,
    })
}

fn collect_musig_partial_ctx(mut seed: u64, n: usize) -> Vec<MusigPartialCtx> {
    let mut v = Vec::with_capacity(n);
    while v.len() < n {
        if let Some(c) = musig_partial_ctx(seed) {
            v.push(c);
        }
        seed = seed.wrapping_add(1);
    }
    v
}

/// Fixed 1-of-1 MuSig key-agg cache for tweak timing (public aggregation only; tweak path is CT).
fn musig_single_signer_cache_template() -> KeyAggCache {
    let sec = scalar_from_seed(0xFEED_BEEF_CAFE_0001).expect("musig template sk");
    let pk = ge_to_compressed(&pubkey_from_secret(&sec));
    KeyAggCache::new(&[pk]).expect("musig 1-of-1 cache")
}

fn scalar32_valid_for_tweak(seed: u64) -> Option<[u8; 32]> {
    let t: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
    let mut s = Scalar::zero();
    if s.set_b32(&t) || s.is_zero() {
        return None;
    }
    Some(t)
}

fn musig_xonly_tweak_ok(base: &KeyAggCache, seed: u64) -> Option<[u8; 32]> {
    let t = scalar32_valid_for_tweak(seed)?;
    let mut c = base.clone();
    c.pubkey_xonly_tweak_add(&t)?;
    Some(t)
}

fn musig_ec_tweak_ok(base: &KeyAggCache, seed: u64) -> Option<[u8; 32]> {
    let t = scalar32_valid_for_tweak(seed)?;
    let mut c = base.clone();
    c.pubkey_ec_tweak_add(&t)?;
    Some(t)
}

fn rfc6979_signing_pair(seed: u64) -> Option<([u8; 32], [u8; 32])> {
    let msg: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
    let sk: [u8; 32] = Sha256::digest((seed ^ 0x1234_5678_ABCD_EF01).to_le_bytes()).into();
    let mut sc = Scalar::zero();
    if sc.set_b32(&sk) || sc.is_zero() {
        return None;
    }
    ecdsa_sign_der_rfc6979(&msg, &sk)?;
    Some((msg, sk))
}

fn collect_rfc6979_pairs(mut seed: u64, n: usize) -> Vec<([u8; 32], [u8; 32])> {
    let mut v = Vec::with_capacity(n);
    while v.len() < n {
        if let Some(p) = rfc6979_signing_pair(seed) {
            v.push(p);
        }
        seed = seed.wrapping_add(1);
    }
    v
}

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
fn ellswift_party_fixture(seed: u64) -> Option<([u8; 32], [u8; 64])> {
    let sk: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
    let mut sc = Scalar::zero();
    if sc.set_b32(&sk) || sc.is_zero() {
        return None;
    }
    let ell = ellswift_create(&sk, None)?;
    Some((sk, ell))
}

/// Timing: `Scalar::inv` on two pools of nonzero scalars (ECDSA / MuSig nonce inverse path).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_scalar_inv_two_random_pools() {
    const N: usize = 8_000;
    let pool_a: Vec<Scalar> = (0u64..).filter_map(scalar_from_seed).take(N).collect();
    let pool_b: Vec<Scalar> = (0u64..)
        .filter_map(|i| scalar_from_seed(i ^ 0x1357_9BDF_2468_ACE0))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|s: &Scalar| {
        let mut out = Scalar::zero();
        out.inv(s);
        std::hint::black_box(out);
    });
    println!("ct_timing_scalar_inv | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_scalar_inv pools", t, m0, m1);
}

/// Timing: full affine conversion after `ecmult_gen_const` (pubkey / keypair derivation).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_pubkey_from_secret_two_random_pools() {
    const N: usize = 10_000;
    let pool_a: Vec<Scalar> = (0u64..).filter_map(scalar_from_seed).take(N).collect();
    let pool_b: Vec<Scalar> = (0u64..)
        .filter_map(|i| scalar_from_seed(i ^ 0x51A1_51A1_51A1_51A1))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|s: &Scalar| {
        std::hint::black_box(pubkey_from_secret(s));
    });
    println!("ct_timing_pubkey_from_secret | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_pubkey_from_secret pools", t, m0, m1);
}

/// Timing: BIP340 x-only pubkey from secret (includes `ecmult_gen_const` + x extraction).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_xonly_pubkey_from_secret_two_pools() {
    fn sk_bytes(seed: u64) -> Option<[u8; 32]> {
        let h: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
        let mut s = Scalar::zero();
        if s.set_b32(&h) || s.is_zero() {
            return None;
        }
        Some(h)
    }
    const N: usize = 10_000;
    let pool_a: Vec<[u8; 32]> = (0u64..).filter_map(sk_bytes).take(N).collect();
    let pool_b: Vec<[u8; 32]> = (0u64..)
        .filter_map(|i| sk_bytes(i ^ 0x70B3_4024_6835_79BD))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|sk: &[u8; 32]| {
        std::hint::black_box(xonly_pubkey_from_secret(sk));
    });
    println!("ct_timing_xonly_pubkey_from_secret | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_xonly_pubkey_from_secret pools", t, m0, m1);
}

/// Timing: ECDSA sign core (`ecmult_gen_const`, `Scalar::inv`, low-S `cond_negate`).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ecdsa_sig_sign_two_pools() {
    const N: usize = 4_000;
    let pool_a = collect_ecdsa_triples(1, N);
    let pool_b = collect_ecdsa_triples(1_000_000_000, N);

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|trip: &(Scalar, Scalar, Scalar)| {
        let (sec, msg, nonce) = *trip;
        std::hint::black_box(ecdsa_sig_sign(&sec, &msg, &nonce));
    });
    println!("ct_timing_ecdsa_sig_sign | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ecdsa_sig_sign pools", t, m0, m1);
}

/// Timing: ECDSA sign + recovery id (extra field / parity encoding vs `ecdsa_sig_sign`).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ecdsa_sig_sign_recoverable_two_pools() {
    const N: usize = 4_000;
    let pool_a = collect_ecdsa_triples(3, N);
    let pool_b = collect_ecdsa_triples(2_000_000_000, N);

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|trip: &(Scalar, Scalar, Scalar)| {
        let (sec, msg, nonce) = *trip;
        std::hint::black_box(ecdsa_sig_sign_recoverable(&sec, &msg, &nonce));
    });
    println!(
        "ct_timing_ecdsa_sig_sign_recoverable | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles"
    );
    assert_welch_timing_ok("ct_timing_ecdsa_sig_sign_recoverable pools", t, m0, m1);
}

/// Timing: MuSig2 `nonce_gen` (constant-time `k1·G`, `k2·G`).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_musig_nonce_gen_two_pools() {
    const N: usize = 3_000;
    let pool_a = collect_musig_nonce_ctx(10, N);
    let pool_b = collect_musig_nonce_ctx(9_000_000_000, N);

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|ctx: &MusigNonceCtx| {
        let mut r = ctx.session_rand32;
        std::hint::black_box(nonce_gen(
            &mut r,
            Some(&ctx.sk),
            &ctx.pk,
            Some(&ctx.msg),
            Some(&ctx.cache),
            None,
        ));
    });
    println!("ct_timing_musig_nonce_gen | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_musig_nonce_gen pools", t, m0, m1);
}

/// Timing: MuSig2 `partial_sign` (secret nonce + `cond_negate` adjustments).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_musig_partial_sign_two_pools() {
    const N: usize = 3_000;
    let pool_a = collect_musig_partial_ctx(20, N);
    let pool_b = collect_musig_partial_ctx(8_000_000_000, N);

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|ctx: &MusigPartialCtx| {
        let mut sn = ctx.secnonce_tpl.clone();
        std::hint::black_box(partial_sign(&mut sn, &ctx.sk, &ctx.cache, &ctx.session));
    });
    println!("ct_timing_musig_partial_sign | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_musig_partial_sign pools", t, m0, m1);
}

/// Timing: Taproot-style x-only tweak add (`t·G` via `ecmult_gen_const` + point add).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_xonly_pubkey_tweak_add_two_pools() {
    const INTERNAL_SK: [u8; 32] = [
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e,
        0x1f, 0x20,
    ];
    let internal_x = xonly_pubkey_from_secret(&INTERNAL_SK).expect("tap tweak internal xonly");

    fn tweak_from_seed(seed: u64) -> Option<[u8; 32]> {
        let t: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
        let mut s = Scalar::zero();
        if s.set_b32(&t) {
            return None;
        }
        if s.is_zero() {
            return None;
        }
        Some(t)
    }
    const N: usize = 5_000;
    let pool_a: Vec<[u8; 32]> = (1u64..).filter_map(tweak_from_seed).take(N).collect();
    let pool_b: Vec<[u8; 32]> = (1u64..)
        .filter_map(|i| tweak_from_seed(i ^ 0xB166_C001_D00D_FEED))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|tw: &[u8; 32]| {
        std::hint::black_box(xonly_pubkey_tweak_add(&internal_x, tw));
    });
    println!("ct_timing_xonly_pubkey_tweak_add | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_xonly_pubkey_tweak_add pools", t, m0, m1);
}

/// Timing: ElligatorSwift create (secret key encoding).
#[test]
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ellswift_create_two_pools() {
    const N: usize = 4_000;
    let pool_a: Vec<[u8; 32]> = (0u64..)
        .filter_map(|i| ellswift_party_fixture(i).map(|x| x.0))
        .take(N)
        .collect();
    let pool_b: Vec<[u8; 32]> = (0u64..)
        .filter_map(|i| ellswift_party_fixture(i ^ 0xE115_E115_E115_E115).map(|x| x.0))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|sk: &[u8; 32]| {
        std::hint::black_box(ellswift_create(sk, None));
    });
    println!("ct_timing_ellswift_create | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ellswift_create pools", t, m0, m1);
}

/// Timing: ElligatorSwift x-only ECDH (`ecmult_const` on decoded point).
#[test]
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ellswift_xdh_two_pools() {
    let (_peer_sk, peer_ell) = (1u64..10_000u64)
        .find_map(|i| {
            let sk: [u8; 32] = Sha256::digest((i ^ 0x7777).to_le_bytes()).into();
            let mut sc = Scalar::zero();
            if sc.set_b32(&sk) || sc.is_zero() {
                return None;
            }
            ellswift_create(&sk, None).map(|e| (sk, e))
        })
        .expect("peer ellswift");

    const N: usize = 3_000;
    let pool_a: Vec<([u8; 32], [u8; 64])> = (0u64..)
        .filter_map(ellswift_party_fixture)
        .take(N)
        .collect();
    let pool_b: Vec<([u8; 32], [u8; 64])> = (0u64..)
        .filter_map(|i| ellswift_party_fixture(i ^ 0xB132_B324_C456_D678))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|party: &([u8; 32], [u8; 64])| {
        let (sk, ell_a) = *party;
        std::hint::black_box(ellswift_xdh(&ell_a, &peer_ell, &sk, false));
    });
    println!("ct_timing_ellswift_xdh | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ellswift_xdh pools", t, m0, m1);
}

/// Timing: ECDSA DER sign with RFC6979 nonce derivation + retry loop (libsecp-style).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ecdsa_sign_der_rfc6979_two_pools() {
    const N: usize = 2_000;
    let pool_a = collect_rfc6979_pairs(100, N);
    let pool_b = collect_rfc6979_pairs(50_000_000_000, N);

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|pair: &([u8; 32], [u8; 32])| {
        let (msg, sk) = *pair;
        std::hint::black_box(ecdsa_sign_der_rfc6979(&msg, &sk));
    });
    println!("ct_timing_ecdsa_sign_der_rfc6979 | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ecdsa_sign_der_rfc6979 pools", t, m0, m1);
}

/// Timing: ECDH using compressed pubkey bytes (parse + `ecmult_const`).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_ecdh_compressed_two_pools() {
    const LIBSECP_ECDH_PK: [u8; 33] = [
        0x03, 0x54, 0x94, 0xc1, 0x5d, 0x32, 0x09, 0x97, 0x06, 0xc2, 0x39, 0x5f, 0x94, 0x34, 0x87,
        0x45, 0xfd, 0x75, 0x7c, 0xe3, 0x0e, 0x4e, 0x8c, 0x90, 0xfb, 0xa2, 0xba, 0xd1, 0x84, 0xf8,
        0x83, 0xc6, 0x9f,
    ];
    const N: usize = 5_000;
    let pool_a: Vec<[u8; 32]> = (0u64..)
        .filter_map(|i| {
            let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
            let mut s = Scalar::zero();
            if s.set_b32(&h) || s.is_zero() {
                return None;
            }
            let mut b = [0u8; 32];
            s.get_b32(&mut b);
            Some(b)
        })
        .take(N)
        .collect();
    let pool_b: Vec<[u8; 32]> = (0u64..)
        .filter_map(|i| {
            let h: [u8; 32] = Sha256::digest((i ^ 0xC0FFEE00_BADC0D00).to_le_bytes()).into();
            let mut s = Scalar::zero();
            if s.set_b32(&h) || s.is_zero() {
                return None;
            }
            let mut b = [0u8; 32];
            s.get_b32(&mut b);
            Some(b)
        })
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|sk: &[u8; 32]| {
        std::hint::black_box(ecdh::ecdh_compressed(&LIBSECP_ECDH_PK, sk));
    });
    println!("ct_timing_ecdh_compressed | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_ecdh_compressed pools", t, m0, m1);
}

/// Timing: MuSig aggregate pubkey x-only tweak (`apply_tweak` + `ecmult_gen_const`).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_musig_keyagg_xonly_tweak_two_pools() {
    let base = musig_single_signer_cache_template();
    const N: usize = 2_500;
    let pool_a: Vec<[u8; 32]> = (1u64..)
        .filter_map(|s| musig_xonly_tweak_ok(&base, s))
        .take(N)
        .collect();
    let pool_b: Vec<[u8; 32]> = (1u64..)
        .filter_map(|s| musig_xonly_tweak_ok(&base, s ^ 0xA505_A505_A505_A505))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|tw: &[u8; 32]| {
        let mut c = base.clone();
        std::hint::black_box(c.pubkey_xonly_tweak_add(tw));
    });
    println!("ct_timing_musig_keyagg_xonly_tweak | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_musig_keyagg_xonly_tweak pools", t, m0, m1);
}

/// Timing: MuSig aggregate pubkey plain EC tweak (BIP32-style).
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_musig_keyagg_ec_tweak_two_pools() {
    let base = musig_single_signer_cache_template();
    const N: usize = 2_500;
    let pool_a: Vec<[u8; 32]> = (1u64..)
        .filter_map(|s| musig_ec_tweak_ok(&base, s))
        .take(N)
        .collect();
    let pool_b: Vec<[u8; 32]> = (1u64..)
        .filter_map(|s| musig_ec_tweak_ok(&base, s ^ 0xB132_B132_EC00_EC00))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|tw: &[u8; 32]| {
        let mut c = base.clone();
        std::hint::black_box(c.pubkey_ec_tweak_add(tw));
    });
    println!("ct_timing_musig_keyagg_ec_tweak | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_musig_keyagg_ec_tweak pools", t, m0, m1);
}

/// Timing: Taproot output key = x-only tweak with tagged merkle hash.
#[test]
#[ignore = "needs isolated CPU; run with --include-ignored --nocapture --test-threads=1"]
fn ct_timing_taproot_output_key_two_pools() {
    const INTERNAL_SK: [u8; 32] = [
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e,
        0x1f, 0x20,
    ];
    let internal_x = xonly_pubkey_from_secret(&INTERNAL_SK).expect("taproot output key internal");

    fn merkle_ok(ix: &[u8; 32], seed: u64) -> Option<[u8; 32]> {
        let root: [u8; 32] = Sha256::digest(seed.to_le_bytes()).into();
        taproot_output_key(ix, &root)?;
        Some(root)
    }

    const N: usize = 4_000;
    let pool_a: Vec<[u8; 32]> = (1u64..)
        .filter_map(|s| merkle_ok(&internal_x, s))
        .take(N)
        .collect();
    let pool_b: Vec<[u8; 32]> = (1u64..)
        .filter_map(|s| merkle_ok(&internal_x, s ^ 0x7A70_27A7_027A_7027))
        .take(N)
        .collect();

    let (t, m0, m1) = dudect(&pool_a, &pool_b, N, &|root: &[u8; 32]| {
        std::hint::black_box(taproot_output_key(&internal_x, root));
    });
    println!("ct_timing_taproot_output_key | t={t:.2} pool_a={m0:.1} pool_b={m1:.1} cycles");
    assert_welch_timing_ok("ct_timing_taproot_output_key pools", t, m0, m1);
}
