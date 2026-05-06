# Timing and constant-time (side-channel) contract

This document describes which APIs in `blvm-secp256k1` are intended to run in *constant time* (no secret-dependent branches, indices, or early returns on secret data) in the same sense as **libsecp256k1** for comparable operations, and which APIs are explicitly **variable time** and safe only for public inputs.

**Security boundary (Bitcoin-relevant):** Any routine that takes a **private scalar** (secret key, ECDSA or Schnorr nonce, ECDH private scalar, MuSig secret nonce / partial signing with secret) must use the **constant-time** path. Routines that only consume **public** keys, **public** points, and **message hashes** may use faster variable-time code.

## Constant-time primitives — implementation status

| Primitive | Location | How CT is achieved |
|-----------|----------|--------------------|
| `ecmult_gen_const` (k·G, signer path) | `ecmult_gen_comb.rs` | Multi-comb table; cmov scan over all 32 entries per block; `add_ge` / `double_ct` (Brier-Joye); no secret-dependent branch |
| `ecmult_const` (q·A, ECDH / verify-of-secret) | `ecmult_const.rs` | GLV + cmov table scan (`ecmult_const_table_get_ge`); unified `add_ge` / `double_ct` |
| `Scalar::div2` / `half_modn` | `scalar.rs` | Branchless: builds `add_mask = 0u64.wrapping_sub(odd_bit)`, adds `N & mask` before shifting — no branch on parity |
| `Scalar::cond_negate` | `scalar.rs` | Mask built via `0u64.wrapping_sub((flag != 0) as u64)`; return value via mask MSB arithmetic — no branch on secret flag |
| Schnorr `schnorr_sign_inner` — R.y parity | `schnorr.rs` | `k.cond_negate(parity)` + `FieldElement::cmov` with `subtle::Choice` — no branch on nonce-derived parity |
| `Scalar::inv` | `scalar.rs` | safegcd (`modinv64`) on x86_64 / aarch64; Fermat fallback on other arches (not timing-safe; see source) |
| `FieldElement::cmov` | `field/layout_5x52.rs` | `subtle::ConditionallySelectable` |

## API matrix

| API / area | CT expectation | Notes |
|------------|---------------|--------|
| `ecdsa::ecdsa_sig_verify` | Variable time OK | Verify workload. |
| `ecdsa::ecdsa_sig_sign`, `ecdsa::ecdsa_sig_sign_recoverable` | **Constant time** (secret) | `ecmult_gen_const` (comb) for nonce k·G; `Scalar::inv` (safegcd) for s⁻¹; `cond_negate` for low-S. |
| `ecdsa::pubkey_from_secret` | **Constant time** (secret) | `ecmult_gen_const` (comb). |
| `schnorr::schnorr_sign`, `schnorr::xonly_pubkey_from_secret` | **Constant time** (secret) | `ecmult_gen_const` for k·G and d·G; branchless R.y parity. |
| `schnorr::schnorr_verify`, `schnorr::schnorr_verify_batch` | Variable time OK | Public verification. |
| `ecdh::ecdh`, `ecdh::ecdh_compressed` | **Constant time** (secret) | `ecmult_const`; branchless `half_modn` in GLV decomp. |
| `ellswift::ellswift_xdh`, `ellswift_create` | **Constant time** (secret) | `ecmult_const` for the EC step. |
| `musig` secret paths | **Constant time** (secret) | `ecmult_gen_const` and `Scalar::inv` where applicable. |
| `ecmult::ecmult` | Variable time | Public / verify. |
| `ecmult::ecmult_gen` | Variable time | Fast path; not for hardened secret k·G. |
| `ecmult::ecmult_const` | **Constant time** (secret q) | Port of libsecp256k1 `ecmult_const` style. |
| `Scalar::inv` | **CT** (x86_64 / aarch64) | safegcd; other archs: see source. |
| `Scalar::inv_var` | Alias of `inv` | Backwards compatible. |

## Callsite summary

| Secret-path API | CT primitive |
|-----------------|-------------|
| `pubkey_from_secret` | `ecmult_gen_const` |
| `ecdsa_sig_sign` / `ecdsa_sig_sign_recoverable` | `ecmult_gen_const`, `Scalar::inv` (nonce), `Scalar::cond_negate` (low-S) |
| `schnorr_sign`, `xonly_pubkey_from_secret` | `ecmult_gen_const`, branchless R.y via `cond_negate` + `FieldElement::cmov` |
| `ecdh` / `ecdh_compressed` | `ecmult_const` (branchless `half_modn` in decomposition) |
| `ellswift::ellswift_xdh` | `ecmult_const` |
| `musig::nonce_gen`, `KeyAggCache::apply_tweak` | `ecmult_gen_const` |
| `taproot::xonly_pubkey_tweak_add` | `ecmult_gen_const` |

Verify and batch paths keep `ecmult` / `ecmult_gen` / `inv` on public-data workloads.

## Verifying CT against the compiler and CPU

Rust/LLVM **do not** promise that “branchless-looking” source stays branchless or constant-latency after optimization. You verify on a **concrete** triple (crate revision, `RUSTFLAGS`, target CPU).

**Suggested local check:** print CPU identity (`lscpu` / `uname -a`), run algebraic tests, then statistical timing pinned to one core:

```text
cargo test --release --test ct_timing -- --test-threads=1 --nocapture
taskset -c 0 cargo test --release --test ct_timing -- \
  --test-threads=1 --include-ignored --nocapture
```

Optionally inspect asm yourself: `RUSTFLAGS='--emit=asm' cargo rustc --release --lib` and review `target/release/deps/*.s` for hot symbols.

**`tests/ct_timing.rs` statistical targets** (with `--include-ignored`): Welch + material effect gate on two input classes of the same recipe — `Scalar::div2`, `cond_negate`, `inv`, `ecmult_gen_const`, `pubkey_from_secret`, `xonly_pubkey_from_secret`, `ecdsa_sig_sign`, `ecdsa_sig_sign_recoverable`, `ecdsa_sign_der_rfc6979`, Schnorr sign (R.y parity split), `ecdh`, `ecdh_compressed`, MuSig `nonce_gen` / `partial_sign` / `KeyAggCache::pubkey_xonly_tweak_add` / `pubkey_ec_tweak_add`, Taproot `xonly_pubkey_tweak_add` / `taproot_output_key`, and (x86_64/aarch64 only) ElligatorSwift `ellswift_create` / `ellswift_xdh`. **Still not** a timing target: `KeyAggCache::new` (public key aggregation only), verify/batch paths, and full BIP workflows that are mostly public (e.g. `nonce_agg` over published nonces).

| Approach | What it catches |
|----------|-----------------|
| **Inspect asm** | `cargo rustc --release -- --emit asm` (or `llvm-objdump -d`) on the hot symbol; look for **conditional jumps** (`jne`, `jb`, …) whose condition could depend on secret-loaded registers. Compare two builds (e.g. with/without a change). |
| **Dudect** | Statistical timing over many inputs; good at spotting **micro-architectural** leakage (not only branches). [google/dudect](https://github.com/google/dudect) is the usual starting point; wrap a C ABI that calls into your Rust `staticlib`/`cdylib` or a thin `#[no_mangle]` harness. |
| **ctgrind / FlowTracker** | Dynamic or static analysis that flags secret-dependent control flow; heavier setup, stronger when it applies. |
| **Same binary, controlled environment** | Disable turbo if you care about stable cycles; pin CPU; huge sample counts for microbenchmarks. |

**Practical hardening** (libsecp-style patterns already used in places here): `core::hint::black_box` on equality flags before table scans; `#[inline(never)]` on barrier functions if you need to stop a specific hoist; **volatile** loads of mask bytes in the hottest cmov loops if a future LLVM revision starts “optimizing away” your intent (rare but documented in CT literature).

**What verification does *not* replace:** cache timing, power analysis, speculative execution — those need different labs and threat models.

## Known limitations

- **`subtle::Choice`-based cmov in `ecmult_const`** is compiler-mediated CT; the `subtle` crate documents its best-effort model. For formal guarantees, use dudect/ctgrind to verify on the target platform.
- **`Scalar::inv` fallback** on non-x86_64/aarch64 is Fermat (VT). Don't use on those targets for secret inputs.
- **`pk_parity_odd` branch** in `schnorr_sign_inner` (`if pk_parity_odd { d_adj.negate(d) }`) is key-fixed (same parity for all signatures under a given key). The effective signing scalar is public given the public key, so this is not a per-nonce timing leak.

## Performance / regression

This repo does not ship Criterion benches or bench scripts. When you change hot paths, compare **before/after** on the same machine (same CPU governor and load) using your own harness — e.g. wall-clock or `perf` around representative `cargo test` subsets, or a private Criterion/binary driver. Treat a sustained **>3%** shift as worth investigating; expect ~0–3% noise from the environment alone.
