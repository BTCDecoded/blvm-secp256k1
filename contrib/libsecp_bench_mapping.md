# bitcoin-core **secp256k1** benchmarks — row semantics

Upstream: [bitcoin-core/secp256k1](https://github.com/bitcoin-core/secp256k1), e.g. clone to `~/src/secp256k1`.

This crate **does not** ship Criterion benches. Use libsecp’s numbers as a **reference**; compare to blvm with your own driver (tests, `perf`, or a private harness) if you need A/B timings.

```bash
cd ~/src/secp256k1
./autogen.sh && ./configure --enable-benchmark
make -j"$(nproc)" bench bench_ecmult
```

Or **`contrib/run_libsecp_parity.sh`** (`LIBSECP_SRC`, `SECP256K1_BENCH_ITERS`).

## Units

| Tool | Output |
|------|--------|
| **libsecp** | **Min / Avg / Max µs** per inner-loop iteration (10 outer rounds). |

## API-level (`./bench …`) — what each row exercises

| libsecp `bench` | Conceptual blvm analogue |
|-----------------|--------------------------|
| **`schnorrsig_sign`** | BIP340 `schnorr_sign`; libsecp uses rotating keypair slots (`bench_impl.h` byte pattern). |
| **`schnorrsig_verify`** | `schnorr_verify` on precomputed sigs/pubkeys. |
| **`ec_keygen`** | Chain of `pubkey_from_secret` / scalar-from-x like `bench.c` `bench_keygen_run`. |
| **`ecdh`** | `ecdh` / `ecdh_compressed` with libsecp’s fixed bench point + scalar (`modules/ecdh/bench_impl.h`). |
| **`ecdsa_sign`** | RFC6979 + DER loop (`ecdsa_sign_der_rfc6979` and inner `ecdsa_sig_sign`). |

**EllSwift:** with EllSwift enabled, `./bench` may include **`ellswift_*`** rows; use the plain **`ecdh`** row for vanilla ECDH.

## Primitives (`./bench_ecmult`)

| libsecp row | Notes |
|-------------|--------|
| **`ecmult_gen`** | CT-style table `k*G` (signer path analogue: `ecmult_gen_const`). |
| **`ecmult_const`** | Constant-time `q*A` (analogue: `ecmult_const`). |
| **`ecmult_0p_g`**, **`ecmult_1p_g`**, **`ecmult_multi_*`**, … | Verify / MSM-style workloads; blvm exposes `ecmult`, `ecmult_gen`, `ecmult_multi`. |

**Not identical:** libsecp uses its own table sizes and scalar streams; expect **ballpark** parity only.

## Iteration count

libsecp: **`SECP256K1_BENCH_ITERS`**. Large values make `bench_ecmult`’s MSM sweep slow; try `500`–`2000` for a quick pass.
