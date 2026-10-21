# blvm-secp256k1 CUDA

BLVM-owned CUDA batch ECDSA / BIP-340 Schnorr verify. No UltrafastSecp256k1 / `libufsecp`.

Only feature `gpu` compiles and links this tree; a default `cargo build` never does, even when `nvcc` is installed, so default binaries do not need `libcudart`. With `gpu`, the build fails if nvcc is missing, and the binary needs the CUDA runtime to load. Runtime uses GPU if a device inits, else CPU.

## Build

```bash
# CPU only (default; no CUDA dependency):
cargo test
# CUDA (nvcc on PATH / CUDA_PATH / /opt/cuda; fails if missing):
cargo test --features gpu
```

Optional env:

| Var | Meaning |
|-----|---------|
| `CUDA_PATH` | Toolkit root (default `/opt/cuda` or `/usr/local/cuda`) |
| `CUDA_LIB_DIR` | `cudart` search path |
| `BLVM_SECP256K1_CUDA_ARCH` | `nvcc -arch=` (default `native`) |

## Runtime (`BLVM_SECP_GPU*`)

Same knobs as Rust `gpu` module: opt-out `BLVM_SECP_GPU=0`, `BLVM_SECP_GPU_BATCH_MIN` (default **1024**), `BLVM_SECP_GPU_CTXS`, etc.

| Extra | Meaning |
|-------|---------|
| `BLVM_SECP_GPU_HOST_BRIDGE=1` | Opt-in debug: CUDA launch still runs, but per-item results are filled by Rust CPU verify (ignore device verdicts). |

**Default:** after CUDA ctx create, an **init KAT** (known ECDSA + BIP340 accept, high-S reject, `r≥p` reject) must pass; otherwise the pool is not created (CPU only). Wire encodings use `fe_set_b32_limit` (reject `x/r ≥ p`). GPU path is public verify only — no secrets.

## Bench

```bash
cargo bench --features gpu --bench gpu_batch
# prints a quick CPU vs GPU table (ungated), then Criterion arms
```

Launch path: skip redundant `cudaSetDevice`, non-blocking stream, 128-thread blocks (heavy EC kernels).
