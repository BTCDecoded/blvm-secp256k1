# blvm-secp256k1 CUDA

BLVM-owned CUDA batch ECDSA / BIP-340 Schnorr verify. No UltrafastSecp256k1 / `libufsecp`.

Default `cargo build` compiles this tree when `nvcc` is found (`cfg(has_cuda_kernels)`). The binary does not link `libcudart`. At runtime it `dlopen`s the runtime and uses the GPU if `cuInit` sees a device; otherwise it stays on the CPU path. Missing `nvcc` skips the kernels and still builds. Build opt-out: `BLVM_SECP256K1_NO_GPU=1`.

## Build

```bash
# Auto when nvcc is on PATH / CUDA_PATH / /opt/cuda (Zeus)
cargo test
# CPU-only even when nvcc is installed:
BLVM_SECP256K1_NO_GPU=1 cargo test
```

Optional env:

| Var | Meaning |
|-----|---------|
| `CUDA_PATH` | Toolkit root (default `/opt/cuda` or `/usr/local/cuda`) |
| `BLVM_SECP256K1_NO_GPU=1` | Do not compile kernels |
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
