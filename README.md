# blvm-secp256k1

[![crates.io](https://img.shields.io/crates/v/blvm-secp256k1.svg)](https://crates.io/crates/blvm-secp256k1)
[![docs.rs](https://docs.rs/blvm-secp256k1/badge.svg)](https://docs.rs/blvm-secp256k1)
[![CI](https://github.com/BTCDecoded/blvm-secp256k1/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/BTCDecoded/blvm-secp256k1/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE)

Pure Rust reimplementation of libsecp256k1 with vendored ASM for hot paths. No FFI to C — all logic in Rust. Performance-optimized secp256k1 for Bitcoin and Ethereum.

Precomputed ecmult tables are committed in `src/ecmult_precomputed.rs` (no build-time generation). To regenerate: copy `precomputed_ecmult.c` from upstream libsecp256k1 `src/` into `build/`, then `cargo run --example regenerate_precomputed`.

## Features

- **Context-free API** — Stateless, no `Secp256k1<C>` context. All functions take raw bytes, points, or scalars directly.
- **Pure Rust + ASM** — Field and scalar assembly for ARM32 and x86_64 where it matters; no C dependencies at runtime.
- **Parity with libsecp256k1** — ECDSA, Schnorr (BIP 340), MuSig2, ElligatorSwift (BIP 324), ECDH, Taproot.

## Side-channel and timing (constant-time) contract

Operations that use **private keys, nonces, or ECDH/MuSig secrets** are implemented with constant-time primitives on supported targets, aligned with [libsecp256k1](https://github.com/bitcoin-core/secp256k1) for comparable entry points. **Verification**, batch verify, and parsing paths may use faster variable-time math on public data only.

- **CT-required (secret inputs):** signing (`ecdsa_sig_sign*`, `schnorr_sign`), public key from secret, ECDH, ElligatorSwift `ellswift_xdh`, MuSig nonce generation / tweaks that multiply by a secret scalar, Taproot x-only tweak (`t*G`).

- **Var-time allowed (public inputs):** `ecdsa_sig_verify`, `schnorr_verify*`, `ecmult` / `ecmult_gen` / `ecmult_multi` when used only with public points and scalars from signatures or hashes.

- **Unsafe use:** never feed **secret** scalars into APIs documented as public-data / variable-time, or use this crate on architectures where `Scalar::inv` is not the constant-time path (see `Scalar::inv` in the `scalar` module).

The full **API matrix** and notes live in [TIMING.md](./TIMING.md) in this crate.

## Installation

```toml
[dependencies]
blvm-secp256k1 = "0.1"
```

## Usage

### ECDSA

```rust
use blvm_secp256k1::ecdsa::{ecdsa_sig_sign, ecdsa_sig_verify, pubkey_from_secret};
use blvm_secp256k1::scalar::Scalar;
use blvm_secp256k1::group::Gej;

let sk_bytes = [1u8; 32];
let msg_hash = [2u8; 32];

let mut pk_gej = Gej::default();
pubkey_from_secret(&sk_bytes, &mut pk_gej);

let sig = ecdsa_sig_sign(&msg_hash, &sk_bytes, None).expect("sign");
assert!(ecdsa_sig_verify(&msg_hash, &sig, &pk_gej));
```

### Schnorr (BIP 340)

```rust
use blvm_secp256k1::schnorr::{schnorr_sign, schnorr_verify, xonly_pubkey_from_secret};

let sk = [1u8; 32];
let msg = [2u8; 32];

let pk = xonly_pubkey_from_secret(&sk).expect("pk");
let sig = schnorr_sign(&msg, &sk, None).expect("sign");
assert!(schnorr_verify(&sig, &msg, &pk));
```

## License

MIT
