//! commons-secp256k1: Rust reimplementation of libsecp256k1
//!
//! Pure Rust implementation with vendored ASM for hot paths (ARM32 field, x86_64 scalar).
//! No FFI to C — all logic in Rust.
//!
//! ## Context-free API
//!
//! This crate is **stateless and context-free**. Unlike rust-secp256k1, there is no
//! `Secp256k1<C>` context object. All functions take raw bytes, points, or scalars
//! directly. No global mutable state. MuSig uses a `Session` struct (created from
//! aggnonce + msg + keyagg_cache) — this is algorithm state, not a global context.

pub mod ecdh;
pub mod ecdsa;
pub mod ecmult;
mod ecmult_const;
mod ecmult_gen_comb;
pub use ecmult_const::ecmult_const;
pub use ecmult_gen_comb::{
    ecmult_gen_const_comb as ecmult_gen_const, precompute as ecmult_gen_precompute,
};
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
pub mod ellswift;
pub mod field;
pub mod group;
#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
pub(crate) mod modinv64;
pub mod musig;
mod rfc6979;
pub mod scalar;
pub mod schnorr;
pub mod schnorr_half_agg;
pub mod taproot;
