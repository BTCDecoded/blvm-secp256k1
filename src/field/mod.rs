//! Field element arithmetic for secp256k1.
//!
//! Field modulus: p = 2^256 - 2^32 - 977 (SEC2 secp256k1)
//!
//! Platform layout:
//! - x86_64, aarch64: 5x52 (u64 limbs), pure Rust with u128 wide mul
//! - arm (32-bit, non-Windows): 10x26 (u32 limbs), ASM for mul/sqr
//! - everything else (including Windows ARM): 5x52 falls back to no u128 path (not yet impl)

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
mod layout_5x52;

#[cfg(all(target_arch = "arm", not(target_os = "windows")))]
mod layout_10x26;

#[cfg(any(target_arch = "x86_64", target_arch = "aarch64"))]
pub use layout_5x52::{FeStorage, FieldElement};

#[cfg(all(target_arch = "arm", not(target_os = "windows")))]
pub use layout_10x26::FieldElement;
