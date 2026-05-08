//! Build script for blvm-secp256k1.
//!
//! - Compiles field_10x26_arm.s on 32-bit ARM targets
//! - Compiles scalar_4x64_x86_64.s on x86_64
//!
//! Precomputed ecmult tables are committed in src/ecmult_precomputed.rs.
//! To regenerate: `cargo run --example regenerate_precomputed`

use std::env;
use std::path::PathBuf;

fn main() {
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    // Respect explicit opt-out (e.g. CI cross-compile to Windows/MinGW).
    let no_asm = env::var("BLVM_SECP256K1_NO_ASM").is_ok_and(|v| v == "1");

    if no_asm || target_os == "windows" {
        return;
    }

    if target_arch == "arm" {
        // Compile the ARM assembly for field mul/sqr.
        // Guard matches field/mod.rs: skip Windows ARM targets (no ELF assembler available).
        if target_os == "windows" {
            return;
        }
        let asm_path = manifest_dir.join("asm").join("field_10x26_arm.s");

        if asm_path.exists() {
            cc::Build::new()
                .file(&asm_path)
                .compile("blvm_secp256k1_field_asm");

            println!("cargo:rerun-if-changed={}", asm_path.display());
        } else {
            panic!(
                "ARM assembly not found at {}. Expected for target_arch=arm.",
                asm_path.display()
            );
        }
    }

    if target_arch == "x86_64" {
        // Compile x86_64 scalar assembly (from libsecp256k1 scalar_4x64_impl.h)
        let scalar_asm_path = manifest_dir.join("asm").join("scalar_4x64_x86_64.s");

        if scalar_asm_path.exists() {
            cc::Build::new()
                .file(&scalar_asm_path)
                .compile("blvm_secp256k1_scalar_asm");

            println!("cargo:rerun-if-changed={}", scalar_asm_path.display());
        }
    }

    println!("cargo:rerun-if-changed=build.rs");
}
