//! Build script for blvm-secp256k1.
//!
//! - Compiles field_10x26_arm.s on 32-bit ARM targets
//! - Compiles scalar_4x64_x86_64.s on x86_64
//! - With feature `gpu`: compiles in-tree CUDA (`cuda/`) via nvcc → libblvm_secp_gpu
//!
//! Precomputed ecmult tables are committed in src/ecmult_precomputed.rs.
//! To regenerate: `cargo run --example regenerate_precomputed`

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    println!("cargo:rerun-if-env-changed=BLVM_SECP256K1_NO_ASM");
    println!("cargo:rerun-if-env-changed=CUDA_PATH");
    println!("cargo:rerun-if-env-changed=CUDA_LIB_DIR");
    println!("cargo:rerun-if-env-changed=BLVM_SECP256K1_CUDA_ARCH");

    if env::var("CARGO_FEATURE_GPU").is_ok() {
        build_cuda(&manifest_dir);
    }

    let no_asm = env::var("BLVM_SECP256K1_NO_ASM").is_ok_and(|v| v == "1");

    if no_asm || target_os == "windows" {
        println!("cargo:rerun-if-changed=build.rs");
        return;
    }

    if target_arch == "arm" {
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

fn find_nvcc() -> PathBuf {
    if let Ok(p) = env::var("CUDA_PATH") {
        let nvcc = PathBuf::from(&p).join("bin").join("nvcc");
        if nvcc.exists() {
            return nvcc;
        }
    }
    for root in ["/opt/cuda", "/usr/local/cuda"] {
        let nvcc = PathBuf::from(root).join("bin").join("nvcc");
        if nvcc.exists() {
            return nvcc;
        }
    }
    if let Ok(out) = Command::new("which").arg("nvcc").output() {
        if out.status.success() {
            let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
            if !s.is_empty() {
                return PathBuf::from(s);
            }
        }
    }
    panic!(
        "feature gpu: nvcc not found. Install CUDA toolkit or set CUDA_PATH \
         (e.g. CUDA_PATH=/opt/cuda). See cuda/README.md."
    );
}

fn build_cuda(manifest_dir: &PathBuf) {
    let nvcc = find_nvcc();
    let cuda_dir = manifest_dir.join("cuda");
    let src = cuda_dir.join("src").join("device.cu");
    let include = cuda_dir.join("include");
    println!("cargo:rerun-if-changed={}", src.display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_secp_gpu.h").display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_field.cuh").display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_scalar.cuh").display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_point.cuh").display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_ecmult.cuh").display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_verify.cuh").display());
    println!("cargo:rerun-if-changed={}", include.join("blvm_sha256.cuh").display());

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let obj = out_dir.join("blvm_secp_gpu.o");
    let lib = out_dir.join("libblvm_secp_gpu.a");

    let arch = env::var("BLVM_SECP256K1_CUDA_ARCH").unwrap_or_else(|_| "native".into());

    let status = Command::new(&nvcc)
        .args([
            "-c",
            "-Xcompiler",
            "-fPIC",
            "-O3",
            "-std=c++17",
            "-arch",
            &arch,
            "-I",
        ])
        .arg(&include)
        .arg(&src)
        .arg("-o")
        .arg(&obj)
        .status()
        .unwrap_or_else(|e| panic!("feature gpu: failed to spawn nvcc: {e}"));
    if !status.success() {
        panic!("feature gpu: nvcc failed compiling {}", src.display());
    }

    let ar = env::var("AR").unwrap_or_else(|_| "ar".into());
    let status = Command::new(&ar)
        .args(["crus"])
        .arg(&lib)
        .arg(&obj)
        .status()
        .unwrap_or_else(|e| panic!("feature gpu: ar failed: {e}"));
    if !status.success() {
        panic!("feature gpu: ar failed creating {}", lib.display());
    }

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=blvm_secp_gpu");

    if let Ok(cuda_lib) = env::var("CUDA_LIB_DIR") {
        println!("cargo:rustc-link-search=native={cuda_lib}");
    } else if PathBuf::from("/opt/cuda/lib64").exists() {
        println!("cargo:rustc-link-search=native=/opt/cuda/lib64");
    } else if PathBuf::from("/usr/local/cuda/lib64").exists() {
        println!("cargo:rustc-link-search=native=/usr/local/cuda/lib64");
    }
    println!("cargo:rustc-link-lib=dylib=cudart");
    println!("cargo:rustc-link-lib=dylib=stdc++");
}
