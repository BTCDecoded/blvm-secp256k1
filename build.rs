//! Build script for blvm-secp256k1.
//!
//! - Compiles field_10x26_arm.s on 32-bit ARM targets
//! - Compiles scalar_4x64_x86_64.s on x86_64
//! - Compiles in-tree CUDA (`cuda/`) when nvcc is found → libblvm_secp_gpu
//!   (`cfg(has_cuda_kernels)`). Does not link libcudart. `BLVM_SECP256K1_NO_GPU=1`
//!   skips the compile. Missing nvcc skips it too. Neither case fails the build.
//!
//! Precomputed ecmult tables are committed in src/ecmult_precomputed.rs.
//! To regenerate: `cargo run --example regenerate_precomputed`

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

fn main() {
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());

    println!("cargo:rerun-if-env-changed=BLVM_SECP256K1_NO_ASM");
    println!("cargo:rerun-if-env-changed=BLVM_SECP256K1_NO_GPU");
    println!("cargo:rerun-if-env-changed=CUDA_PATH");
    println!("cargo:rerun-if-env-changed=CUDA_LIB_DIR");
    println!("cargo:rerun-if-env-changed=BLVM_SECP256K1_CUDA_ARCH");
    println!("cargo:rustc-check-cfg=cfg(has_cuda_kernels)");
    println!("cargo:rustc-check-cfg=cfg(blvm_secp_gpu)");

    try_enable_cuda(&manifest_dir);

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

fn gpu_build_disabled() -> bool {
    env::var("BLVM_SECP256K1_NO_GPU").is_ok_and(|v| v == "1")
}

fn try_enable_cuda(manifest_dir: &Path) {
    if gpu_build_disabled() {
        return;
    }

    // Host nvcc does not produce objects for a cross target. Those binaries
    // stay CPU-only. Native Linux is where the dynamic loader has ELF stubs.
    let host = env::var("HOST").unwrap_or_default();
    let target = env::var("TARGET").unwrap_or_default();
    if !host.is_empty() && host != target {
        return;
    }
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    if target_os != "linux" || !matches!(target_arch.as_str(), "x86_64" | "aarch64") {
        return;
    }

    let Some(nvcc) = find_nvcc() else {
        return;
    };

    match build_cuda(manifest_dir, &nvcc) {
        Ok(()) => {
            println!("cargo:rustc-cfg=has_cuda_kernels");
            // Same condition, previous name. In-tree cfgs use has_cuda_kernels.
            println!("cargo:rustc-cfg=blvm_secp_gpu");
        }
        Err(e) => {
            println!("cargo:warning=blvm-secp256k1: CUDA skipped ({e}); CPU verify only");
        }
    }
}

fn find_nvcc() -> Option<PathBuf> {
    if let Ok(p) = env::var("CUDA_PATH") {
        let nvcc = PathBuf::from(&p).join("bin").join("nvcc");
        if nvcc.exists() {
            return Some(nvcc);
        }
    }
    for root in ["/opt/cuda", "/usr/local/cuda"] {
        let nvcc = PathBuf::from(root).join("bin").join("nvcc");
        if nvcc.exists() {
            return Some(nvcc);
        }
    }
    if let Ok(out) = Command::new("which").arg("nvcc").output() {
        if out.status.success() {
            let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
            if !s.is_empty() {
                return Some(PathBuf::from(s));
            }
        }
    }
    None
}

fn build_cuda(manifest_dir: &Path, nvcc: &Path) -> Result<(), String> {
    let cuda_dir = manifest_dir.join("cuda");
    let src = cuda_dir.join("src").join("device.cu");
    let include = cuda_dir.join("include");
    println!("cargo:rerun-if-changed={}", src.display());
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_secp_gpu.h").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_field.cuh").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_scalar.cuh").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_point.cuh").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_ecmult.cuh").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_verify.cuh").display()
    );
    println!(
        "cargo:rerun-if-changed={}",
        include.join("blvm_sha256.cuh").display()
    );

    let out_dir = PathBuf::from(env::var("OUT_DIR").map_err(|e| e.to_string())?);
    let obj = out_dir.join("blvm_secp_gpu.o");
    let lib = out_dir.join("libblvm_secp_gpu.a");

    let arch = env::var("BLVM_SECP256K1_CUDA_ARCH").unwrap_or_else(|_| "native".into());

    let status = Command::new(nvcc)
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
        .map_err(|e| format!("failed to spawn nvcc: {e}"))?;
    if !status.success() {
        return Err(format!("nvcc failed compiling {}", src.display()));
    }

    let ar = env::var("AR").unwrap_or_else(|_| "ar".into());
    let status = Command::new(&ar)
        .args(["crus"])
        .arg(&lib)
        .arg(&obj)
        .status()
        .map_err(|e| format!("ar failed: {e}"))?;
    if !status.success() {
        return Err(format!("ar failed creating {}", lib.display()));
    }

    let symbols = cuda_runtime_symbols(&obj)?;
    write_cudart_stubs(&out_dir, &target_arch_for_stubs()?, &symbols)?;
    compile_stubs(&out_dir)?;

    let cc_obj = out_dir.join("blvm_cudart_stub.o");
    let asm_obj = out_dir.join("blvm_cudart_stub.s.o");
    let status = Command::new(&ar)
        .args(["crus"])
        .arg(&lib)
        .arg(&obj)
        .arg(&cc_obj)
        .arg(&asm_obj)
        .status()
        .map_err(|e| format!("ar failed: {e}"))?;
    if !status.success() {
        return Err(format!("ar failed adding stubs to {}", lib.display()));
    }

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=blvm_secp_gpu");
    // C++ kernel host code. Not libcudart.
    println!("cargo:rustc-link-lib=dylib=stdc++");
    println!("cargo:rustc-link-lib=dylib=dl");
    Ok(())
}

fn target_arch_for_stubs() -> Result<String, String> {
    let arch = env::var("CARGO_CFG_TARGET_ARCH").unwrap_or_default();
    match arch.as_str() {
        "x86_64" | "aarch64" => Ok(arch),
        other => Err(format!("no libcudart stubs for arch {other}")),
    }
}

fn cuda_runtime_symbols(obj: &Path) -> Result<Vec<String>, String> {
    let out = Command::new("nm")
        .arg("-u")
        .arg(obj)
        .output()
        .map_err(|e| format!("nm failed: {e}"))?;
    if !out.status.success() {
        return Err("nm -u failed on the CUDA object".into());
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let mut syms = Vec::new();
    for line in text.lines() {
        let mut parts = line.split_whitespace();
        let Some(first) = parts.next() else {
            continue;
        };
        let name = if first == "U" {
            parts.next().unwrap_or("")
        } else {
            ""
        };
        if name.contains("cuda") && !name.contains('@') {
            syms.push(name.to_string());
        }
    }
    syms.sort();
    syms.dedup();
    if syms.is_empty() {
        return Err("CUDA object has no runtime symbols to stub".into());
    }
    Ok(syms)
}

fn write_cudart_stubs(out_dir: &Path, arch: &str, symbols: &[String]) -> Result<(), String> {
    let mut c = String::from(
        r#"#define _GNU_SOURCE
#include <stddef.h>
#include <dlfcn.h>

/* Filled by the constructor. Null means the stub returns 0 and does not call. */
"#,
    );
    for sym in symbols {
        c.push_str(&format!("void *blvm_cudart_{sym};\n"));
    }
    c.push_str(
        r#"
static void *open_cudart(void) {
    const char *names[] = {
        "libcudart.so",
        "libcudart.so.13",
        "libcudart.so.12",
        "libcudart.so.11.0",
        NULL
    };
    for (int i = 0; names[i] != NULL; i++) {
        void *handle = dlopen(names[i], RTLD_NOW | RTLD_GLOBAL);
        if (handle) return handle;
    }
    return NULL;
}

/* Earlier than the default fatbin constructor (priority 65535). */
static void blvm_bind_cudart(void) __attribute__((constructor(101)));
static void blvm_bind_cudart(void) {
    void *handle = open_cudart();
    if (!handle) return;
"#,
    );
    for sym in symbols {
        c.push_str(&format!(
            "    blvm_cudart_{sym} = dlsym(handle, \"{sym}\");\n"
        ));
    }
    c.push_str("}\n");
    let c_path = out_dir.join("blvm_cudart_stub.c");
    std::fs::write(&c_path, c).map_err(|e| format!("write stub c: {e}"))?;

    let mut asm = String::from("/* Jump to the real libcudart symbol, or return 0. */\n.text\n");
    for sym in symbols {
        if arch == "x86_64" {
            asm.push_str(&format!(
                r#"
    .globl {sym}
    .hidden {sym}
    .type {sym}, @function
{sym}:
    movq blvm_cudart_{sym}(%rip), %r11
    testq %r11, %r11
    jz 1f
    jmp *%r11
1:
    xorl %eax, %eax
    ret
    .size {sym}, .-{sym}
"#
            ));
        } else {
            asm.push_str(&format!(
                r#"
    .globl {sym}
    .hidden {sym}
    .type {sym}, %function
{sym}:
    adrp x16, blvm_cudart_{sym}
    ldr x16, [x16, :lo12:blvm_cudart_{sym}]
    cbz x16, 1f
    br x16
1:
    mov w0, #0
    ret
    .size {sym}, .-{sym}
"#
            ));
        }
    }
    let s_path = out_dir.join("blvm_cudart_stub.s");
    std::fs::write(&s_path, asm).map_err(|e| format!("write stub asm: {e}"))?;
    Ok(())
}

fn compile_stubs(out_dir: &Path) -> Result<(), String> {
    let cc_bin = env::var("CC").unwrap_or_else(|_| "cc".into());
    let c_obj = out_dir.join("blvm_cudart_stub.o");
    let s_obj = out_dir.join("blvm_cudart_stub.s.o");
    run_cc(
        Command::new(&cc_bin)
            .args(["-c", "-fPIC", "-O2", "-fvisibility=hidden"])
            .arg(out_dir.join("blvm_cudart_stub.c"))
            .arg("-o")
            .arg(&c_obj),
    )?;
    run_cc(
        Command::new(&cc_bin)
            .args(["-c", "-fPIC"])
            .arg(out_dir.join("blvm_cudart_stub.s"))
            .arg("-o")
            .arg(&s_obj),
    )?;
    Ok(())
}

fn run_cc(cmd: &mut Command) -> Result<(), String> {
    let out = cmd
        .output()
        .map_err(|e| format!("failed to spawn compiler: {e}"))?;
    if out.status.success() {
        return Ok(());
    }
    let err = String::from_utf8_lossy(&out.stderr);
    let tail: Vec<&str> = err.lines().rev().take(8).collect();
    Err(tail.into_iter().rev().collect::<Vec<_>>().join(" | "))
}
