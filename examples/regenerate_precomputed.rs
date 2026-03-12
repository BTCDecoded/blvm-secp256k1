//! Regenerate src/ecmult_precomputed.rs from upstream libsecp256k1.
//!
//! Run: cargo run --example regenerate_precomputed
//!
//! First copy precomputed_ecmult.c from upstream libsecp256k1 src/ into
//! build/, then run this example.

use std::env;
use std::fs;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let c_path = manifest_dir.join("build").join("precomputed_ecmult.c");
    let out_path = manifest_dir.join("src").join("ecmult_precomputed.rs");

    let content = fs::read_to_string(&c_path).unwrap_or_else(|e| {
        eprintln!("Failed to read {}: {}", c_path.display(), e);
        eprintln!("Copy precomputed_ecmult.c from upstream libsecp256k1 src/");
        std::process::exit(1);
    });

    let (pre_g, pre_g_128) = parse_precomputed(&content);
    let rust_code = format_output(&pre_g, &pre_g_128);

    fs::write(&out_path, rust_code).unwrap_or_else(|e| {
        eprintln!("Failed to write {}: {}", out_path.display(), e);
        std::process::exit(1);
    });

    println!(
        "Wrote {} ({} PRE_G, {} PRE_G_128)",
        out_path.display(),
        pre_g.len(),
        pre_g_128.len()
    );
}

fn parse_precomputed(content: &str) -> (Vec<String>, Vec<String>) {
    let mut pre_g = Vec::new();
    let mut pre_g_128 = Vec::new();
    let mut in_pre_g_128 = false;

    for line in content.lines() {
        let line = line.trim();
        if line.contains("pre_g_128[") {
            in_pre_g_128 = true;
            continue;
        }
        let rest = if let Some(pos) = line.find("S(") {
            &line[pos + 2..]
        } else {
            continue;
        };
        let rest = rest.trim_end_matches(')');
        if !rest.is_empty() {
            let parts: Vec<&str> = rest.split(',').collect();
            if parts.len() != 16 {
                continue;
            }
            let vals: Vec<u32> = parts
                .iter()
                .map(|s| u32::from_str_radix(s.trim(), 16).unwrap_or(0))
                .collect();
            let x_n0 = (vals[7] as u64) | ((vals[6] as u64) << 32);
            let x_n1 = (vals[5] as u64) | ((vals[4] as u64) << 32);
            let x_n2 = (vals[3] as u64) | ((vals[2] as u64) << 32);
            let x_n3 = (vals[1] as u64) | ((vals[0] as u64) << 32);
            let y_n0 = (vals[15] as u64) | ((vals[14] as u64) << 32);
            let y_n1 = (vals[13] as u64) | ((vals[12] as u64) << 32);
            let y_n2 = (vals[11] as u64) | ((vals[10] as u64) << 32);
            let y_n3 = (vals[9] as u64) | ((vals[8] as u64) << 32);

            let entry = format!(
                "    GeStorage {{ x: FeStorage {{ n: [0x{:x}, 0x{:x}, 0x{:x}, 0x{:x}] }}, y: FeStorage {{ n: [0x{:x}, 0x{:x}, 0x{:x}, 0x{:x}] }} }}",
                x_n0, x_n1, x_n2, x_n3, y_n0, y_n1, y_n2, y_n3
            );
            if in_pre_g_128 {
                pre_g_128.push(entry);
            } else {
                pre_g.push(entry);
            }
        }
    }

    (pre_g, pre_g_128)
}

fn format_output(pre_g: &[String], pre_g_128: &[String]) -> String {
    format!(
        r#"// Precomputed ecmult tables for secp256k1 (G and 2^128*G).
// Committed directly; no build-time generation.
// To regenerate from upstream libsecp256k1: cargo run --example regenerate_precomputed

use crate::field::FeStorage;
use crate::group::GeStorage;

#[allow(clippy::large_const_arrays)]
pub const PRE_G: [GeStorage; {}] = [
{}
];

#[allow(clippy::large_const_arrays)]
pub const PRE_G_128: [GeStorage; {}] = [
{}
];
"#,
        pre_g.len(),
        pre_g.join(",\n"),
        pre_g_128.len(),
        pre_g_128.join(",\n")
    )
}
