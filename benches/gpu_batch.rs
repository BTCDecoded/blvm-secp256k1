//! CPU vs GPU batch verify (`cfg(blvm_secp_gpu)` / feature `gpu`).
//!
//! ```bash
//! cargo bench --features gpu --bench gpu_batch
//! # quick:
//! cargo bench --features gpu --bench gpu_batch -- --sample-size 20 --warm-up-time 1 --measurement-time 2
//! ```

#![cfg(any(feature = "gpu", blvm_secp_gpu))]

use std::hint::black_box;
use std::time::Instant;

use blvm_secp256k1::ecdsa::{
    ecdsa_sign_compact_rfc6979, ecdsa_verify_one_compact, ge_to_compressed, pubkey_from_secret,
};
use blvm_secp256k1::gpu::{
    gpu_available, try_ecdsa_verify_batch_ungated, try_schnorr_verify_batch_ungated,
};
use blvm_secp256k1::scalar::Scalar;
use blvm_secp256k1::schnorr::{schnorr_sign, schnorr_verify, xonly_pubkey_from_secret};
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

fn scalar_from_b32(b: &[u8; 32]) -> Scalar {
    let mut s = Scalar::zero();
    s.set_b32(b);
    s
}

struct EcdsaBatch {
    sigs: Vec<[u8; 64]>,
    msgs: Vec<[u8; 32]>,
    pks: Vec<[u8; 33]>,
}

struct SchnorrBatch {
    sigs: Vec<[u8; 64]>,
    msgs: Vec<[u8; 32]>,
    pks: Vec<[u8; 32]>,
}

fn make_ecdsa(n: usize) -> EcdsaBatch {
    let mut sigs = Vec::with_capacity(n);
    let mut msgs = Vec::with_capacity(n);
    let mut pks = Vec::with_capacity(n);
    for i in 0..n {
        let mut sk = [0u8; 32];
        sk[30] = (i / 256) as u8;
        sk[31] = (i % 256) as u8;
        if sk == [0u8; 32] {
            sk[31] = 1;
        }
        let mut msg = [0u8; 32];
        msg[30] = ((i * 3) / 256) as u8;
        msg[31] = ((i * 3) % 256) as u8;
        let pk = ge_to_compressed(&pubkey_from_secret(&scalar_from_b32(&sk)));
        let sig = ecdsa_sign_compact_rfc6979(&msg, &sk).expect("sign");
        sigs.push(sig);
        msgs.push(msg);
        pks.push(pk);
    }
    EcdsaBatch { sigs, msgs, pks }
}

fn make_schnorr(n: usize) -> SchnorrBatch {
    let aux = [0u8; 32];
    let mut sigs = Vec::with_capacity(n);
    let mut msgs = Vec::with_capacity(n);
    let mut pks = Vec::with_capacity(n);
    for i in 0..n {
        let mut sk = [0u8; 32];
        sk[30] = (i / 256) as u8;
        sk[31] = (i % 256) as u8;
        if sk == [0u8; 32] {
            sk[31] = 1;
        }
        let mut msg = [0u8; 32];
        msg[30] = ((i * 5) / 256) as u8;
        msg[31] = ((i * 5) % 256) as u8;
        let pk = xonly_pubkey_from_secret(&sk).expect("pk");
        let sig = schnorr_sign(&sk, &msg, &aux).expect("sign");
        sigs.push(sig);
        msgs.push(msg);
        pks.push(pk);
    }
    SchnorrBatch { sigs, msgs, pks }
}

fn cpu_ecdsa_per_item(b: &EcdsaBatch) -> usize {
    let mut ok = 0usize;
    for i in 0..b.sigs.len() {
        if ecdsa_verify_one_compact(&b.sigs[i], &b.msgs[i], &b.pks[i]) {
            ok += 1;
        }
    }
    ok
}

fn cpu_schnorr_per_item(b: &SchnorrBatch) -> usize {
    let mut ok = 0usize;
    for i in 0..b.sigs.len() {
        if schnorr_verify(&b.sigs[i], &b.msgs[i], &b.pks[i]) {
            ok += 1;
        }
    }
    ok
}

/// One-shot wall-clock summary printed before Criterion (handy for a quick look).
fn print_quick_summary(sizes: &[usize]) {
    if !gpu_available() {
        eprintln!("[gpu_batch] GPU unavailable — Criterion GPU arms will skip/fail.");
        return;
    }
    eprintln!("\n=== Quick CPU vs GPU batch verify (release) ===");
    eprintln!(
        "{:>8}  {:>12}  {:>12}  {:>8}  {:>12}  {:>12}  {:>8}",
        "n", "ecdsa_cpu", "ecdsa_gpu", "x", "schnorr_cpu", "schnorr_gpu", "x"
    );
    for &n in sizes {
        let e = make_ecdsa(n);
        let s = make_schnorr(n);
        let warmup = 2usize;
        let iters = 8usize;

        for _ in 0..warmup {
            black_box(cpu_ecdsa_per_item(&e));
            black_box(try_ecdsa_verify_batch_ungated(&e.sigs, &e.msgs, &e.pks));
            black_box(cpu_schnorr_per_item(&s));
            let refs: Vec<&[u8]> = s.msgs.iter().map(|m| m.as_slice()).collect();
            black_box(try_schnorr_verify_batch_ungated(&s.sigs, &refs, &s.pks));
        }

        let t0 = Instant::now();
        for _ in 0..iters {
            black_box(cpu_ecdsa_per_item(&e));
        }
        let e_cpu = t0.elapsed() / iters as u32;

        let t0 = Instant::now();
        for _ in 0..iters {
            black_box(try_ecdsa_verify_batch_ungated(&e.sigs, &e.msgs, &e.pks).expect("gpu ecdsa"));
        }
        let e_gpu = t0.elapsed() / iters as u32;

        let t0 = Instant::now();
        for _ in 0..iters {
            black_box(cpu_schnorr_per_item(&s));
        }
        let s_cpu = t0.elapsed() / iters as u32;

        let t0 = Instant::now();
        for _ in 0..iters {
            let refs: Vec<&[u8]> = s.msgs.iter().map(|m| m.as_slice()).collect();
            black_box(
                try_schnorr_verify_batch_ungated(&s.sigs, &refs, &s.pks).expect("gpu schnorr"),
            );
        }
        let s_gpu = t0.elapsed() / iters as u32;

        let e_x = e_cpu.as_secs_f64() / e_gpu.as_secs_f64().max(1e-12);
        let s_x = s_cpu.as_secs_f64() / s_gpu.as_secs_f64().max(1e-12);
        eprintln!(
            "{n:>8}  {:>10.2?}  {:>10.2?}  {e_x:>7.2}x  {:>10.2?}  {:>10.2?}  {s_x:>7.2}x",
            e_cpu, e_gpu, s_cpu, s_gpu
        );
        eprintln!(
            "         ({:.1} ksig/s cpu, {:.1} ksig/s gpu ecdsa | {:.1} / {:.1} schnorr)",
            n as f64 / e_cpu.as_secs_f64() / 1e3,
            n as f64 / e_gpu.as_secs_f64() / 1e3,
            n as f64 / s_cpu.as_secs_f64() / 1e3,
            n as f64 / s_gpu.as_secs_f64() / 1e3,
        );
    }
    eprintln!();
}

fn bench_ecdsa_batch(c: &mut Criterion) {
    if !gpu_available() {
        eprintln!("skip ecdsa GPU benches: no GPU");
        return;
    }
    let mut group = c.benchmark_group("ecdsa_batch_verify");
    for n in [64usize, 256, 1024, 4096] {
        let batch = make_ecdsa(n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("cpu_per_item", n), &batch, |b, batch| {
            b.iter(|| black_box(cpu_ecdsa_per_item(batch)))
        });
        group.bench_with_input(BenchmarkId::new("gpu", n), &batch, |b, batch| {
            b.iter(|| {
                black_box(
                    try_ecdsa_verify_batch_ungated(&batch.sigs, &batch.msgs, &batch.pks)
                        .expect("gpu"),
                )
            })
        });
    }
    group.finish();
}

fn bench_schnorr_batch(c: &mut Criterion) {
    if !gpu_available() {
        eprintln!("skip schnorr GPU benches: no GPU");
        return;
    }
    let mut group = c.benchmark_group("schnorr_batch_verify");
    for n in [64usize, 256, 1024, 4096] {
        let batch = make_schnorr(n);
        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::new("cpu_per_item", n), &batch, |b, batch| {
            b.iter(|| black_box(cpu_schnorr_per_item(batch)))
        });
        group.bench_with_input(BenchmarkId::new("gpu", n), &batch, |b, batch| {
            b.iter(|| {
                let refs: Vec<&[u8]> = batch.msgs.iter().map(|m| m.as_slice()).collect();
                black_box(
                    try_schnorr_verify_batch_ungated(&batch.sigs, &refs, &batch.pks).expect("gpu"),
                )
            })
        });
    }
    group.finish();
}

fn benches(c: &mut Criterion) {
    // Production BATCH_MIN stays 1024; ungated APIs measure GPU at small n too.
    print_quick_summary(&[64, 256, 1024, 4096, 16384]);
    bench_ecdsa_batch(c);
    bench_schnorr_batch(c);
}

criterion_group!(gpu_benches, benches);
criterion_main!(gpu_benches);
