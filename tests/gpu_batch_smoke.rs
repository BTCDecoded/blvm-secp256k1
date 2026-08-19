//! Smoke test for optional `gpu` feature (in-tree CUDA).
//! Skips when no GPU context can be created.

#![cfg(any(feature = "gpu", blvm_secp_gpu))]

use blvm_secp256k1::ecdsa::{
    ecdsa_sign_compact_rfc6979, ecdsa_verify_one_compact, ge_to_compressed, pubkey_from_secret,
};
use blvm_secp256k1::gpu::{
    gpu_available, try_ecdsa_verify_batch, try_schnorr_verify_batch, GPU_BATCH_MIN,
};
use blvm_secp256k1::scalar::Scalar;
use blvm_secp256k1::schnorr::{schnorr_sign, schnorr_verify, xonly_pubkey_from_secret};

fn scalar_from_b32(b: &[u8; 32]) -> Scalar {
    let mut s = Scalar::zero();
    s.set_b32(b);
    s
}

#[test]
fn gpu_ecdsa_batch_matches_cpu_when_available() {
    if !gpu_available() {
        eprintln!("skip: native CUDA context unavailable");
        return;
    }
    let n = GPU_BATCH_MIN.max(128);
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

    let gpu = try_ecdsa_verify_batch(&sigs, &msgs, &pks).expect("GPU batch should run");
    // Per-item CPU (not batch_results — that would re-enter GPU).
    let cpu: Vec<bool> = (0..n)
        .map(|i| ecdsa_verify_one_compact(&sigs[i], &msgs[i], &pks[i]))
        .collect();
    assert_eq!(gpu.len(), n);
    assert!(gpu.iter().all(|&ok| ok), "GPU rejected valid batch");
    assert_eq!(gpu, cpu);

    sigs[7][40] ^= 0x01;
    let gpu_mixed = try_ecdsa_verify_batch(&sigs, &msgs, &pks).expect("GPU mixed");
    assert!(!gpu_mixed[7], "GPU must reject mutated sig");
    assert_eq!(
        gpu_mixed.iter().filter(|&&ok| ok).count(),
        n - 1,
        "only one failure expected"
    );
}

#[test]
fn gpu_schnorr_batch_matches_cpu_when_available() {
    if !gpu_available() {
        eprintln!("skip: native CUDA context unavailable");
        return;
    }
    let n = GPU_BATCH_MIN.max(128);
    let mut sigs = Vec::with_capacity(n);
    let mut msgs = Vec::with_capacity(n);
    let mut pks = Vec::with_capacity(n);
    let aux = [0u8; 32];
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
        let pk = xonly_pubkey_from_secret(&sk).expect("xonly pk");
        let sig = schnorr_sign(&sk, &msg, &aux).expect("sign");
        sigs.push(sig);
        msgs.push(msg);
        pks.push(pk);
    }
    let msg_refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();

    let gpu = try_schnorr_verify_batch(&sigs, &msg_refs, &pks).expect("GPU Schnorr batch");
    // Per-item CPU (batch size 1 < GPU_BATCH_MIN → no GPU re-entry).
    let cpu: Vec<bool> = (0..n)
        .map(|i| schnorr_verify(&sigs[i], &msgs[i], &pks[i]))
        .collect();
    assert_eq!(gpu.len(), n);
    assert!(gpu.iter().all(|&ok| ok), "GPU rejected valid Schnorr batch");
    assert_eq!(gpu, cpu);

    sigs[11][40] ^= 0x01;
    let msg_refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
    let gpu_mixed = try_schnorr_verify_batch(&sigs, &msg_refs, &pks).expect("GPU Schnorr mixed");
    assert!(!gpu_mixed[11], "GPU must reject mutated Schnorr sig");
    assert_eq!(
        gpu_mixed.iter().filter(|&&ok| ok).count(),
        n - 1,
        "only one Schnorr failure expected"
    );
}

/// secp256k1 field prime p as big-endian bytes.
fn field_p_be() -> [u8; 32] {
    // FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFF FFFFFFFE FFFFFC2F
    [
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xff, 0xff,
        0xfc, 0x2f,
    ]
}

#[test]
fn gpu_rejects_field_overflow_encodings() {
    if !gpu_available() {
        eprintln!("skip: native CUDA context unavailable");
        return;
    }
    let n = GPU_BATCH_MIN.max(64);
    let p = field_p_be();
    let aux = [0u8; 32];

    // Schnorr: pad with valid, inject r=p and pk=p.
    let mut sigs = Vec::with_capacity(n);
    let mut msgs = Vec::with_capacity(n);
    let mut pks = Vec::with_capacity(n);
    for i in 0..n {
        let mut sk = [0u8; 32];
        sk[31] = (i as u8).saturating_add(1);
        let mut msg = [0u8; 32];
        msg[31] = i as u8;
        let pk = xonly_pubkey_from_secret(&sk).expect("xonly");
        let sig = schnorr_sign(&sk, &msg, &aux).expect("sign");
        sigs.push(sig);
        msgs.push(msg);
        pks.push(pk);
    }
    sigs[3][..32].copy_from_slice(&p);
    pks[5] = p;
    let msg_refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
    let gpu = try_schnorr_verify_batch(&sigs, &msg_refs, &pks).expect("GPU overflow batch");
    assert!(!gpu[3], "GPU must reject Schnorr r=p");
    assert!(!gpu[5], "GPU must reject Schnorr pk x=p");
    assert!(!schnorr_verify(&sigs[3], &msgs[3], &pks[3]));
    assert!(!schnorr_verify(&sigs[5], &msgs[5], &pks[5]));

    // ECDSA: inject compressed pk with x=p.
    let mut esigs = Vec::with_capacity(n);
    let mut emsgs = Vec::with_capacity(n);
    let mut epks = Vec::with_capacity(n);
    for i in 0..n {
        let mut sk = [0u8; 32];
        sk[31] = (i as u8).saturating_add(1);
        let mut msg = [0u8; 32];
        msg[31] = (i as u8).wrapping_mul(3);
        let pk = ge_to_compressed(&pubkey_from_secret(&scalar_from_b32(&sk)));
        let sig = ecdsa_sign_compact_rfc6979(&msg, &sk).expect("sign");
        esigs.push(sig);
        emsgs.push(msg);
        epks.push(pk);
    }
    let mut bad_pk = [0u8; 33];
    bad_pk[0] = 0x02;
    bad_pk[1..].copy_from_slice(&p);
    epks[9] = bad_pk;
    let gpu_e = try_ecdsa_verify_batch(&esigs, &emsgs, &epks).expect("GPU ECDSA overflow");
    assert!(!gpu_e[9], "GPU must reject ECDSA pk x=p");
    assert!(!ecdsa_verify_one_compact(&esigs[9], &emsgs[9], &epks[9]));
}
