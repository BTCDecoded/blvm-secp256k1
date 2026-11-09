//! A batch large enough for the GPU path must still verify on the CPU when GPU
//! use is turned off. Own process so the env is set before the GPU pool is built.

use blvm_secp256k1::ecdsa::{
    ecdsa_sign_compact_rfc6979, ecdsa_verify_batch_results, ge_to_compressed, pubkey_from_secret,
};
use blvm_secp256k1::gpu::{GPU_BATCH_MIN, try_ecdsa_verify_batch};
use blvm_secp256k1::scalar::Scalar;

fn scalar_from_b32(b: &[u8; 32]) -> Scalar {
    let mut s = Scalar::zero();
    s.set_b32(b);
    s
}

#[test]
fn large_batch_verifies_on_cpu_when_gpu_is_off() {
    // SAFETY: set_var is safe on the 1.85 toolchain and unsafe on newer ones.
    // This test binary does not read the variable before this call.
    unsafe { std::env::set_var("BLVM_SECP_GPU", "0") };

    let n = GPU_BATCH_MIN;
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
        msg[31] = (i % 251) as u8;
        let pk = ge_to_compressed(&pubkey_from_secret(&scalar_from_b32(&sk)));
        let sig = ecdsa_sign_compact_rfc6979(&msg, &sk).expect("sign");
        sigs.push(sig);
        msgs.push(msg);
        pks.push(pk);
    }

    assert!(try_ecdsa_verify_batch(&sigs, &msgs, &pks).is_none());
    assert_eq!(
        ecdsa_verify_batch_results(&sigs, &msgs, &pks),
        vec![true; n]
    );

    sigs[3][40] ^= 0x01;
    let mixed = ecdsa_verify_batch_results(&sigs, &msgs, &pks);
    assert!(!mixed[3]);
    assert_eq!(mixed.iter().filter(|ok| **ok).count(), n - 1);
}
