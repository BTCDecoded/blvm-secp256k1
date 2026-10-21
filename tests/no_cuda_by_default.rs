//! A default build must not link CUDA. `build.rs` sets `cfg(blvm_secp_gpu)` exactly
//! when it compiled `cuda/` and emitted `rustc-link-lib=dylib=cudart`, which puts
//! `libcudart` in the binary's NEEDED list. That must only happen with feature `gpu`.
//!
//! On a host with nvcc (CUDA_PATH set), this failed before CUDA became opt-in.

#[test]
fn cuda_links_only_with_feature_gpu() {
    assert_eq!(
        cfg!(blvm_secp_gpu),
        cfg!(feature = "gpu"),
        "cfg(blvm_secp_gpu) must match feature `gpu`"
    );
}
