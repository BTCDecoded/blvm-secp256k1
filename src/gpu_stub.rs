//! CPU-only stubs when CUDA was not linked (`cfg(not(blvm_secp_gpu))`).
//!
//! Same public surface as [`crate::gpu`] so callers can always try GPU and
//! fall back when this module returns `None` / `false`.

/// Default minimum batch size (matches the CUDA module).
pub const GPU_BATCH_MIN: usize = 1024;

/// Runtime batch-size floor.
pub fn gpu_batch_min() -> usize {
    GPU_BATCH_MIN
}

/// Context-pool size (unused without CUDA).
pub fn gpu_ctx_count() -> usize {
    1
}

/// Timer snapshot (always zeros without CUDA).
pub struct GpuTimerSnapshot {
    pub calls: u64,
    pub sigs: u64,
    pub lock_ns: u64,
    pub pack_ns: u64,
    pub kernel_ns: u64,
}

/// Snapshot GPU timer counters (zeros).
pub fn gpu_timer_snapshot() -> GpuTimerSnapshot {
    GpuTimerSnapshot {
        calls: 0,
        sigs: 0,
        lock_ns: 0,
        pack_ns: 0,
        kernel_ns: 0,
    }
}

/// Reset GPU timer counters (no-op).
pub fn gpu_timer_reset() {}

/// Enqueue an owned ECDSA job. Always `None` without CUDA.
pub fn enqueue_ecdsa_job(
    _msgs: Vec<[u8; 32]>,
    _pubkeys: Vec<[u8; 33]>,
    _sigs: Vec<[u8; 64]>,
) -> Option<std::sync::mpsc::Receiver<Option<Vec<bool>>>> {
    None
}

/// Per-item ECDSA verify on GPU. Always `None` without CUDA.
pub fn try_ecdsa_verify_batch(
    _sigs: &[[u8; 64]],
    _msgs: &[[u8; 32]],
    _pubkeys: &[[u8; 33]],
) -> Option<Vec<bool>> {
    None
}

/// Ungated ECDSA GPU path (benches). Always `None` without CUDA.
#[doc(hidden)]
pub fn try_ecdsa_verify_batch_ungated(
    _sigs: &[[u8; 64]],
    _msgs: &[[u8; 32]],
    _pubkeys: &[[u8; 33]],
) -> Option<Vec<bool>> {
    None
}

/// Per-item Schnorr verify on GPU. Always `None` without CUDA.
pub fn try_schnorr_verify_batch(
    _sigs: &[[u8; 64]],
    _msgs: &[&[u8]],
    _pubkeys: &[[u8; 32]],
) -> Option<Vec<bool>> {
    None
}

/// Ungated Schnorr GPU path (benches). Always `None` without CUDA.
#[doc(hidden)]
pub fn try_schnorr_verify_batch_ungated(
    _sigs: &[[u8; 64]],
    _msgs: &[&[u8]],
    _pubkeys: &[[u8; 32]],
) -> Option<Vec<bool>> {
    None
}

/// True if at least one CUDA GPU context is live.
pub fn gpu_available() -> bool {
    false
}
