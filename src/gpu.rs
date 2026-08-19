//! In-tree CUDA GPU batch verify (`cfg(blvm_secp_gpu)`; auto when nvcc is found).
//!
//! Links static `libblvm_secp_gpu` built from `cuda/` via nvcc (see `build.rs`).
//! No UltrafastSecp256k1 / `libufsecp`.
//!
//! ## Runtime policy
//! - **Default on** when CUDA + ctx create succeed.
//! - Opt out: `BLVM_SECP_GPU=0` (also `false` / `off` / `no`).
//! - Batch threshold: `BLVM_SECP_GPU_BATCH_MIN` (default **1024**; below this GPU loses to CPU).
//! - Context pool: `BLVM_SECP_GPU_CTXS` (default **1**).
//! - `BLVM_SECP_GPU_TRY_LOCK=1`: busy slots → `None` (CPU fallback).
//! - `BLVM_SECP_GPU_TIMERS=1`: accumulate lock/pack/kernel ns; log periodically.
//! - `BLVM_SECP_GPU_SUBMITTERS=1`: dedicated submitter threads own ctx mutexes.
//! - `BLVM_SECP_GPU_HOST_BRIDGE=1`: ignore device verdicts; fill via CPU (debug only).
//! - Pool create runs a device init KAT; failure → no GPU pool (CPU only).
//! - Compact SoA (`[[u8;32]]` / `33` / `64`) passed as contiguous pointers.

use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard, OnceLock};
use std::time::Instant;

use std::os::raw::{c_char, c_int, c_void};

const BLVM_SECP_GPU_OK: c_int = 0;

/// Default minimum batch size for GPU launch (overridable via `BLVM_SECP_GPU_BATCH_MIN`).
/// Tuned above the CPU↔GPU crossover (~1k sigs on current Zeus benches).
pub const GPU_BATCH_MIN: usize = 1024;

type GpuCtx = c_void;

#[link(name = "blvm_secp_gpu", kind = "static")]
unsafe extern "C" {
    fn blvm_secp_gpu_is_available() -> c_int;
    fn blvm_secp_gpu_device_count() -> u32;
    fn blvm_secp_gpu_ctx_create(ctx_out: *mut *mut GpuCtx, device_index: u32) -> c_int;
    fn blvm_secp_gpu_ctx_destroy(ctx: *mut GpuCtx);
    fn blvm_secp_gpu_ctx_set_pre_g(
        ctx: *mut GpuCtx,
        pre_g: *const u8,
        pre_g_128: *const u8,
        n_entries: usize,
    ) -> c_int;
    fn blvm_secp_gpu_ecdsa_verify_batch(
        ctx: *mut GpuCtx,
        msg_hashes32: *const u8,
        pubkeys33: *const u8,
        sigs64: *const u8,
        count: usize,
        out_results: *mut u8,
    ) -> c_int;
    fn blvm_secp_gpu_schnorr_verify_batch(
        ctx: *mut GpuCtx,
        msg_hashes32: *const u8,
        pubkeys_x32: *const u8,
        sigs64: *const u8,
        count: usize,
        out_results: *mut u8,
    ) -> c_int;
    #[allow(dead_code)]
    fn blvm_secp_gpu_error_str(err: c_int) -> *const c_char;
}

struct GpuState {
    ctx: *mut GpuCtx,
    device_index: u32,
}

// SAFETY: each gpu_ctx is single-thread; we serialize per slot via Mutex.
unsafe impl Send for GpuState {}

impl Drop for GpuState {
    fn drop(&mut self) {
        if !self.ctx.is_null() {
            unsafe { blvm_secp_gpu_ctx_destroy(self.ctx) };
            self.ctx = std::ptr::null_mut();
        }
    }
}

struct GpuPool {
    slots: Vec<Mutex<Option<GpuState>>>,
    rr: AtomicUsize,
}

fn env_disabled(name: &str) -> bool {
    match std::env::var(name) {
        Ok(v) => matches!(
            v.trim().to_ascii_lowercase().as_str(),
            "0" | "false" | "off" | "no"
        ),
        Err(_) => false,
    }
}

fn env_truthy(name: &str) -> bool {
    matches!(
        std::env::var(name)
            .ok()
            .as_deref()
            .map(str::trim)
            .map(str::to_ascii_lowercase)
            .as_deref(),
        Some("1") | Some("true") | Some("yes") | Some("on")
    )
}

/// Runtime batch-size floor (default [`GPU_BATCH_MIN`]).
pub fn gpu_batch_min() -> usize {
    static MIN: OnceLock<usize> = OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("BLVM_SECP_GPU_BATCH_MIN")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|&n| n >= 1)
            .unwrap_or(GPU_BATCH_MIN)
    })
}

/// Number of GPU contexts in the pool (default 1, clamp 1..=16).
pub fn gpu_ctx_count() -> usize {
    static N: OnceLock<usize> = OnceLock::new();
    *N.get_or_init(|| {
        std::env::var("BLVM_SECP_GPU_CTXS")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .unwrap_or(1)
            .clamp(1, 16)
    })
}

fn gpu_try_lock_from_env() -> bool {
    env_truthy("BLVM_SECP_GPU_TRY_LOCK")
}

fn gpu_timers_enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| env_truthy("BLVM_SECP_GPU_TIMERS"))
}

fn gpu_submitters_enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| env_truthy("BLVM_SECP_GPU_SUBMITTERS"))
}

/// Cheap Phase-0 / soak GPU timers (opt-in via `BLVM_SECP_GPU_TIMERS=1`).
pub struct GpuTimerSnapshot {
    pub calls: u64,
    pub sigs: u64,
    pub lock_ns: u64,
    pub pack_ns: u64,
    pub kernel_ns: u64,
}

struct GpuTimers {
    calls: AtomicU64,
    sigs: AtomicU64,
    lock_ns: AtomicU64,
    pack_ns: AtomicU64,
    kernel_ns: AtomicU64,
    last_log_calls: AtomicU64,
}

impl GpuTimers {
    const fn new() -> Self {
        Self {
            calls: AtomicU64::new(0),
            sigs: AtomicU64::new(0),
            lock_ns: AtomicU64::new(0),
            pack_ns: AtomicU64::new(0),
            kernel_ns: AtomicU64::new(0),
            last_log_calls: AtomicU64::new(0),
        }
    }

    fn record(&self, n: usize, lock_ns: u64, pack_ns: u64, kernel_ns: u64) {
        self.calls.fetch_add(1, Ordering::Relaxed);
        self.sigs.fetch_add(n as u64, Ordering::Relaxed);
        self.lock_ns.fetch_add(lock_ns, Ordering::Relaxed);
        self.pack_ns.fetch_add(pack_ns, Ordering::Relaxed);
        self.kernel_ns.fetch_add(kernel_ns, Ordering::Relaxed);
        let calls = self.calls.load(Ordering::Relaxed);
        let last = self.last_log_calls.load(Ordering::Relaxed);
        if calls == 1 || calls.saturating_sub(last) >= 64 {
            if self
                .last_log_calls
                .compare_exchange(last, calls, Ordering::Relaxed, Ordering::Relaxed)
                .is_ok()
            {
                let snap = self.snapshot();
                let total = snap.lock_ns + snap.pack_ns + snap.kernel_ns;
                let ms = |ns: u64| ns as f64 / 1_000_000.0;
                eprintln!(
                    "[BLVM_SECP_GPU_TIMERS] calls={} sigs={} lock_ms={:.1} pack_ms={:.1} \
                     kernel_ms={:.1} total_ms={:.1} avg_batch={:.0}",
                    snap.calls,
                    snap.sigs,
                    ms(snap.lock_ns),
                    ms(snap.pack_ns),
                    ms(snap.kernel_ns),
                    ms(total),
                    if snap.calls > 0 {
                        snap.sigs as f64 / snap.calls as f64
                    } else {
                        0.0
                    }
                );
            }
        }
    }

    fn snapshot(&self) -> GpuTimerSnapshot {
        GpuTimerSnapshot {
            calls: self.calls.load(Ordering::Relaxed),
            sigs: self.sigs.load(Ordering::Relaxed),
            lock_ns: self.lock_ns.load(Ordering::Relaxed),
            pack_ns: self.pack_ns.load(Ordering::Relaxed),
            kernel_ns: self.kernel_ns.load(Ordering::Relaxed),
        }
    }
}

static GPU_TIMERS: GpuTimers = GpuTimers::new();

/// Snapshot GPU timer counters (zeros if timers never enabled / no calls).
pub fn gpu_timer_snapshot() -> GpuTimerSnapshot {
    GPU_TIMERS.snapshot()
}

/// Reset GPU timer counters (for A/B cell isolation).
pub fn gpu_timer_reset() {
    GPU_TIMERS.calls.store(0, Ordering::Relaxed);
    GPU_TIMERS.sigs.store(0, Ordering::Relaxed);
    GPU_TIMERS.lock_ns.store(0, Ordering::Relaxed);
    GPU_TIMERS.pack_ns.store(0, Ordering::Relaxed);
    GPU_TIMERS.kernel_ns.store(0, Ordering::Relaxed);
    GPU_TIMERS.last_log_calls.store(0, Ordering::Relaxed);
}

fn upload_pre_g(ctx: *mut GpuCtx) -> Result<(), c_int> {
    let (pre_g, pre_g_128) = crate::ecmult::pre_g_table_bytes();
    const N: usize = 8192;
    debug_assert_eq!(pre_g.len(), N * 64);
    debug_assert_eq!(pre_g_128.len(), N * 64);
    let rc = unsafe {
        blvm_secp_gpu_ctx_set_pre_g(ctx, pre_g.as_ptr(), pre_g_128.as_ptr(), N)
    };
    if rc != BLVM_SECP_GPU_OK {
        return Err(rc);
    }
    Ok(())
}

fn try_create_one(device_index: u32) -> Option<GpuState> {
    unsafe {
        let mut ctx: *mut GpuCtx = std::ptr::null_mut();
        let rc = blvm_secp_gpu_ctx_create(&mut ctx, device_index);
        if rc != BLVM_SECP_GPU_OK || ctx.is_null() {
            eprintln!(
                "[BLVM_SECP_GPU] native CUDA ctx create failed device={device_index} rc={rc}"
            );
            return None;
        }
        if let Err(rc) = upload_pre_g(ctx) {
            eprintln!(
                "[BLVM_SECP_GPU] PRE_G upload failed device={device_index} rc={rc}"
            );
            blvm_secp_gpu_ctx_destroy(ctx);
            return None;
        }
        Some(GpuState { ctx, device_index })
    }
}

/// Device self-test before trusting kernel verdicts. Fail-closed → no GPU pool.
fn device_init_kat(ctx: *mut GpuCtx) -> Result<(), &'static str> {
    // ECDSA: RFC6979 sk=1 / msg=...07 (accept) + high-S (reject).
    let ecdsa_msg: [u8; 32] = [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 7,
    ];
    let ecdsa_pk: [u8; 33] = [
        0x02, 0x79, 0xbe, 0x66, 0x7e, 0xf9, 0xdc, 0xbb, 0xac, 0x55, 0xa0, 0x62, 0x95, 0xce, 0x87,
        0x0b, 0x07, 0x02, 0x9b, 0xfc, 0xdb, 0x2d, 0xce, 0x28, 0xd9, 0x59, 0xf2, 0x81, 0x5b, 0x16,
        0xf8, 0x17, 0x98,
    ];
    let ecdsa_sig: [u8; 64] = [
        0xcb, 0x65, 0x47, 0xdc, 0x10, 0xa2, 0x6c, 0x2b, 0xb5, 0xc8, 0xbc, 0xac, 0xdd, 0x86, 0x9a,
        0xa1, 0x04, 0x7b, 0x6f, 0x6a, 0x55, 0xa2, 0x32, 0x00, 0x3b, 0x53, 0x53, 0x38, 0x28, 0xcb,
        0xa1, 0x03, 0x1a, 0x95, 0xcd, 0xa9, 0xa2, 0xcd, 0xd0, 0x5b, 0xbd, 0x61, 0x15, 0x09, 0x74,
        0x58, 0xa7, 0x41, 0xaa, 0xf0, 0xb0, 0x37, 0x58, 0xfc, 0x2a, 0x5e, 0x4d, 0x9f, 0xe5, 0x39,
        0xc8, 0x62, 0x1e, 0x6f,
    ];
    let mut high_s_sig = ecdsa_sig;
    {
        let mut s = crate::scalar::Scalar::zero();
        let mut s_bytes = [0u8; 32];
        s_bytes.copy_from_slice(&ecdsa_sig[32..]);
        s.set_b32(&s_bytes);
        if s.is_high() {
            return Err("ECDSA init KAT fixture unexpectedly high-S");
        }
        let mut sn = crate::scalar::Scalar::zero();
        sn.negate(&s);
        if !sn.is_high() {
            return Err("ECDSA init KAT failed to build high-S");
        }
        sn.get_b32((&mut high_s_sig[32..]).try_into().unwrap());
    }

    let msgs = [ecdsa_msg, ecdsa_msg];
    let pks = [ecdsa_pk, ecdsa_pk];
    let sigs = [ecdsa_sig, high_s_sig];
    let mut out = [0u8; 2];
    let rc = unsafe {
        blvm_secp_gpu_ecdsa_verify_batch(
            ctx,
            msgs.as_ptr() as *const u8,
            pks.as_ptr() as *const u8,
            sigs.as_ptr() as *const u8,
            2,
            out.as_mut_ptr(),
        )
    };
    if rc != BLVM_SECP_GPU_OK {
        return Err("ECDSA init KAT launch failed");
    }
    if out[0] == 0 {
        return Err("ECDSA init KAT rejected known-good vector");
    }
    if out[1] != 0 {
        return Err("ECDSA init KAT accepted high-S");
    }

    // Schnorr: BIP340 vector 0 (accept) + r = p (reject).
    let sch_msg = [0u8; 32];
    let sch_pk: [u8; 32] = [
        0xf9, 0x30, 0x8a, 0x01, 0x92, 0x58, 0xc3, 0x10, 0x49, 0x34, 0x4f, 0x85, 0xf8, 0x9d, 0x52,
        0x29, 0xb5, 0x31, 0xc8, 0x45, 0x83, 0x6f, 0x99, 0xb0, 0x86, 0x01, 0xf1, 0x13, 0xbc, 0xe0,
        0x36, 0xf9,
    ];
    let sch_sig: [u8; 64] = [
        0xe9, 0x07, 0x83, 0x1f, 0x80, 0x84, 0x8d, 0x10, 0x69, 0xa5, 0x37, 0x1b, 0x40, 0x24, 0x10,
        0x36, 0x4b, 0xdf, 0x1c, 0x5f, 0x83, 0x07, 0xb0, 0x08, 0x4c, 0x55, 0xf1, 0xce, 0x2d, 0xca,
        0x82, 0x15, 0x25, 0xf6, 0x6a, 0x4a, 0x85, 0xea, 0x8b, 0x71, 0xe4, 0x82, 0xa7, 0x4f, 0x38,
        0x2d, 0x2c, 0xe5, 0xeb, 0xee, 0xe8, 0xfd, 0xb2, 0x17, 0x2f, 0x47, 0x7d, 0xf4, 0x90, 0x0d,
        0x31, 0x05, 0x36, 0xc0,
    ];
    let mut bad_r = sch_sig;
    bad_r[..32].copy_from_slice(&[
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
        0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xff, 0xff,
        0xfc, 0x2f,
    ]);
    let msgs2 = [sch_msg, sch_msg];
    let pks2 = [sch_pk, sch_pk];
    let sigs2 = [sch_sig, bad_r];
    let mut out2 = [0u8; 2];
    let rc = unsafe {
        blvm_secp_gpu_schnorr_verify_batch(
            ctx,
            msgs2.as_ptr() as *const u8,
            pks2.as_ptr() as *const u8,
            sigs2.as_ptr() as *const u8,
            2,
            out2.as_mut_ptr(),
        )
    };
    if rc != BLVM_SECP_GPU_OK {
        return Err("Schnorr init KAT launch failed");
    }
    if out2[0] == 0 {
        return Err("Schnorr init KAT rejected BIP340 vector0");
    }
    if out2[1] != 0 {
        return Err("Schnorr init KAT accepted r≥p");
    }
    Ok(())
}

fn try_create_pool() -> Option<GpuPool> {
    if env_disabled("BLVM_SECP_GPU") {
        eprintln!("[BLVM_SECP_GPU] disabled via BLVM_SECP_GPU — using CPU");
        return None;
    }
    unsafe {
        if blvm_secp_gpu_is_available() == 0 {
            eprintln!("[BLVM_SECP_GPU] CUDA unavailable — using CPU");
            return None;
        }
        let n_dev = blvm_secp_gpu_device_count();
        if n_dev == 0 {
            eprintln!("[BLVM_SECP_GPU] no CUDA devices — using CPU");
            return None;
        }
        let want = gpu_ctx_count();
        let mut slots = Vec::with_capacity(want);
        let mut created = 0u32;
        let mut first_ctx: *mut GpuCtx = std::ptr::null_mut();
        for i in 0..want {
            let dev = (i as u32) % n_dev;
            match try_create_one(dev) {
                Some(st) => {
                    if first_ctx.is_null() {
                        first_ctx = st.ctx;
                    }
                    created += 1;
                    slots.push(Mutex::new(Some(st)));
                }
                None => {
                    slots.push(Mutex::new(None));
                }
            }
        }
        if created == 0 {
            eprintln!("[BLVM_SECP_GPU] all CUDA ctx creates failed — using CPU");
            return None;
        }
        if let Err(reason) = device_init_kat(first_ctx) {
            eprintln!("[BLVM_SECP_GPU] init KAT failed ({reason}) — using CPU");
            // Drop destroys all ctxs.
            drop(slots);
            return None;
        }
        let min = gpu_batch_min();
        let submitters = gpu_submitters_enabled();
        let host_bridge = env_truthy("BLVM_SECP_GPU_HOST_BRIDGE");
        eprintln!(
            "[BLVM_SECP_GPU] native CUDA ready ctxs={created}/{want} devices={n_dev} \
             batch_min={min} submitters={} timers={} device_crypto={} init_kat=ok — \
             ECDSA/Schnorr batches ≥{min} use GPU (failover to CPU on launch failure)",
            if submitters { "on" } else { "off" },
            if gpu_timers_enabled() { "on" } else { "off" },
            if host_bridge { "host-bridge" } else { "on" }
        );
        Some(GpuPool {
            slots,
            rr: AtomicUsize::new(0),
        })
    }
}

fn gpu_pool() -> Option<&'static GpuPool> {
    static POOL: OnceLock<Option<GpuPool>> = OnceLock::new();
    POOL.get_or_init(try_create_pool).as_ref()
}

fn disable_slot(guard: &mut Option<GpuState>, reason: &str) {
    let dev = guard.as_ref().map(|s| s.device_index);
    eprintln!("[BLVM_SECP_GPU] {reason} device={dev:?} — sticky failover for this ctx slot");
    *guard = None;
}

/// Acquire a live GPU slot. Prefers free slots (try_lock round-robin); otherwise
/// blocks on the next slot unless `BLVM_SECP_GPU_TRY_LOCK=1`.
/// Note: sticky per-thread affinity dens (H11) was REVERT (−3.4% wall vs H6).
fn lock_gpu_slot() -> Option<(MutexGuard<'static, Option<GpuState>>, u64)> {
    let pool = gpu_pool()?;
    let n = pool.slots.len();
    if n == 0 {
        return None;
    }
    let t0 = Instant::now();
    let start = pool.rr.fetch_add(1, Ordering::Relaxed) % n;
    for i in 0..n {
        let idx = (start + i) % n;
        if let Ok(guard) = pool.slots[idx].try_lock() {
            if guard.is_some() {
                return Some((guard, t0.elapsed().as_nanos() as u64));
            }
        }
    }
    if gpu_try_lock_from_env() {
        return None;
    }
    for i in 0..n {
        let idx = (start + i) % n;
        if let Ok(guard) = pool.slots[idx].lock() {
            if guard.is_some() {
                return Some((guard, t0.elapsed().as_nanos() as u64));
            }
        }
    }
    None
}

/// After a successful CUDA launch, fill verdicts via CPU-only verify (no GPU re-entry).
/// Used when device crypto self-test failed at pool init (host-bridge mode).
fn ecdsa_host_bridge_fill(
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
    sigs: &[[u8; 64]],
    n: usize,
) -> Vec<bool> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(crate::ecdsa::ecdsa_verify_one_compact(
            &sigs[i], &msgs[i], &pubkeys[i],
        ));
    }
    out
}

fn schnorr_host_bridge_fill(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
) -> Vec<bool> {
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        out.push(crate::schnorr::schnorr_verify(
            &sigs[i], msgs[i], &pubkeys[i],
        ));
    }
    out
}

/// Trust device kernel verdicts unless `BLVM_SECP_GPU_HOST_BRIDGE=1`.
fn device_crypto_ok() -> bool {
    static OK: OnceLock<bool> = OnceLock::new();
    *OK.get_or_init(|| !env_truthy("BLVM_SECP_GPU_HOST_BRIDGE"))
}

fn run_ecdsa_on_ctx(
    ctx: *mut GpuCtx,
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
    sigs: &[[u8; 64]],
    n: usize,
) -> Result<Vec<bool>, c_int> {
    let mut out = vec![0u8; n];
    let rc = unsafe {
        blvm_secp_gpu_ecdsa_verify_batch(
            ctx,
            msgs.as_ptr() as *const u8,
            pubkeys.as_ptr() as *const u8,
            sigs.as_ptr() as *const u8,
            n,
            out.as_mut_ptr(),
        )
    };
    if rc != BLVM_SECP_GPU_OK {
        return Err(rc);
    }
    if device_crypto_ok() {
        Ok(out.into_iter().map(|b| b != 0).collect())
    } else {
        Ok(ecdsa_host_bridge_fill(msgs, pubkeys, sigs, n))
    }
}

/// Owned ECDSA batch job for async submitters / wave merge.
pub struct EcdsaGpuJob {
    pub msgs: Vec<[u8; 32]>,
    pub pubkeys: Vec<[u8; 33]>,
    pub sigs: Vec<[u8; 64]>,
    pub reply: std::sync::mpsc::SyncSender<Option<Vec<bool>>>,
}

struct SubmitterHub {
    tx: crossbeam_channel::Sender<EcdsaGpuJob>,
}

fn submitter_hub() -> Option<&'static SubmitterHub> {
    if !gpu_submitters_enabled() {
        return None;
    }
    static HUB: OnceLock<Option<SubmitterHub>> = OnceLock::new();
    HUB.get_or_init(|| {
        let pool = match gpu_pool() {
            Some(p) => p,
            None => return None,
        };
        let (tx, rx) = crossbeam_channel::unbounded::<EcdsaGpuJob>();
        let n_workers = pool.slots.len().max(1);
        for w in 0..n_workers {
            let rx = rx.clone();
            if std::thread::Builder::new()
                .name(format!("blvm-secp-gpu-{w}"))
                .spawn(move || {
                    while let Ok(job) = rx.recv() {
                        let n = job
                            .sigs
                            .len()
                            .min(job.msgs.len())
                            .min(job.pubkeys.len());
                        if n == 0 {
                            let _ = job.reply.send(Some(Vec::new()));
                            continue;
                        }
                        let t_lock = Instant::now();
                        let locked = lock_gpu_slot();
                        let lock_ns = t_lock.elapsed().as_nanos() as u64;
                        let Some((mut guard, _)) = locked else {
                            let _ = job.reply.send(None);
                            continue;
                        };
                        let Some(st) = guard.as_mut() else {
                            let _ = job.reply.send(None);
                            continue;
                        };
                        let t_k = Instant::now();
                        let result = run_ecdsa_on_ctx(
                            st.ctx,
                            &job.msgs[..n],
                            &job.pubkeys[..n],
                            &job.sigs[..n],
                            n,
                        );
                        let kernel_ns = t_k.elapsed().as_nanos() as u64;
                        match result {
                            Ok(v) => {
                                if gpu_timers_enabled() {
                                    GPU_TIMERS.record(n, lock_ns, 0, kernel_ns);
                                }
                                let _ = job.reply.send(Some(v));
                            }
                            Err(rc) => {
                                disable_slot(
                                    &mut guard,
                                    &format!("ECDSA batch launch failed rc={rc}"),
                                );
                                let _ = job.reply.send(None);
                            }
                        }
                    }
                })
                .is_err()
            {
                eprintln!("[BLVM_SECP_GPU] failed to spawn submitter {w}");
                return None;
            }
        }
        eprintln!("[BLVM_SECP_GPU] submitter threads={n_workers} started");
        Some(SubmitterHub { tx })
    })
    .as_ref()
}

/// Enqueue an owned ECDSA job to submitter threads. Returns `None` if submitters off / GPU down.
pub fn enqueue_ecdsa_job(
    msgs: Vec<[u8; 32]>,
    pubkeys: Vec<[u8; 33]>,
    sigs: Vec<[u8; 64]>,
) -> Option<std::sync::mpsc::Receiver<Option<Vec<bool>>>> {
    let hub = submitter_hub()?;
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n < gpu_batch_min() && n > 0 {
        return None;
    }
    let (rtx, rrx) = std::sync::mpsc::sync_channel(1);
    let job = EcdsaGpuJob {
        msgs,
        pubkeys,
        sigs,
        reply: rtx,
    };
    hub.tx.send(job).ok()?;
    Some(rrx)
}

fn ecdsa_verify_batch_ungated(
    sigs: &[[u8; 64]],
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
    n: usize,
) -> Option<Vec<bool>> {
    if gpu_submitters_enabled() {
        let rx = enqueue_ecdsa_job(msgs[..n].to_vec(), pubkeys[..n].to_vec(), sigs[..n].to_vec())?;
        return rx.recv().ok().flatten();
    }

    let timers = gpu_timers_enabled();
    let t_lock = Instant::now();
    let (mut guard, lock_ns_pass) = lock_gpu_slot()?;
    let lock_ns = if timers {
        t_lock.elapsed().as_nanos() as u64
    } else {
        lock_ns_pass
    };
    let ctx = guard.as_mut()?.ctx;

    let pack_ns = 0u64;
    let t_k = Instant::now();
    match run_ecdsa_on_ctx(ctx, &msgs[..n], &pubkeys[..n], &sigs[..n], n) {
        Ok(v) => {
            let kernel_ns = t_k.elapsed().as_nanos() as u64;
            if timers {
                GPU_TIMERS.record(n, lock_ns, pack_ns, kernel_ns);
            }
            Some(v)
        }
        Err(rc) => {
            disable_slot(&mut guard, &format!("ECDSA batch launch failed rc={rc}"));
            None
        }
    }
}

/// Per-item ECDSA verify on GPU. Returns `None` if GPU unavailable or launch failed.
pub fn try_ecdsa_verify_batch(
    sigs: &[[u8; 64]],
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
) -> Option<Vec<bool>> {
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n == 0 {
        return Some(Vec::new());
    }
    if n < gpu_batch_min() {
        return None;
    }
    ecdsa_verify_batch_ungated(sigs, msgs, pubkeys, n)
}

/// Like [`try_ecdsa_verify_batch`] but ignores `GPU_BATCH_MIN` (benches / microbenchmarks).
#[doc(hidden)]
pub fn try_ecdsa_verify_batch_ungated(
    sigs: &[[u8; 64]],
    msgs: &[[u8; 32]],
    pubkeys: &[[u8; 33]],
) -> Option<Vec<bool>> {
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n == 0 {
        return Some(Vec::new());
    }
    ecdsa_verify_batch_ungated(sigs, msgs, pubkeys, n)
}

fn schnorr_verify_batch_ungated(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
    n: usize,
) -> Option<Vec<bool>> {
    let (mut guard, lock_ns) = lock_gpu_slot()?;
    let ctx = guard.as_mut()?.ctx;

    let t_pack = Instant::now();
    let mut msg_buf = Vec::with_capacity(n * 32);
    for m in msgs.iter().take(n) {
        msg_buf.extend_from_slice(m);
    }
    let pack_ns = t_pack.elapsed().as_nanos() as u64;
    let mut out = vec![0u8; n];
    let t_k = Instant::now();
    let rc = unsafe {
        blvm_secp_gpu_schnorr_verify_batch(
            ctx,
            msg_buf.as_ptr(),
            pubkeys.as_ptr() as *const u8,
            sigs.as_ptr() as *const u8,
            n,
            out.as_mut_ptr(),
        )
    };
    let kernel_ns = t_k.elapsed().as_nanos() as u64;
    if rc != BLVM_SECP_GPU_OK {
        disable_slot(&mut guard, &format!("Schnorr batch launch failed rc={rc}"));
        return None;
    }
    if gpu_timers_enabled() {
        GPU_TIMERS.record(n, lock_ns, pack_ns, kernel_ns);
    }
    if device_crypto_ok() {
        Some(out.into_iter().map(|b| b != 0).collect())
    } else {
        Some(schnorr_host_bridge_fill(sigs, msgs, pubkeys, n))
    }
}

/// Per-item Schnorr verify on GPU. `msgs` must each be exactly 32 bytes (BIP-340).
pub fn try_schnorr_verify_batch(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
) -> Option<Vec<bool>> {
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n == 0 {
        return Some(Vec::new());
    }
    if n < gpu_batch_min() {
        return None;
    }
    if msgs.iter().take(n).any(|m| m.len() != 32) {
        return None;
    }
    schnorr_verify_batch_ungated(sigs, msgs, pubkeys, n)
}

/// Like [`try_schnorr_verify_batch`] but ignores `GPU_BATCH_MIN` (benches).
#[doc(hidden)]
pub fn try_schnorr_verify_batch_ungated(
    sigs: &[[u8; 64]],
    msgs: &[&[u8]],
    pubkeys: &[[u8; 32]],
) -> Option<Vec<bool>> {
    let n = sigs.len().min(msgs.len()).min(pubkeys.len());
    if n == 0 {
        return Some(Vec::new());
    }
    if msgs.iter().take(n).any(|m| m.len() != 32) {
        return None;
    }
    schnorr_verify_batch_ungated(sigs, msgs, pubkeys, n)
}

/// True if at least one CUDA GPU context is live (lazy).
pub fn gpu_available() -> bool {
    let pool = match gpu_pool() {
        Some(p) => p,
        None => return false,
    };
    for slot in &pool.slots {
        if let Ok(g) = slot.lock() {
            if g.is_some() {
                return true;
            }
        }
    }
    false
}
