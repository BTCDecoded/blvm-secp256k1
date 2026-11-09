//! Runtime CUDA probe. Compiled only with `cfg(has_cuda_kernels)`.
//!
//! The kernel archive is linked into the binary. `libcudart` is not. Host
//! stubs (see `build.rs`) `dlopen` the runtime before C++ constructors register
//! the fatbin. This module is the dispatch gate: one `OnceLock` probe, then
//! either the GPU entry points or the CPU path.
//!
//! `cuInit` and `cuDeviceGetCount` are driver-API symbols. On Linux they live
//! in `libcuda.so.1`, which `libcudart.so` loads. Windows exposes them from
//! `nvcuda.dll`. macOS exposes them from `libcuda.dylib`.

use std::sync::OnceLock;

use libloading::Library;

/// Loaded driver and runtime. Held so `dlclose` does not run while kernels do.
pub(crate) struct CudaHandle {
    _runtime: Library,
    _driver: Library,
}

/// `CUDA_ERROR_ALREADY_INITIALIZED`. The process constructor may have called it.
const CUDA_ERROR_ALREADY_INITIALIZED: i32 = 3;

/// True when the driver initialized and reported at least one device.
pub(crate) fn driver_ready(init_rc: i32, device_count: i32) -> bool {
    (init_rc == 0 || init_rc == CUDA_ERROR_ALREADY_INITIALIZED) && device_count > 0
}

pub(crate) fn handle() -> Option<&'static CudaHandle> {
    static SLOT: OnceLock<Option<CudaHandle>> = OnceLock::new();
    SLOT.get_or_init(probe).as_ref()
}

fn runtime_names() -> &'static [&'static str] {
    if cfg!(target_os = "windows") {
        &["cudart64_13.dll", "cudart64_12.dll", "cudart64_110.dll"]
    } else if cfg!(target_os = "macos") {
        &["libcudart.dylib"]
    } else {
        &[
            "libcudart.so",
            "libcudart.so.13",
            "libcudart.so.12",
            "libcudart.so.11.0",
        ]
    }
}

fn driver_names() -> &'static [&'static str] {
    if cfg!(target_os = "windows") {
        &["nvcuda.dll"]
    } else if cfg!(target_os = "macos") {
        &["libcuda.dylib"]
    } else {
        &["libcuda.so.1", "libcuda.so"]
    }
}

fn open_first(names: &[&str]) -> Option<Library> {
    for name in names {
        // SAFETY: the library is a CUDA runtime or driver. Its constructors
        // only register CUDA state. A missing file returns an error.
        if let Ok(lib) = unsafe { Library::new(name) } {
            return Some(lib);
        }
    }
    None
}

fn probe() -> Option<CudaHandle> {
    let runtime = open_first(runtime_names())?;
    let driver = open_first(driver_names())?;

    type CuInit = unsafe extern "C" fn(u32) -> i32;
    type CuDeviceGetCount = unsafe extern "C" fn(*mut i32) -> i32;

    // SAFETY: symbols are the CUDA driver ABI. `cuInit(0)` and
    // `cuDeviceGetCount` write a single i32. A missing symbol is `Err`.
    unsafe {
        let cu_init: libloading::Symbol<CuInit> = driver.get(b"cuInit\0").ok()?;
        let rc = cu_init(0);
        let cu_count: libloading::Symbol<CuDeviceGetCount> =
            driver.get(b"cuDeviceGetCount\0").ok()?;
        let mut n = 0i32;
        if cu_count(&mut n) != 0 {
            return None;
        }
        if !driver_ready(rc, n) {
            return None;
        }
    }

    Some(CudaHandle {
        _runtime: runtime,
        _driver: driver,
    })
}

#[cfg(test)]
mod tests {
    #[test]
    fn missing_library_is_none() {
        assert!(super::open_first(&["libblvm_no_such_cudart.so"]).is_none());
    }

    #[test]
    fn probe_twice_agrees() {
        let first = super::handle().is_some();
        let second = super::handle().is_some();
        assert_eq!(first, second);
    }

    #[test]
    fn driver_ready_requires_a_device() {
        assert!(super::driver_ready(0, 1));
        assert!(super::driver_ready(0, 2));
        assert!(super::driver_ready(
            super::CUDA_ERROR_ALREADY_INITIALIZED,
            1
        ));
        assert!(!super::driver_ready(0, 0));
        assert!(!super::driver_ready(0, -1));
        // `CUDA_ERROR_NO_DEVICE`.
        assert!(!super::driver_ready(100, 1));
        assert!(!super::driver_ready(100, 0));
    }
}
