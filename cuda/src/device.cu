/* BLVM secp GPU: context, grow-only buffers, batch launches. */
#include "../include/blvm_secp_gpu.h"
#include "../include/blvm_verify.cuh"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/* Device-visible PRE_G pointers (per CUDA device; synced before each launch). */
__device__ const GeStorage* blvm_d_pre_g = nullptr;
__device__ const GeStorage* blvm_d_pre_g_128 = nullptr;

__device__ const GeStorage* blvm_pre_g_get(void) { return blvm_d_pre_g; }
__device__ const GeStorage* blvm_pre_g_128_get(void) { return blvm_d_pre_g_128; }

struct blvm_secp_gpu_ctx {
    int device;
    cudaStream_t stream;
    /* Device buffers (grow-only). */
    uint8_t* d_msgs;
    uint8_t* d_pks;
    uint8_t* d_sigs;
    uint8_t* d_out;
    size_t cap_msgs;
    size_t cap_pks;
    size_t cap_sigs;
    size_t cap_out;
    /* PRE_G / PRE_G_128 (WINDOW_G=15); optional until set_pre_g. */
    GeStorage* d_pre_g;
    GeStorage* d_pre_g_128;
    size_t pre_g_entries;
};

static thread_local int tls_cuda_device = -1;

static int ensure_device(int device) {
    if (tls_cuda_device == device) return BLVM_SECP_GPU_OK;
    if (cudaSetDevice(device) != cudaSuccess) return BLVM_SECP_GPU_ERR_DEVICE;
    tls_cuda_device = device;
    return BLVM_SECP_GPU_OK;
}

static int ensure_dev_buf(uint8_t** p, size_t* cap, size_t need) {
    if (*cap >= need) return BLVM_SECP_GPU_OK;
    if (*p) cudaFree(*p);
    *p = nullptr;
    *cap = 0;
    if (cudaMalloc(p, need) != cudaSuccess) return BLVM_SECP_GPU_ERR_MEMORY;
    *cap = need;
    return BLVM_SECP_GPU_OK;
}

extern "C" int blvm_secp_gpu_is_available(void) {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) return 0;
    return n > 0 ? 1 : 0;
}

extern "C" uint32_t blvm_secp_gpu_device_count(void) {
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess) return 0;
    return n > 0 ? (uint32_t)n : 0;
}

extern "C" int blvm_secp_gpu_ctx_create(blvm_secp_gpu_ctx** ctx_out, uint32_t device_index) {
    if (!ctx_out) return BLVM_SECP_GPU_ERR_BAD_INPUT;
    *ctx_out = nullptr;
    int n = 0;
    if (cudaGetDeviceCount(&n) != cudaSuccess || n <= 0) return BLVM_SECP_GPU_ERR_UNAVAILABLE;
    if ((int)device_index >= n) return BLVM_SECP_GPU_ERR_DEVICE;
    if (ensure_device((int)device_index) != BLVM_SECP_GPU_OK) return BLVM_SECP_GPU_ERR_DEVICE;
    auto* ctx = new blvm_secp_gpu_ctx();
    std::memset(ctx, 0, sizeof(*ctx));
    ctx->device = (int)device_index;
    if (cudaStreamCreateWithFlags(&ctx->stream, cudaStreamNonBlocking) != cudaSuccess) {
        delete ctx;
        return BLVM_SECP_GPU_ERR_DEVICE;
    }
    *ctx_out = ctx;
    return BLVM_SECP_GPU_OK;
}

extern "C" int blvm_secp_gpu_ctx_set_pre_g(blvm_secp_gpu_ctx* ctx, const uint8_t* pre_g,
                                           const uint8_t* pre_g_128, size_t n_entries) {
    if (!ctx || !pre_g || !pre_g_128) return BLVM_SECP_GPU_ERR_BAD_INPUT;
    if (n_entries != BLVM_SECP_GPU_PRE_G_ENTRIES) return BLVM_SECP_GPU_ERR_BAD_INPUT;
    if (ensure_device(ctx->device) != BLVM_SECP_GPU_OK) return BLVM_SECP_GPU_ERR_DEVICE;

    const size_t bytes = n_entries * BLVM_SECP_GPU_GE_STORAGE_SIZE;
    GeStorage* d_g = nullptr;
    GeStorage* d_g128 = nullptr;
    if (cudaMalloc(&d_g, bytes) != cudaSuccess) return BLVM_SECP_GPU_ERR_MEMORY;
    if (cudaMalloc(&d_g128, bytes) != cudaSuccess) {
        cudaFree(d_g);
        return BLVM_SECP_GPU_ERR_MEMORY;
    }
    if (cudaMemcpy(d_g, pre_g, bytes, cudaMemcpyHostToDevice) != cudaSuccess ||
        cudaMemcpy(d_g128, pre_g_128, bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
        cudaFree(d_g);
        cudaFree(d_g128);
        return BLVM_SECP_GPU_ERR_MEMORY;
    }

    if (ctx->d_pre_g) cudaFree(ctx->d_pre_g);
    if (ctx->d_pre_g_128) cudaFree(ctx->d_pre_g_128);
    ctx->d_pre_g = d_g;
    ctx->d_pre_g_128 = d_g128;
    ctx->pre_g_entries = n_entries;

    /* Bind symbols now so a subsequent KAT without launch_batch still sees tables. */
    if (cudaMemcpyToSymbol(blvm_d_pre_g, &d_g, sizeof(d_g)) != cudaSuccess ||
        cudaMemcpyToSymbol(blvm_d_pre_g_128, &d_g128, sizeof(d_g128)) != cudaSuccess)
        return BLVM_SECP_GPU_ERR_MEMORY;
    return BLVM_SECP_GPU_OK;
}

extern "C" void blvm_secp_gpu_ctx_destroy(blvm_secp_gpu_ctx* ctx) {
    if (!ctx) return;
    ensure_device(ctx->device);
    if (ctx->d_msgs) cudaFree(ctx->d_msgs);
    if (ctx->d_pks) cudaFree(ctx->d_pks);
    if (ctx->d_sigs) cudaFree(ctx->d_sigs);
    if (ctx->d_out) cudaFree(ctx->d_out);
    if (ctx->d_pre_g) cudaFree(ctx->d_pre_g);
    if (ctx->d_pre_g_128) cudaFree(ctx->d_pre_g_128);
    if (ctx->stream) cudaStreamDestroy(ctx->stream);
    delete ctx;
}

static int launch_batch(blvm_secp_gpu_ctx* ctx, const uint8_t* msg_hashes32, size_t msg_stride,
                        const uint8_t* pubkeys, size_t pk_stride, const uint8_t* sigs64,
                        size_t count, uint8_t* out_results, bool ecdsa) {
    if (ensure_device(ctx->device) != BLVM_SECP_GPU_OK) return BLVM_SECP_GPU_ERR_DEVICE;

    int rc;
    if ((rc = ensure_dev_buf(&ctx->d_msgs, &ctx->cap_msgs, count * 32))) return rc;
    if ((rc = ensure_dev_buf(&ctx->d_pks, &ctx->cap_pks, count * pk_stride))) return rc;
    if ((rc = ensure_dev_buf(&ctx->d_sigs, &ctx->cap_sigs, count * 64))) return rc;
    if ((rc = ensure_dev_buf(&ctx->d_out, &ctx->cap_out, count))) return rc;
    /* PRE_G passed as kernel args (H12) — no per-launch MemcpyToSymbol. */

    /* Direct pageable H2D — pinned staging added host memcpy and slowed 1k–4k batches. */
    if (cudaMemcpyAsync(ctx->d_msgs, msg_hashes32, count * 32, cudaMemcpyHostToDevice, ctx->stream) !=
        cudaSuccess)
        return BLVM_SECP_GPU_ERR_MEMORY;
    if (cudaMemcpyAsync(ctx->d_pks, pubkeys, count * pk_stride, cudaMemcpyHostToDevice,
                        ctx->stream) != cudaSuccess)
        return BLVM_SECP_GPU_ERR_MEMORY;
    if (cudaMemcpyAsync(ctx->d_sigs, sigs64, count * 64, cudaMemcpyHostToDevice, ctx->stream) !=
        cudaSuccess)
        return BLVM_SECP_GPU_ERR_MEMORY;

    /* 128 threads: heavy register EC kernels lose occupancy at 256. */
    const int threads = 128;
    const int blocks = (int)((count + threads - 1) / threads);
    if (ecdsa) {
        blvm_ecdsa_verify_batch_kernel<<<blocks, threads, 0, ctx->stream>>>(
            ctx->d_msgs, ctx->d_pks, ctx->d_sigs, ctx->d_out, count, ctx->d_pre_g,
            ctx->d_pre_g_128);
    } else {
        blvm_schnorr_verify_batch_kernel<<<blocks, threads, 0, ctx->stream>>>(
            ctx->d_msgs, ctx->d_pks, ctx->d_sigs, ctx->d_out, count);
    }
    if (cudaGetLastError() != cudaSuccess) return BLVM_SECP_GPU_ERR_LAUNCH;

    if (cudaMemcpyAsync(out_results, ctx->d_out, count, cudaMemcpyDeviceToHost, ctx->stream) !=
        cudaSuccess)
        return BLVM_SECP_GPU_ERR_MEMORY;
    if (cudaStreamSynchronize(ctx->stream) != cudaSuccess) return BLVM_SECP_GPU_ERR_LAUNCH;
    (void)msg_stride;
    return BLVM_SECP_GPU_OK;
}

extern "C" int blvm_secp_gpu_ecdsa_verify_batch(blvm_secp_gpu_ctx* ctx, const uint8_t* msg_hashes32,
                                                const uint8_t* pubkeys33, const uint8_t* sigs64,
                                                size_t count, uint8_t* out_results) {
    if (!ctx || (!msg_hashes32 && count) || (!pubkeys33 && count) || (!sigs64 && count) ||
        !out_results)
        return BLVM_SECP_GPU_ERR_BAD_INPUT;
    if (count == 0) return BLVM_SECP_GPU_OK;
    return launch_batch(ctx, msg_hashes32, 32, pubkeys33, 33, sigs64, count, out_results, true);
}

extern "C" int blvm_secp_gpu_schnorr_verify_batch(blvm_secp_gpu_ctx* ctx, const uint8_t* msg_hashes32,
                                                  const uint8_t* pubkeys_x32, const uint8_t* sigs64,
                                                  size_t count, uint8_t* out_results) {
    if (!ctx || (!msg_hashes32 && count) || (!pubkeys_x32 && count) || (!sigs64 && count) ||
        !out_results)
        return BLVM_SECP_GPU_ERR_BAD_INPUT;
    if (count == 0) return BLVM_SECP_GPU_OK;
    return launch_batch(ctx, msg_hashes32, 32, pubkeys_x32, 32, sigs64, count, out_results, false);
}

extern "C" const char* blvm_secp_gpu_error_str(int err) {
    switch (err) {
        case BLVM_SECP_GPU_OK:
            return "ok";
        case BLVM_SECP_GPU_ERR_UNAVAILABLE:
            return "cuda unavailable";
        case BLVM_SECP_GPU_ERR_DEVICE:
            return "cuda device error";
        case BLVM_SECP_GPU_ERR_LAUNCH:
            return "cuda launch error";
        case BLVM_SECP_GPU_ERR_MEMORY:
            return "cuda memory error";
        case BLVM_SECP_GPU_ERR_BAD_INPUT:
            return "bad input";
        default:
            return "unknown error";
    }
}
