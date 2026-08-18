/* BLVM-owned CUDA secp256k1 GPU C ABI (feature `gpu`).
 *
 * Compact SoA layouts only:
 *   ECDSA:  msg32 | pubkey33 | sig64(r||s)  → out_results[i] = 0/1
 *   Schnorr: msg32 | pubkey_x32 | sig64      → out_results[i] = 0/1
 *
 * Each context is single-threaded; serialize externally (see Rust gpu.rs pool).
 */
#ifndef BLVM_SECP_GPU_H
#define BLVM_SECP_GPU_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BLVM_SECP_GPU_OK                 0
#define BLVM_SECP_GPU_ERR_UNAVAILABLE  100
#define BLVM_SECP_GPU_ERR_DEVICE       101
#define BLVM_SECP_GPU_ERR_LAUNCH       102
#define BLVM_SECP_GPU_ERR_MEMORY       103
#define BLVM_SECP_GPU_ERR_BAD_INPUT    104

/** PRE_G / PRE_G_128 entry count (WINDOW_G=15 → 2^(15-2)). */
#define BLVM_SECP_GPU_PRE_G_ENTRIES    8192
/** Bytes per GeStorage (2 × FeStorage × 4 × u64). */
#define BLVM_SECP_GPU_GE_STORAGE_SIZE  64

typedef struct blvm_secp_gpu_ctx blvm_secp_gpu_ctx;

/** 1 if CUDA runtime sees at least one device. */
int blvm_secp_gpu_is_available(void);

/** Number of CUDA devices (0 if unavailable). */
uint32_t blvm_secp_gpu_device_count(void);

int blvm_secp_gpu_ctx_create(blvm_secp_gpu_ctx** ctx_out, uint32_t device_index);
void blvm_secp_gpu_ctx_destroy(blvm_secp_gpu_ctx* ctx);

/** Upload host PRE_G / PRE_G_128 tables (each n_entries × 64 bytes, n_entries == 8192).
 *  Layout must match Rust GeStorage / FeStorage LE u64 limbs (== CUDA Fe.d[4]).
 *  Call once after ctx_create; required for WINDOW_G device Strauss. */
int blvm_secp_gpu_ctx_set_pre_g(
    blvm_secp_gpu_ctx* ctx,
    const uint8_t* pre_g,
    const uint8_t* pre_g_128,
    size_t n_entries);

int blvm_secp_gpu_ecdsa_verify_batch(
    blvm_secp_gpu_ctx* ctx,
    const uint8_t* msg_hashes32,
    const uint8_t* pubkeys33,
    const uint8_t* sigs64,
    size_t count,
    uint8_t* out_results);

int blvm_secp_gpu_schnorr_verify_batch(
    blvm_secp_gpu_ctx* ctx,
    const uint8_t* msg_hashes32,
    const uint8_t* pubkeys_x32,
    const uint8_t* sigs64,
    size_t count,
    uint8_t* out_results);

const char* blvm_secp_gpu_error_str(int err);

#ifdef __cplusplus
}
#endif

#endif /* BLVM_SECP_GPU_H */
