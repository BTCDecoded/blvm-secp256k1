/* Host+device annotation for KATs and kernels. */
#pragma once
#ifdef __CUDACC__
#define BLVM_HD __host__ __device__
#define BLVM_HD_INLINE __host__ __device__ __forceinline__
#else
#define BLVM_HD
#define BLVM_HD_INLINE inline
#endif
