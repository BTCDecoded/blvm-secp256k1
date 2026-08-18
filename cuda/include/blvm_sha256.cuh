/* Compact SHA-256 for BIP-340 (host+device). */
#pragma once
#include "blvm_hd.h"
#include <stddef.h>
#include <stdint.h>

BLVM_HD_INLINE uint32_t blvm_rotr32(uint32_t x, uint32_t n) {
    return (x >> n) | (x << (32 - n));
}

BLVM_HD void blvm_sha256_fixed(const uint8_t* msg, int len, uint8_t out[32]) {
    const uint32_t K[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u,
        0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu,
        0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu,
        0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau, 0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
        0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu,
        0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
        0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u, 0x19a4c116u,
        0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u,
        0xc67178f2u};
    uint32_t h0 = 0x6a09e667u, h1 = 0xbb67ae85u, h2 = 0x3c6ef372u, h3 = 0xa54ff53au;
    uint32_t h4 = 0x510e527fu, h5 = 0x9b05688cu, h6 = 0x1f83d9abu, h7 = 0x5be0cd19u;
    int blocks = (len + 9 + 63) / 64;
    for (int bi = 0; bi < blocks; ++bi) {
        uint8_t block[64];
        for (int i = 0; i < 64; ++i) block[i] = 0;
        int start = bi * 64;
        for (int i = 0; i < 64; ++i) {
            int idx = start + i;
            if (idx < len) block[i] = msg[idx];
            else if (idx == len) block[i] = 0x80;
        }
        if (bi == blocks - 1) {
            uint64_t bitlen = (uint64_t)len * 8ull;
            for (int i = 0; i < 8; ++i) block[63 - i] = (uint8_t)(bitlen >> (8 * i));
        }
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = ((uint32_t)block[4 * i] << 24) | ((uint32_t)block[4 * i + 1] << 16) |
                   ((uint32_t)block[4 * i + 2] << 8) | (uint32_t)block[4 * i + 3];
        for (int i = 16; i < 64; ++i) {
            uint32_t s0 = blvm_rotr32(w[i - 15], 7) ^ blvm_rotr32(w[i - 15], 18) ^ (w[i - 15] >> 3);
            uint32_t s1 = blvm_rotr32(w[i - 2], 17) ^ blvm_rotr32(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h0, b = h1, c = h2, d = h3, e = h4, f = h5, g = h6, hh = h7;
        for (int i = 0; i < 64; ++i) {
            uint32_t S1 = blvm_rotr32(e, 6) ^ blvm_rotr32(e, 11) ^ blvm_rotr32(e, 25);
            uint32_t ch = (e & f) ^ ((~e) & g);
            uint32_t t1 = hh + S1 + ch + K[i] + w[i];
            uint32_t S0 = blvm_rotr32(a, 2) ^ blvm_rotr32(a, 13) ^ blvm_rotr32(a, 22);
            uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t t2 = S0 + maj;
            hh = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        h0 += a;
        h1 += b;
        h2 += c;
        h3 += d;
        h4 += e;
        h5 += f;
        h6 += g;
        h7 += hh;
    }
    uint32_t hv[8] = {h0, h1, h2, h3, h4, h5, h6, h7};
    for (int i = 0; i < 8; ++i) {
        out[4 * i] = (uint8_t)(hv[i] >> 24);
        out[4 * i + 1] = (uint8_t)(hv[i] >> 16);
        out[4 * i + 2] = (uint8_t)(hv[i] >> 8);
        out[4 * i + 3] = (uint8_t)hv[i];
    }
}
