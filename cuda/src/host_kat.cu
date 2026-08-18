/* Host KATs for BLVM CUDA secp math (compiled with nvcc, runs on CPU). */
#include "../include/blvm_verify.cuh"
#include <cstdio>
#include <cstring>

/* Null PRE_G stubs so device HD paths link; host_kat never uploads tables. */
__device__ const GeStorage* blvm_d_pre_g = nullptr;
__device__ const GeStorage* blvm_d_pre_g_128 = nullptr;
__device__ const GeStorage* blvm_pre_g_get(void) { return blvm_d_pre_g; }
__device__ const GeStorage* blvm_pre_g_128_get(void) { return blvm_d_pre_g_128; }

static int fail = 0;
#define CHECK(cond, msg)                         \
    do {                                         \
        if (!(cond)) {                           \
            std::printf("FAIL: %s\n", msg);       \
            fail = 1;                            \
        } else {                                 \
            std::printf("ok: %s\n", msg);         \
        }                                        \
    } while (0)

static void test_fe_basic() {
    Fe a, b, c, two, three, six;
    fe_set0(&a);
    a.d[0] = 2;
    fe_set0(&b);
    b.d[0] = 3;
    fe_mul(&c, &a, &b);
    fe_set0(&six);
    six.d[0] = 6;
    CHECK(fe_eq(&c, &six), "2*3=6");
    fe_add(&two, &a, &a);
    Fe four;
    fe_set0(&four);
    four.d[0] = 4;
    CHECK(fe_eq(&two, &four), "2+2=4");
    (void)three;
    Fe inv2, one;
    fe_inv(&inv2, &a);
    fe_mul(&c, &a, &inv2);
    fe_set0(&one);
    one.d[0] = 1;
    CHECK(fe_eq(&c, &one), "2*inv(2)=1");
}

static void test_g_decompress() {
    const uint8_t g33[33] = {
        0x02, 0x79, 0xBE, 0x66, 0x7E, 0xF9, 0xDC, 0xBB, 0xAC, 0x55, 0xA0, 0x62, 0x95, 0xCE, 0x87, 0x0B,
        0x07, 0x02, 0x9B, 0xFC, 0xDB, 0x2D, 0xCE, 0x28, 0xD9, 0x59, 0xF2, 0x81, 0x5B, 0x16, 0xF8, 0x17,
        0x98};
    Gej g;
    CHECK(ge_from_compressed(&g, g33), "decompress G");
    Gej g2;
    gej_set_g(&g2);
    Fe x1, x2;
    gej_get_x(&x1, &g);
    gej_get_x(&x2, &g2);
    CHECK(fe_eq(&x1, &x2), "G.x matches set_g");
}

/* RFC6979-style vector: secret=1, msg=zeros → verify via our sign-less known path.
 * Use pubkey = G compressed, and a valid sig produced offline... 
 * Instead: k*G with k=1 is G; craft trivial check that 1*G decompresses.
 */
static int hexbyte(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}
static void unhex(const char* hex, uint8_t* out, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        out[i] = (uint8_t)((hexbyte(hex[2 * i]) << 4) | hexbyte(hex[2 * i + 1]));
    }
}

static void test_ecdsa_self() {
    uint8_t msg[32] = {0};
    msg[31] = 1;
    const uint8_t pk0[33] = {
        0x02, 0x79, 0xBE, 0x66, 0x7E, 0xF9, 0xDC, 0xBB, 0xAC, 0x55, 0xA0, 0x62, 0x95, 0xCE, 0x87, 0x0B,
        0x07, 0x02, 0x9B, 0xFC, 0xDB, 0x2D, 0xCE, 0x28, 0xD9, 0x59, 0xF2, 0x81, 0x5B, 0x16, 0xF8, 0x17,
        0x98};
    uint8_t sig0[64] = {0};
    CHECK(blvm_ecdsa_verify_one(msg, pk0, sig0, nullptr, nullptr) == 0, "zero sig rejects");

    /* sk=1, msg=...07 — produced by blvm-secp256k1 RFC6979 sign */
    uint8_t msg2[32], pk2[33], sig2[64];
    unhex("0000000000000000000000000000000000000000000000000000000000000007", msg2, 32);
    unhex("0279be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", pk2, 33);
    unhex("cb6547dc10a26c2bb5c8bcacdd869aa1047b6f6a55a232003b53533828cba103"
          "1a95cda9a2cdd05bbd6115097458a741aaf0b03758fc2a5e4d9fe539c8621e6f",
          sig2, 64);
    CHECK(blvm_ecdsa_verify_one(msg2, pk2, sig2, nullptr, nullptr) == 1,
          "RFC6979 sk=1 vector verifies");
    sig2[40] ^= 1;
    CHECK(blvm_ecdsa_verify_one(msg2, pk2, sig2, nullptr, nullptr) == 0, "mutated sig rejects");
}

/* BIP-340 test vector 0 (verify-only). */
static void test_schnorr_bip340_v0() {
    uint8_t msg[32], pk[32], sig[64];
    unhex("0000000000000000000000000000000000000000000000000000000000000000", msg, 32);
    unhex("F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9", pk, 32);
    unhex("E907831F80848D1069A5371B402410364BDF1C5F8307B0084C55F1CE2DCA8215"
          "25F66A4A85EA8B71E482A74F382D2CE5EBEEE8FDB2172F477DF4900D310536C0",
          sig, 64);
    CHECK(blvm_schnorr_verify_one(msg, pk, sig) == 1, "BIP340 vector0 verifies");
    sig[40] ^= 1;
    CHECK(blvm_schnorr_verify_one(msg, pk, sig) == 0, "BIP340 mutated rejects");
}

/* Wire encodings with x or r equal to field prime p must reject (no reduce-mod-p). */
static void test_encoding_ge_p() {
    uint8_t p32[32];
    unhex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", p32, 32);
    Fe fe;
    CHECK(fe_set_b32_limit(&fe, p32) == 0, "fe_set_b32_limit rejects p");

    Gej g;
    CHECK(ge_lift_x(&g, p32) == 0, "lift_x rejects x=p");

    uint8_t pk33[33];
    pk33[0] = 0x02;
    std::memcpy(pk33 + 1, p32, 32);
    CHECK(ge_from_compressed(&g, pk33) == 0, "compressed pk x=p rejects");

    uint8_t msg[32] = {0};
    uint8_t pk[32];
    uint8_t sig[64];
    unhex("F9308A019258C31049344F85F89D5229B531C845836F99B08601F113BCE036F9", pk, 32);
    unhex("E907831F80848D1069A5371B402410364BDF1C5F8307B0084C55F1CE2DCA8215"
          "25F66A4A85EA8B71E482A74F382D2CE5EBEEE8FDB2172F477DF4900D310536C0",
          sig, 64);
    std::memcpy(sig, p32, 32); /* r = p */
    CHECK(blvm_schnorr_verify_one(msg, pk, sig) == 0, "BIP340 r=p rejects");
    std::memcpy(pk, p32, 32);
    unhex("E907831F80848D1069A5371B402410364BDF1C5F8307B0084C55F1CE2DCA8215"
          "25F66A4A85EA8B71E482A74F382D2CE5EBEEE8FDB2172F477DF4900D310536C0",
          sig, 64);
    CHECK(blvm_schnorr_verify_one(msg, pk, sig) == 0, "BIP340 pk x=p rejects");

    uint8_t ecdsa_pk[33];
    ecdsa_pk[0] = 0x02;
    std::memcpy(ecdsa_pk + 1, p32, 32);
    uint8_t ecdsa_sig[64];
    unhex("cb6547dc10a26c2bb5c8bcacdd869aa1047b6f6a55a232003b53533828cba103"
          "1a95cda9a2cdd05bbd6115097458a741aaf0b03758fc2a5e4d9fe539c8621e6f",
          ecdsa_sig, 64);
    uint8_t ecdsa_msg[32];
    unhex("0000000000000000000000000000000000000000000000000000000000000007", ecdsa_msg, 32);
    CHECK(blvm_ecdsa_verify_one(ecdsa_msg, ecdsa_pk, ecdsa_sig, nullptr, nullptr) == 0,
          "ECDSA pk x=p rejects");
}

int main() {
    test_fe_basic();
    test_g_decompress();
    test_ecdsa_self();
    test_schnorr_bip340_v0();
    test_encoding_ge_p();
    if (fail) {
        std::printf("HOST_KAT FAILED\n");
        return 1;
    }
    std::printf("HOST_KAT PASSED\n");
    return 0;
}
