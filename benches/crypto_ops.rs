//! Performance baselines aligned with bitcoin-core **secp256k1** `bench` / `bench_ecmult`
//! workloads where a public analogue exists. See `contrib/libsecp_bench_mapping.md`.
//!
//! Run: `cargo bench --profile release --bench crypto_ops`

use std::cell::RefCell;
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};

use blvm_secp256k1::ecdsa::{
    ecdsa_sig_sign, ecdsa_sign_der_rfc6979, ge_from_compressed, ge_to_compressed,
    pubkey_from_secret, verify_ecdsa_direct,
};
use blvm_secp256k1::ecmult::{ecmult, ecmult_gen, ecmult_multi};
use blvm_secp256k1::group::{Ge, Gej, generator_g};
use blvm_secp256k1::scalar::Scalar;
use blvm_secp256k1::schnorr::{
    Keypair, schnorr_sign, schnorr_sign_with_keypair, schnorr_verify, xonly_pubkey_from_secret,
};
use blvm_secp256k1::{ecdh, ecmult_const, ecmult_gen_const};
use criterion::{Criterion, black_box, criterion_group, criterion_main};
use sha2::{Digest, Sha256};

/// Compressed pubkey and scalar from `src/modules/ecdh/bench_impl.h` (libsecp ECDH bench).
const LIBSECP_ECDH_COMPRESSED: [u8; 33] = [
    0x03, 0x54, 0x94, 0xc1, 0x5d, 0x32, 0x09, 0x97, 0x06, 0xc2, 0x39, 0x5f, 0x94, 0x34, 0x87, 0x45,
    0xfd, 0x75, 0x7c, 0xe3, 0x0e, 0x4e, 0x8c, 0x90, 0xfb, 0xa2, 0xba, 0xd1, 0x84, 0xf8, 0x83, 0xc6,
    0x9f,
];

fn libsecp_ecdh_scalar_bytes() -> [u8; 32] {
    let mut s = [0u8; 32];
    for (i, slot) in s.iter_mut().enumerate() {
        *slot = (i as u8) + 1;
    }
    s
}

struct SchnorrSlot {
    sk: [u8; 32],
    msg: [u8; 32],
    aux: [u8; 32],
    sig64: [u8; 64],
    pk_xonly: [u8; 32],
    keypair: Keypair,
}

/// Same per-index byte pattern as `modules/schnorrsig/bench_impl.h` (keypair + msg + sign_custom).
fn schnorr_slots() -> &'static [SchnorrSlot] {
    static SLOTS: OnceLock<Vec<SchnorrSlot>> = OnceLock::new();
    SLOTS.get_or_init(|| {
        const N: usize = 4096;
        let mut v = Vec::with_capacity(N);
        for i in 0..N {
            let mut sk = [0u8; 32];
            sk[0] = i as u8;
            sk[1] = (i >> 8) as u8;
            sk[2] = (i >> 16) as u8;
            sk[3] = (i >> 24) as u8;
            sk[4..].fill(b's');
            let mut msg = [0u8; 32];
            msg[..4].copy_from_slice(&sk[..4]);
            msg[4..].fill(b'm');
            let aux = [0u8; 32];
            let pk_xonly = xonly_pubkey_from_secret(&sk).expect("slot pk");
            let sig64 = schnorr_sign(&sk, &msg, &aux).expect("slot sig");
            let keypair = Keypair::from_seckey(&sk).expect("slot keypair");
            v.push(SchnorrSlot {
                sk,
                msg,
                aux,
                sig64,
                pk_xonly,
                keypair,
            });
        }
        v
    })
}

/// Secret-key chain from `bench.c` `bench_keygen_run`: start bytes `i+65`, then x-only bytes of each
/// compressed pubkey (same as `memcpy(data->key, pub33 + 1, 32)`).
fn keygen_secret_chain() -> &'static [[u8; 32]] {
    static CHAIN: OnceLock<Vec<[u8; 32]>> = OnceLock::new();
    CHAIN.get_or_init(|| {
        const LEN: usize = 65_536;
        let mut k = [0u8; 32];
        for (i, slot) in k.iter_mut().enumerate() {
            *slot = (i as u8) + 65;
        }
        let mut v = Vec::with_capacity(LEN);
        for _ in 0..LEN {
            v.push(k);
            let mut s = Scalar::zero();
            if s.set_b32(&k) || s.is_zero() {
                break;
            }
            let ge = pubkey_from_secret(&s);
            let c = ge_to_compressed(&ge);
            k.copy_from_slice(&c[1..33]);
        }
        v
    })
}

static SCHNORR_SIGN_IDX: AtomicUsize = AtomicUsize::new(0);
static SCHNORR_VERIFY_IDX: AtomicUsize = AtomicUsize::new(0);
static KEYGEN_IDX: AtomicUsize = AtomicUsize::new(0);
static ECMULT_HASH_IDX: AtomicUsize = AtomicUsize::new(0);
static ECMULT_CONST_PAIR_IDX: AtomicUsize = AtomicUsize::new(0);

fn bench_schnorr_sign_libsecp_pattern(c: &mut Criterion) {
    let slots = schnorr_slots();
    // Mirrors libsecp `bench_schnorrsig_sign` which uses a `secp256k1_keypair` with cached
    // pubkey, so the per-iteration work skips `d*G` and only does `k*G` for the nonce.
    c.bench_function("schnorr_sign_libsecp_pattern", |b| {
        b.iter(|| {
            let i = SCHNORR_SIGN_IDX.fetch_add(1, Ordering::Relaxed) % slots.len();
            let s = &slots[i];
            black_box(schnorr_sign_with_keypair(
                black_box(&s.keypair),
                black_box(&s.msg),
                black_box(&s.aux),
            ))
            .unwrap()
        })
    });
}

/// Backward-compatible bench: `schnorr_sign(sk, msg, aux)` — derives pubkey internally with
/// a CT `d*G`, so the per-iteration cost is **2** `ecmult_gen_const` calls instead of 1.
/// Use this only to track the legacy entry point; pair the libsecp comparison row with
/// `schnorr_sign_libsecp_pattern` (cached-pubkey variant).
fn bench_schnorr_sign_no_keypair(c: &mut Criterion) {
    let slots = schnorr_slots();
    c.bench_function("schnorr_sign_no_keypair", |b| {
        b.iter(|| {
            let i = SCHNORR_SIGN_IDX.fetch_add(1, Ordering::Relaxed) % slots.len();
            let s = &slots[i];
            black_box(schnorr_sign(
                black_box(&s.sk),
                black_box(&s.msg),
                black_box(&s.aux),
            ))
            .unwrap()
        })
    });
}

fn bench_schnorr_verify_libsecp_pattern(c: &mut Criterion) {
    let slots = schnorr_slots();
    c.bench_function("schnorr_verify_libsecp_pattern", |b| {
        b.iter(|| {
            let i = SCHNORR_VERIFY_IDX.fetch_add(1, Ordering::Relaxed) % slots.len();
            let s = &slots[i];
            black_box(schnorr_verify(
                black_box(&s.sig64),
                black_box(&s.msg),
                black_box(&s.pk_xonly),
            ))
        })
    });
}

/// Mirrors libsecp `bench_verify_run`: per iteration the last 3 bytes of the DER signature
/// are XORed with the iteration counter, so most iterations exercise the failure path
/// through the full verify pipeline. We accept either Some(false) (parsed but invalid) or
/// None (parse rejection from a corrupted DER) — both represent the same upstream cost
/// shape (DER parse + ecmult). The first iteration with `i==0` always verifies cleanly.
fn bench_ecdsa_verify_libsecp_pattern(c: &mut Criterion) {
    let seckey: [u8; 32] = {
        let mut k = [0u8; 32];
        for (i, slot) in k.iter_mut().enumerate() {
            *slot = (i as u8) + 65;
        }
        k
    };
    let msg: [u8; 32] = {
        let mut m = [0u8; 32];
        for (i, slot) in m.iter_mut().enumerate() {
            *slot = (i as u8) + 1;
        }
        m
    };
    let mut s = Scalar::zero();
    let _ = s.set_b32(&seckey);
    let pk_ge = pubkey_from_secret(&s);
    let pk33 = ge_to_compressed(&pk_ge);
    let der = ecdsa_sign_der_rfc6979(&msg, &seckey).expect("seed sig");
    let der_len = der.len();
    let mut sig_buf = vec![0u8; der_len];
    sig_buf.copy_from_slice(&der);

    let counter = AtomicUsize::new(0);
    c.bench_function("ecdsa_verify_libsecp_pattern", |b| {
        b.iter(|| {
            let i = counter.fetch_add(1, Ordering::Relaxed);
            // Restore baseline DER, then mangle last 3 bytes with iteration counter (libsecp
            // pattern: data->sig[len-1] ^= i; [-2] ^= i>>8; [-3] ^= i>>16). Most iterations
            // fail verify which is the realistic workload; a clean iteration happens when
            // the XOR mask is 0 (i.e. i==0).
            sig_buf.copy_from_slice(&der);
            sig_buf[der_len - 1] ^= (i & 0xff) as u8;
            sig_buf[der_len - 2] ^= ((i >> 8) & 0xff) as u8;
            sig_buf[der_len - 3] ^= ((i >> 16) & 0xff) as u8;
            black_box(verify_ecdsa_direct(
                black_box(&sig_buf),
                black_box(&pk33),
                black_box(&msg),
                false,
                false,
            ))
        })
    });
}

fn bench_ecdsa_sign_explicit_nonce(c: &mut Criterion) {
    let mut sec = Scalar::zero();
    sec.set_int(0x1234_5678);
    let mut msg = Scalar::zero();
    msg.set_int(0xdead_beef);
    let mut nonce = Scalar::zero();
    nonce.set_int(0x1111_2222);
    c.bench_function("ecdsa_sig_sign_explicit_nonce", |b| {
        b.iter(|| {
            black_box(ecdsa_sig_sign(
                black_box(&sec),
                black_box(&msg),
                black_box(&nonce),
            ))
            .unwrap()
        })
    });
}

/// Matches `bench.c` `bench_sign_run`: RFC6979 + DER each iteration, then `msg[0..32]` and
/// `key[0..32]` are replaced by `sig_buf[0..32]` and `sig_buf[32..64]` (same as libsecp’s
/// `unsigned char sig[74]` scratch).
fn bench_ecdsa_sign_der_rfc6979_libsecp_bench_loop(c: &mut Criterion) {
    let msg = RefCell::new({
        let mut m = [0u8; 32];
        for (i, slot) in m.iter_mut().enumerate() {
            *slot = i as u8 + 1;
        }
        m
    });
    let key = RefCell::new({
        let mut k = [0u8; 32];
        for (i, slot) in k.iter_mut().enumerate() {
            *slot = i as u8 + 65;
        }
        k
    });
    let mut sig_scratch = [0u8; 74];
    c.bench_function("ecdsa_sign_der_rfc6979_libsecp_bench_loop", |b| {
        b.iter(|| {
            let m = *msg.borrow();
            let k = *key.borrow();
            let der = ecdsa_sign_der_rfc6979(black_box(&m), black_box(&k)).expect("ecdsa sign");
            let n = der.len().min(74);
            sig_scratch[..n].copy_from_slice(&der[..n]);
            let mut mg = msg.borrow_mut();
            let mut kg = key.borrow_mut();
            for j in 0..32 {
                mg[j] = sig_scratch[j];
                kg[j] = sig_scratch[j + 32];
            }
            black_box(n)
        })
    });
}

fn bench_pubkey_from_secret_chained(c: &mut Criterion) {
    let chain = keygen_secret_chain();
    c.bench_function("pubkey_from_secret_chained", |b| {
        b.iter(|| {
            let i = KEYGEN_IDX.fetch_add(1, Ordering::Relaxed) % chain.len();
            let mut s = Scalar::zero();
            let _ = s.set_b32(&chain[i]);
            black_box(pubkey_from_secret(black_box(&s)))
        })
    });
}

fn bench_ecdh_libsecp_point(c: &mut Criterion) {
    let pk = ge_from_compressed(&LIBSECP_ECDH_COMPRESSED).expect("libsecp ecdh bench point");
    let mut sec = Scalar::zero();
    let _ = sec.set_b32(&libsecp_ecdh_scalar_bytes());
    c.bench_function("ecdh_libsecp_point", |b| {
        b.iter(|| black_box(ecdh::ecdh(black_box(&pk), black_box(&sec)).unwrap()))
    });
}

fn bench_ecmult_gen_const(c: &mut Criterion) {
    let mut k = Scalar::zero();
    k.set_int(0xabc);
    let mut out = Gej::default();
    c.bench_function("ecmult_gen_const", |b| {
        b.iter(|| {
            ecmult_gen_const(black_box(&mut out), black_box(&k));
        })
    });
}

/// Pseudorandom scalar per iteration (SHA256(counter)) — closer to libsecp `bench_ecmult`’s
/// large precomputed scalar table than a single fixed `k`.
fn bench_ecmult_gen_const_hashed_scalars(c: &mut Criterion) {
    let mut out = Gej::default();
    c.bench_function("ecmult_gen_const_hashed_scalars", |b| {
        b.iter(|| {
            let i = ECMULT_HASH_IDX.fetch_add(1, Ordering::Relaxed);
            let h: [u8; 32] = Sha256::digest(i.to_le_bytes()).into();
            let mut k = Scalar::zero();
            let _ = k.set_b32(&h);
            ecmult_gen_const(black_box(&mut out), black_box(&k));
        })
    });
}

fn bench_ecmult_const_affine(c: &mut Criterion) {
    let g = generator_g();
    let mut k = Scalar::zero();
    k.set_int(0xdef);
    let mut out = Gej::default();
    c.bench_function("ecmult_const_affine", |b| {
        b.iter(|| {
            ecmult_const(black_box(&mut out), black_box(&g), black_box(&k));
        })
    });
}

fn bench_ecmult_const_affine_hashed_scalars(c: &mut Criterion) {
    let g = generator_g();
    let mut out = Gej::default();
    c.bench_function("ecmult_const_affine_hashed_scalars", |b| {
        b.iter(|| {
            let i = ECMULT_HASH_IDX.fetch_add(1, Ordering::Relaxed);
            let h: [u8; 32] = Sha256::digest((i ^ 0x9e37_79b9).to_le_bytes()).into();
            let mut k = Scalar::zero();
            let _ = k.set_b32(&h);
            ecmult_const(black_box(&mut out), black_box(&g), black_box(&k));
        })
    });
}

/// Independent rotating affine base and scalar (4096 each), mimicking libsecp `bench_ecmult_const`
/// which indexes **`pubkeys[(offset1+i)%N]`** and **`scalars[(offset2+i)%N]`** separately.
fn ecmult_const_table_pair() -> &'static (Vec<Ge>, Vec<Scalar>) {
    static T: OnceLock<(Vec<Ge>, Vec<Scalar>)> = OnceLock::new();
    T.get_or_init(|| {
        const N: usize = 4096;
        let mut bases = Vec::with_capacity(N);
        let mut scalars = Vec::with_capacity(N);
        for i in 0..N {
            let mut sg = Scalar::zero();
            let _ = sg.set_b32(&Sha256::digest((i as u64).to_le_bytes()).into());
            let mut r = Gej::default();
            ecmult_gen(&mut r, &sg);
            let mut ge = Ge::default();
            ge.set_gej_var(&r);
            bases.push(ge);
            let mut sc = Scalar::zero();
            let _ = sc.set_b32(&Sha256::digest((i as u64 ^ 0x51f4_a552).to_le_bytes()).into());
            scalars.push(sc);
        }
        (bases, scalars)
    })
}

fn bench_ecmult_const_varying_affine_hashed_k(c: &mut Criterion) {
    let (bases, scalars) = ecmult_const_table_pair();
    let n = bases.len();
    let mut out = Gej::default();
    c.bench_function("ecmult_const_varying_affine_hashed_k", |b| {
        b.iter(|| {
            let i = ECMULT_CONST_PAIR_IDX.fetch_add(1, Ordering::Relaxed);
            let i1 = i % n;
            let i2 = i.wrapping_mul(0x9e37_79b9) % n;
            ecmult_const(
                black_box(&mut out),
                black_box(&bases[i1]),
                black_box(&scalars[i2]),
            );
        })
    });
}

/// Variable-time table-based `k*G` (blvm fast path). Compare to libsecp **`ecmult_gen`** only
/// qualitatively — libsecp’s public `ecmult_gen` in `bench_ecmult` is the signing/CT-style table.
fn bench_ecmult_gen_wnaf_table(c: &mut Criterion) {
    let mut k = Scalar::zero();
    k.set_int(0x135);
    let mut out = Gej::default();
    c.bench_function("ecmult_gen_wnaf_table", |b| {
        b.iter(|| {
            ecmult_gen(black_box(&mut out), black_box(&k));
        })
    });
}

/// `k*G` via `ecmult(..., R, inf, 0, Some(k))` — matches libsecp **`ecmult_0p_g`** (verify-style Strauss).
fn bench_ecmult_g_via_ecmult(c: &mut Criterion) {
    let mut k = Scalar::zero();
    k.set_int(0x135);
    let mut inf = Gej::default();
    inf.set_infinity();
    let zero = Scalar::zero();
    let mut out = Gej::default();
    c.bench_function("ecmult_g_via_ecmult", |b| {
        b.iter(|| {
            ecmult(
                black_box(&mut out),
                black_box(&inf),
                black_box(&zero),
                Some(black_box(&k)),
            );
        })
    });
}

/// `R = g_scalar * G + sum_i scalars[i] * points[i]` with **7** affine points — matches libsecp
/// **`ecmult_multi_7p_g`** (seven non-G terms plus one G coefficient in their callback layout).
fn bench_ecmult_multi_7p_g(c: &mut Criterion) {
    let mut g_scalar = Scalar::zero();
    g_scalar.set_int(101);
    let mut points = [Ge::default(); 7];
    let mut scalars = [Scalar::zero(); 7];
    for i in 0..7 {
        let mut sg = Scalar::zero();
        sg.set_int((i + 3) as u32);
        let mut r = Gej::default();
        ecmult_gen(&mut r, &sg);
        points[i].set_gej_var(&r);
        scalars[i].set_int((i * 7 + 11) as u32);
    }
    let mut acc = Gej::default();
    c.bench_function("ecmult_multi_7p_g", |b| {
        b.iter(|| {
            ecmult_multi(
                black_box(&mut acc),
                black_box(&g_scalar),
                black_box(&scalars),
                black_box(&points),
            );
        })
    });
}

/// `na * A + ng * G` with **A ≠ G** (Jacobian), matching libsecp **`ecmult_1p_g`** shape.
fn bench_ecmult_strauss_1p_g(c: &mut Criterion) {
    let mut base_s = Scalar::zero();
    base_s.set_int(19);
    let mut aj = Gej::default();
    ecmult_gen(&mut aj, &base_s);
    let mut na = Scalar::zero();
    na.set_int(0x9999);
    let mut ng = Scalar::zero();
    ng.set_int(0xaaaau32);
    let mut out = Gej::default();
    c.bench_function("ecmult_strauss_1p_g", |b| {
        b.iter(|| {
            ecmult(
                black_box(&mut out),
                black_box(&aj),
                black_box(&na),
                Some(black_box(&ng)),
            );
        })
    });
}

criterion_group!(
    benches,
    bench_schnorr_sign_libsecp_pattern,
    bench_schnorr_sign_no_keypair,
    bench_schnorr_verify_libsecp_pattern,
    bench_ecdsa_verify_libsecp_pattern,
    bench_ecdsa_sign_explicit_nonce,
    bench_ecdsa_sign_der_rfc6979_libsecp_bench_loop,
    bench_pubkey_from_secret_chained,
    bench_ecdh_libsecp_point,
    bench_ecmult_gen_const,
    bench_ecmult_gen_const_hashed_scalars,
    bench_ecmult_const_affine,
    bench_ecmult_const_affine_hashed_scalars,
    bench_ecmult_const_varying_affine_hashed_k,
    bench_ecmult_gen_wnaf_table,
    bench_ecmult_g_via_ecmult,
    bench_ecmult_multi_7p_g,
    bench_ecmult_strauss_1p_g,
);
criterion_main!(benches);
