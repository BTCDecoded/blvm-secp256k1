#!/usr/bin/env bash
# Build and run bitcoin-core secp256k1 benchmarks (reference timings for this crate).
# See libsecp_bench_mapping.md for column/units and pairing notes.
#
# Usage:
#   LIBSECP_SRC=~/src/secp256k1 SECP256K1_BENCH_ITERS=1000 ./contrib/run_libsecp_parity.sh
#
set -euo pipefail

ROOT="${LIBSECP_SRC:-${HOME}/src/secp256k1}"
export SECP256K1_BENCH_ITERS="${SECP256K1_BENCH_ITERS:-2000}"

if [[ ! -d "${ROOT}/.git" && ! -f "${ROOT}/configure.ac" ]]; then
  echo "Expected libsecp256k1 sources at: ${ROOT}" >&2
  echo "Clone with: git clone https://github.com/bitcoin-core/secp256k1.git ${ROOT}" >&2
  exit 1
fi

cd "${ROOT}"

if [[ ! -f Makefile ]]; then
  ./autogen.sh
  ./configure --enable-benchmark
fi

make -j"$(nproc)" bench bench_ecmult

echo ""
echo "=== libsecp bench (µs/op min, avg, max) — API-level rows ==="
echo "    SECP256K1_BENCH_ITERS=${SECP256K1_BENCH_ITERS}"
# Selective runs avoid pulling in unrelated groups; ellswift_ecdh only appears if you pass ellswift.
./bench ecdsa_sign ecdsa_verify ec_keygen ecdh schnorrsig_sign schnorrsig_verify

echo ""
echo "=== libsecp bench_ecmult — primitive parity (longer; lower SECP256K1_BENCH_ITERS for a quick pass) ==="
./bench_ecmult

echo ""
echo "Done. Compare to your own blvm timing harness if you maintain one."
echo "See contrib/libsecp_bench_mapping.md for row semantics."
