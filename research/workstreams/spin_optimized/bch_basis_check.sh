#!/usr/bin/env bash
set -euo pipefail
root=${1:?root}
cd "$root"
result=workstreams/spin_optimized/measurements/bch-basis
mkdir -p "$result"
vendor=constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation/vendor
for mode in global local8 local16 local32 single reload; do
 build="build-bch-$mode"
 g++ -std=c++20 -O2 -march=x86-64-v3 -mpclmul -mvpclmulqdq \
  -I"$build/generated" -I"$vendor" workstreams/spin_optimized/bch_basis_test.cpp \
  "$build/libspin_half_transpose.a" -o "$build/spin_bch_basis_standalone"
 flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  "$build/spin_bch_basis_standalone" | tee "$result/$mode.log"
done
