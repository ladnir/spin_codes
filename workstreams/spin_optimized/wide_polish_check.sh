#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/wide_polish
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-wide-polish \
  -C workstreams/spin_optimized/cmake/ForwardRecommended.cmake \
  -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
# The polished generator must reproduce the measured circuit byte-for-byte.
cmp build-wide-polish/bidirectional/wide/generated/WideCircuit.h \
  build-wide-schedule-dfs/bidirectional/wide/generated/WideCircuit.h
cmake --build build-wide-polish --target spin_wide_test spin_bidirectional_test \
  spin_k16_bidirectional_test spin_wide_benchmark -j2 > "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-wide-polish -R '^(wide_.*|bidirectional|k16_bidirectional)$' \
  --output-on-failure -j1 | tee "$result/tests.log"
python3 workstreams/spin_optimized/isa_check.py build-wide-polish --optional-only | tee "$result/isa.log"
