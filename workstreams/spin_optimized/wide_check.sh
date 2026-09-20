#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/wide
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-wide -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_WIDE=ON -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
cmake --build build-wide -j2 > "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-wide --output-on-failure -j1 | tee "$result/tests.log"
python3 workstreams/spin_optimized/isa_check.py build-wide | tee "$result/isa.log"
for repeat in 1 2 3; do
  for m in 14 16 18 20; do
    build-wide/spin_wide_benchmark "$m" 7 > "$result/m$m-r$repeat.csv"
  done
done
sha256sum build-wide/spin_wide_benchmark > "$result/binary.sha256"
cp build-wide/bidirectional/WIDE_IMPORT.json "$result/"
cp build-wide/bidirectional/wide/generated/WIDE_MANIFEST.json "$result/"
