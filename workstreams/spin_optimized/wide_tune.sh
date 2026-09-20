#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/wide_tuning
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-wide -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_WIDE=ON -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
cmake --build build-wide -j2 > "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-wide --output-on-failure -j1 | tee "$result/tests.log"
python3 workstreams/spin_optimized/isa_check.py build-wide | tee "$result/isa.log"
for repeat in 1 2 3; do
  for cfg in 12819 6412r2; do
    build-wide/spin_wide_benchmark 16 7 "$cfg" > "$result/k16-$cfg-r$repeat.csv"
  done
done
for repeat in 1 2; do
  tiles=(8 16 32 64 128 256 512 1024 2048 4096)
  if [[ $repeat == 2 ]]; then tiles=(4096 2048 1024 512 256 128 64 32 16 8); fi
  for tile in "${tiles[@]}"; do
    build-wide/spin_wide_benchmark 20 3 12819 "$tile" 2 > "$result/screen-t$tile-r$repeat.csv"
  done
done
sha256sum build-wide/spin_wide_benchmark > "$result/binary.sha256"
cp build-wide/bidirectional/WIDE_IMPORT.json "$result/"
cp build-wide/bidirectional/wide/generated/WIDE_MANIFEST.json "$result/"
