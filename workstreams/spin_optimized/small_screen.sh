#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements
for mode in 0 1 2; do
 cmake -S workstreams/spin_optimized -B "build-small-$mode" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_SMALL_EXPERIMENT="$mode" \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin"
 cmake --build "build-small-$mode" -j2
 flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "build-small-$mode" --output-on-failure -j1 | tee "$result/small-$mode-tests.log"
done
for mode in 0 1 2; do
 for kind in spin_benchmark spin_bidirectional_benchmark; do
  "build-small-$mode/$kind" 16 auto 101 1 > "$result/screen-$mode-$kind.json"
 done
done
for tile in 32 64 128 256 512; do
 for layout in 0 1; do
  build-small-0/spin_benchmark 16 auto 101 1 "$tile" "$layout" > "$result/screen-tile$tile-layout$layout.json"
 done
done
for m in 16 18 20; do
 for backend in avx2 auto; do
  build-small-0/spin_bidirectional_benchmark "$m" "$backend" 101 1 > "$result/bidir-m$m-$backend.json"
 done
done
