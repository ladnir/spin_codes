#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
cmake -S workstreams/spin_optimized -B build-integrated-tuned -DCMAKE_BUILD_TYPE=Release \
 -DSPIN_BCH_AVX512=ON -DSPIN_BUILD_BENCHMARK=ON -DSPIN_TUNE=znver4
cmake --build build-integrated-tuned -j2
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir build-integrated-tuned --output-on-failure -j1 | tee workstreams/spin_optimized/measurements/test-tuned.log
for m in 16 18 20; do
 for repeat in 1 2 3; do
  backends=(avx2 auto)
  if [[ $repeat == 2 ]]; then backends=(auto avx2); fi
  for backend in "${backends[@]}"; do
   build-integrated-tuned/spin_benchmark "$m" "$backend" 101 1 > "workstreams/spin_optimized/measurements/tuned-$backend-m$m-r$repeat.json"
  done
 done
done
sha256sum build-integrated-tuned/spin_benchmark > workstreams/spin_optimized/measurements/binary-tuned.sha256
