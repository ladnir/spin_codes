#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements
for mode in 0 1 2 3 4; do
 cmake -S workstreams/spin_optimized -B "build-small-$mode" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_SMALL_EXPERIMENT="$mode" \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin"
 cmake --build "build-small-$mode" -j2
 flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "build-small-$mode" --output-on-failure -j1 | tee "$result/final-small-$mode-tests.log"
done
for seed in 1 17; do
 for repeat in 1 2 3; do
  modes=(0 1 2 3 4)
  if [[ $repeat == 2 ]]; then modes=(4 3 2 1 0); fi
  for mode in "${modes[@]}"; do
   for kind in spin_benchmark spin_bidirectional_benchmark; do
    "build-small-$mode/$kind" 16 auto 101 "$seed" > "$result/small-$mode-$kind-s$seed-r$repeat.json"
   done
  done
 done
done
for m in 18 20; do
 for repeat in 1 2 3; do
  for mode in 0 1; do
   for kind in spin_benchmark spin_bidirectional_benchmark; do
    "build-small-$mode/$kind" "$m" auto 101 1 > "$result/large-$mode-$kind-m$m-r$repeat.json"
   done
  done
 done
done
for mode in 0 1 2 3 4; do
 find "build-small-$mode/generated" "build-small-$mode/bidirectional" -type f -print0 | sort -z | xargs -0 sha256sum > "$result/small-$mode-sources.sha256"
 sha256sum "build-small-$mode/spin_benchmark" "build-small-$mode/spin_bidirectional_benchmark" > "$result/small-$mode-binaries.sha256"
done
