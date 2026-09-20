#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/k20-batch
mkdir -p "$result"
for mode in 6 7; do
 build="build-k20-$mode"
 cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" -DSPIN_K20_EXPERIMENT="$mode" > "$result/build-$mode.log" 2>&1
 cmake --build "$build" -j2 >> "$result/build-$mode.log" 2>&1
 flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/test-$mode.log"
done
for repeat in 1 2; do
 modes=(0 6 7)
 if [[ $repeat == 2 ]]; then modes=(7 6 0); fi
 for mode in "${modes[@]}"; do
  for kind in spin_benchmark spin_bidirectional_benchmark; do
   "build-k20-$mode/$kind" 20 auto 51 1 > "$result/$kind-mode$mode-r$repeat.json"
  done
 done
done
