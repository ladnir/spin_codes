#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/forward_expansion
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-forward-expansion -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_INNER_PROBE=ON \
  -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON -DSPIN_FORWARD_SCHEDULE=dfs -DSPIN_FORWARD_EXPANSION=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
cmake --build build-forward-expansion --target spin_single_benchmark spin_bidirectional_test spin_k16_bidirectional_test spin_inner_probe_benchmark -j2 > "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-forward-expansion -R '^(bidirectional|k16_bidirectional)$' --output-on-failure -j1 | tee "$result/tests.log"
for repeat in 1 2 3; do
  variants=(control circuit)
  [[ $repeat != 2 ]] || variants=(circuit control)
  for m in 18 20; do
    for variant in "${variants[@]}"; do
      build=build-forward-schedule-dfs
      [[ $variant != circuit ]] || build=build-forward-expansion
      "$build/spin_single_benchmark" "$m" forward 101 256 > "$result/m$m-$variant-r$repeat.json"
    done
    for mode in 0 1 2; do
      build-forward-expansion/spin_inner_probe_benchmark "$m" 101 "$mode" > "$result/probe-m$m-mode$mode-r$repeat.csv"
    done
  done
done
