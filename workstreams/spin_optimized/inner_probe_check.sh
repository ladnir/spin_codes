#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/inner_probe_isolated
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-inner-probe -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_INNER_PROBE=ON \
  -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON -DSPIN_FORWARD_SCHEDULE=dfs \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
cmake --build build-inner-probe --target spin_inner_probe_benchmark -j2 > "$result/build.log" 2>&1
for repeat in 1 2 3; do
  for m in 18 20; do
    modes=(0 1 2)
    [[ $repeat != 2 ]] || modes=(2 1 0)
    for mode in "${modes[@]}"; do
      build-inner-probe/spin_inner_probe_benchmark "$m" 101 "$mode" > "$result/m$m-mode$mode-r$repeat.csv"
    done
  done
done
