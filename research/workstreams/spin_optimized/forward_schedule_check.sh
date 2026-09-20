#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/forward_schedule
mkdir -p "$result"
for mode in global dfs local16 local32; do
  build="build-forward-schedule-$mode"
  cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
    -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON -DSPIN_FORWARD_SCHEDULE="$mode" \
    -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/$mode-configure.log" 2>&1
  cmake --build "$build" --target spin_single_benchmark spin_bidirectional_test spin_k16_bidirectional_test -j2 > "$result/$mode-build.log" 2>&1
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "$build" -R '^(bidirectional|k16_bidirectional)$' --output-on-failure -j1 | tee "$result/$mode-tests.log"
  objdump -d -C "$build/CMakeFiles/spin_bidir512.dir/bidirectional/generated/BchForward512.cpp.o" > "$result/$mode-assembly.txt"
done
for repeat in 1 2 3; do
  modes=(global dfs local16 local32)
  [[ $repeat != 2 ]] || modes=(local32 local16 dfs global)
  for m in 16 18 20; do
    for mode in "${modes[@]}"; do
      "build-forward-schedule-$mode/spin_single_benchmark" "$m" forward 101 256 > "$result/m$m-$mode-r$repeat.json"
    done
  done
done
