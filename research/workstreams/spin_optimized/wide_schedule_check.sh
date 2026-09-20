#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/wide_schedule
mkdir -p "$result"
for mode in global dfs; do
  build="build-wide-schedule-$mode"
  cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_WIDE=ON -DSPIN_BUILD_BENCHMARK=ON \
    -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON -DSPIN_FORWARD_SCHEDULE=dfs \
    -DSPIN_WIDE_SCHEDULE="$mode" -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/$mode-configure.log" 2>&1
  cmake --build "$build" --target spin_wide_test spin_bidirectional_test spin_k16_bidirectional_test spin_wide_benchmark -j2 > "$result/$mode-build.log" 2>&1
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "$build" -R '^(wide_forward|bidirectional|k16_bidirectional)$' --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
for repeat in 1 2 3; do
  modes=(global dfs)
  [[ $repeat != 2 ]] || modes=(dfs global)
  for mode in "${modes[@]}"; do
    for m in 16 18 20; do
      cfg=12819; [[ $m != 16 ]] || cfg=6412r2
      for lanes in 1 2 4; do
        tile=256
        [[ $m != 20 || $lanes == 1 ]] || tile=512
        "build-wide-schedule-$mode/spin_wide_benchmark" "$m" 7 "$cfg" "$tile" "$lanes" > "$result/$mode-m$m-w$lanes-r$repeat.csv"
      done
    done
  done
done
