#!/usr/bin/env bash
# Validate and compare the optional four-row forward kernel. Never time in parallel.
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/forward_four
mkdir -p "$result"
for mode in control four; do
  option=OFF
  [[ $mode != four ]] || option=ON
  build="build-forward-$mode"
  cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_WIDE=ON \
    -DSPIN_BUILD_BENCHMARK=ON -DSPIN_FORWARD_FOUR="$option" \
    -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/$mode-configure.log" 2>&1
  cmake --build "$build" --target spin_bidirectional_test spin_k16_bidirectional_test spin_wide_test spin_wide_benchmark -j2 > "$result/$mode-build.log" 2>&1
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "$build" -R '^(bidirectional|k16_bidirectional|wide_forward)$' --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
for repeat in 1 2 3; do
  modes=(control four)
  [[ $repeat != 2 ]] || modes=(four control)
  for mode in "${modes[@]}"; do
    for m in 16 18 20; do
      cfg=12819; tile=512
      [[ $m != 16 ]] || { cfg=6412r2; tile=256; }
      "build-forward-$mode/spin_wide_benchmark" "$m" 7 "$cfg" "$tile" 1 > "$result/$mode-m$m-r$repeat.csv"
    done
  done
done
