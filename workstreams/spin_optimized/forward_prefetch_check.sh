#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/forward_prefetch
mkdir -p "$result"
build=build-forward-pf64
cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_WIDE=ON \
  -DSPIN_BUILD_BENCHMARK=ON -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_PREFETCH=64 \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
cmake --build "$build" --target spin_bidirectional_test spin_k16_bidirectional_test spin_wide_test spin_wide_benchmark -j2 > "$result/build.log" 2>&1
for mode in control four pf64; do
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "build-forward-$mode" -R '^(bidirectional|k16_bidirectional|wide_forward)$' --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
for repeat in 1 2 3; do
  modes=(four pf64)
  [[ $repeat != 2 ]] || modes=(pf64 four)
  for mode in "${modes[@]}"; do
    for m in 16 18 20; do
      cfg=12819; tile=512
      [[ $m != 16 ]] || { cfg=6412r2; tile=256; }
      for lanes in 1 2; do
        "build-forward-$mode/spin_wide_benchmark" "$m" 7 "$cfg" "$tile" "$lanes" > "$result/$mode-m$m-w$lanes-r$repeat.csv"
      done
    done
  done
done
