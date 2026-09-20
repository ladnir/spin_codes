#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/single_direct
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-single-direct -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/configure.log" 2>&1
cmake --build build-single-direct --target spin_single_benchmark spin_bidirectional_test spin_k16_bidirectional_test -j2 > "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-single-direct -R '^(bidirectional|k16_bidirectional)$' --output-on-failure -j1 | tee "$result/tests.log"
for repeat in 1 2 3; do
  modes=(timing direct)
  [[ $repeat != 2 ]] || modes=(direct timing)
  for m in 16 18; do
    tile=256; [[ $m != 16 ]] || tile=64
    for mode in "${modes[@]}"; do
      "build-single-$mode/spin_single_benchmark" "$m" forward 101 "$tile" > "$result/m$m-$mode-r$repeat.json"
    done
  done
done
