#!/usr/bin/env bash
set -euo pipefail
root=${1:?root}
cd "$root"
result=workstreams/spin_optimized/measurements/bch-local-screen
mkdir -p "$result"
for mode in global local8 local16 local32 single; do
 build="build-bch-$mode"
 cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BCH_SCHEDULE="$mode" -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/build-$mode.log" 2>&1
 cmake --build "$build" -j2 >> "$result/build-$mode.log" 2>&1
 flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/test-$mode.log"
done
for repeat in 1 2; do
 modes=(global local8 local16 local32 single)
 if [[ $repeat == 2 ]]; then modes=(single local32 local16 local8 global); fi
 for mode in "${modes[@]}"; do
  for kind in spin_benchmark spin_bidirectional_benchmark; do
   "build-bch-$mode/$kind" 20 auto 51 1 > "$result/$kind-$mode-r$repeat.json"
  done
 done
done
for mode in global local8 local16 local32 single; do
 objdump -d -C "build-bch-$mode/CMakeFiles/spin_bch512.dir/generated/generated/BchAvx512.cpp.o" > "$result/assembly-$mode.txt"
done
