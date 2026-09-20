#!/usr/bin/env bash
set -euo pipefail
root=${1:?root}
cd "$root"
result=workstreams/spin_optimized/measurements/bch-reload-screen
mkdir -p "$result"
build=build-bch-reload
cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
 -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
 -DSPIN_BCH_SCHEDULE=reload -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/build-reload.log" 2>&1
cmake --build "$build" -j2 >> "$result/build-reload.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/test-reload.log"
for repeat in 1 2; do
 modes=(global reload)
 if [[ $repeat == 2 ]]; then modes=(reload global); fi
 for mode in "${modes[@]}"; do
  for kind in spin_benchmark spin_bidirectional_benchmark; do
   "build-bch-$mode/$kind" 20 auto 51 1 > "$result/$kind-$mode-r$repeat.json"
  done
 done
done
objdump -d -C build-bch-reload/CMakeFiles/spin_bch512.dir/generated/generated/BchAvx512.cpp.o > "$result/assembly-reload.txt"
