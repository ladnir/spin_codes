#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/lengths
mkdir -p "$result"
cmake -S workstreams/spin_optimized -B build-length-no-range -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_FORWARD_FOUR=OFF \
  -DSPIN_GENERAL_LENGTHS=ON -DSPIN_BCH_AVX512=ON -DSPIN_DIRECT_MAX_K=0 > "$result/no-range-configure.log" 2>&1
cmake --build build-length-no-range -j2 > "$result/no-range-build.log" 2>&1
# Deliberately test beyond the selected 458752 crossover to show the loss at 2^19.
cmake -S workstreams/spin_optimized -B build-length-range -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_FORWARD_FOUR=OFF \
  -DSPIN_GENERAL_LENGTHS=ON -DSPIN_BCH_AVX512=ON -DSPIN_DIRECT_MAX_K=524288 > "$result/range-configure.log" 2>&1
cmake --build build-length-range -j2 > "$result/range-build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-length-range --output-on-failure -j1 | tee "$result/range-tests.log"
for repeat in 1 2 3; do
  for k in 16384 49152 65536 81920 163840 245760 262144 278528 327680 393216 458752 524288 540672; do
    for mode in general range; do
      build="build-length-$mode"
      if [[ $mode == general ]]; then build=build-length-no-range; fi
      "$build/spin_benchmark" 20 auto 31 1 0 0 12819 "$k" > "$result/range-$mode-k$k-r$repeat.json"
    done
  done
done
