#!/usr/bin/env bash
# Correctness builds first; never overlap a benchmark with builds or other runs.
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements
mkdir -p "$result"
for mode in tiled off masked sanitize; do
 options=(-DSPIN_BCH_AVX512=ON -DSPIN_K16_DIRECT=ON)
 if [[ $mode == tiled ]]; then options+=(-DSPIN_K16_DIRECT=OFF); fi
 if [[ $mode == off ]]; then options+=(-DSPIN_BCH_AVX512=OFF); fi
 if [[ $mode == masked ]]; then options+=(-DSPIN_TEST_NO_AVX512=ON); fi
 if [[ $mode == sanitize ]]; then
  options+=('-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined')
 fi
 build="build-k16-integrated-$mode"
 cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=OFF \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" "${options[@]}" > "$result/k16-$mode-build.log" 2>&1
 cmake --build "$build" -j2 >> "$result/k16-$mode-build.log" 2>&1
 ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/k16-$mode-tests.log"
done
python3 workstreams/spin_optimized/isa_check.py build-k16-integrated | tee "$result/k16-isa.log"
for seed in 1 17; do
 for repeat in 1 2 3; do
  configurations=(12819 6412)
  if [[ $repeat == 2 ]]; then configurations=(6412 12819); fi
  for cfg in "${configurations[@]}"; do
   for kind in spin_benchmark spin_bidirectional_benchmark; do
    "build-k16-integrated/$kind" 16 auto 101 "$seed" 0 0 "$cfg" > "$result/k16-integrated-$kind-$cfg-s$seed-r$repeat.json"
   done
  done
 done
done
find build-k16-integrated/generated build-k16-integrated/bidirectional -type f -print0 | sort -z | xargs -0 sha256sum > "$result/k16-integrated-sources.sha256"
sha256sum build-k16-integrated/spin_benchmark build-k16-integrated/spin_bidirectional_benchmark > "$result/k16-integrated-binaries.sha256"
