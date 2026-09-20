#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/generic
mkdir -p "$result"
for mode in release avx2 sanitize; do
  options=(-DSPIN_BCH_AVX512=ON)
  if [[ $mode == avx2 ]]; then options=(-DSPIN_BCH_AVX512=OFF); fi
  if [[ $mode == sanitize ]]; then options+=('-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer'
    '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined'); fi
  cmake -S workstreams/spin_optimized -B "build-generic-$mode" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_FORWARD_FOUR=OFF "${options[@]}" > "$result/$mode-configure.log" 2>&1
  cmake --build "build-generic-$mode" -j2 > "$result/$mode-build.log" 2>&1
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
    flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "build-generic-$mode" --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
python3 workstreams/spin_optimized/isa_check.py build-generic-release > "$result/isa.log"
# Original optimized path remains a separate executable/kernel; serial regression.
for repeat in 1 2 3; do
  for m in 16 18 20; do
    cfg=12819
    if [[ $m == 16 ]]; then cfg=6412r2; fi
    for build in build-length-general build-generic-release; do
      "$build/spin_benchmark" "$m" auto 51 1 0 0 "$cfg" > "$result/$build-m$m-r$repeat.json"
    done
  done
done
