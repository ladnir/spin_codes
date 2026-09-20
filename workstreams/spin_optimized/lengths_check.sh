#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/lengths
mkdir -p "$result"
for mode in control general masked sanitize avx2; do
  options=(-DSPIN_GENERAL_LENGTHS=ON -DSPIN_BCH_AVX512=ON -DSPIN_TEST_NO_AVX512=OFF -DSPIN_DIRECT_MAX_K=458752)
  if [[ $mode == control ]]; then options+=(-DSPIN_GENERAL_LENGTHS=OFF); fi
  if [[ $mode == masked ]]; then options+=(-DSPIN_TEST_NO_AVX512=ON); fi
  if [[ $mode == avx2 ]]; then options+=(-DSPIN_BCH_AVX512=OFF); fi
  if [[ $mode == sanitize ]]; then options+=('-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer'
      '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined'); fi
  cmake -S workstreams/spin_optimized -B "build-length-$mode" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON -DSPIN_FORWARD_FOUR=OFF "${options[@]}" > "$result/$mode-configure.log" 2>&1
  cmake --build "build-length-$mode" -j2 > "$result/$mode-build.log" 2>&1
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
    flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "build-length-$mode" --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
python3 workstreams/spin_optimized/isa_check.py build-length-general > "$result/isa.log"
for k in 8404992 67108864; do
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    build-length-general/spin_lengths_test "$k" > "$result/large-$k.log"
done
# Only now measure. No compilation or other benchmark runs concurrently.
for repeat in 1 2 3; do
  modes=(control general)
  if [[ $repeat == 2 ]]; then modes=(general control); fi
  for m in 16 18 20; do
    cfg=12819
    if [[ $m == 16 ]]; then cfg=6412r2; fi
    for mode in "${modes[@]}"; do
      "build-length-$mode/spin_benchmark" "$m" auto 51 1 0 0 "$cfg" > "$result/reg-$mode-m$m-r$repeat.json"
    done
  done
done
for k in 16384 32768 49152 65536 81920 98304 114688 131072 147456 163840 180224 196608 212992 229376 245760 262144 278528 393216 524288 786432 1032192 1048576 1064960; do
  build-length-general/spin_benchmark 20 auto 21 1 0 0 12819 "$k" > "$result/grid-k$k.json"
done
