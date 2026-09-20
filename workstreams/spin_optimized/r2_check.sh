#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/r2-integration
mkdir -p "$result"
for mode in direct tiled off masked sanitize; do
 build="build-k16-r2-integrated${mode/direct/}"
 if [[ $mode != direct ]]; then build="build-k16-r2-integrated-$mode"; fi
 options=(-DSPIN_BCH_AVX512=ON -DSPIN_K16_DIRECT=ON -DSPIN_BUILD_BENCHMARK=OFF)
 if [[ $mode == direct ]]; then options+=(-DSPIN_BUILD_BENCHMARK=ON); fi
 if [[ $mode == tiled ]]; then options+=(-DSPIN_K16_DIRECT=OFF); fi
 if [[ $mode == off ]]; then options+=(-DSPIN_BCH_AVX512=OFF); fi
 if [[ $mode == masked ]]; then options+=(-DSPIN_TEST_NO_AVX512=ON); fi
 if [[ $mode == sanitize ]]; then
  options+=('-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined')
 fi
 cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_TUNE=znver4 -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" "${options[@]}" > "$result/$mode-build.log" 2>&1
 cmake --build "$build" -j2 >> "$result/$mode-build.log" 2>&1
 ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/$mode-tests.log"
 cp "$build/Testing/Temporary/LastTest.log" "$result/$mode-details.log"
done
python3 workstreams/spin_optimized/isa_check.py build-k16-r2-integrated | tee "$result/isa.log"
for seed in 1 17; do
 for repeat in 1 2 3; do
  configs=(12819 6412 6412r2)
  if [[ $repeat == 2 ]]; then configs=(6412r2 6412 12819); fi
  for config in "${configs[@]}"; do
   for kind in spin_benchmark spin_bidirectional_benchmark; do
    "build-k16-r2-integrated/$kind" 16 auto 101 "$seed" 0 0 "$config" > "$result/$kind-$config-s$seed-r$repeat.json"
   done
  done
 done
done
find build-k16-r2-integrated/generated build-k16-r2-integrated/bidirectional -type f -print0 | sort -z | xargs -0 sha256sum > "$result/sources.sha256"
sha256sum build-k16-r2-integrated/spin_benchmark build-k16-r2-integrated/spin_bidirectional_benchmark > "$result/binaries.sha256"
