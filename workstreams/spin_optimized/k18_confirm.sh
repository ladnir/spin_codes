#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/k18-confirm
mkdir -p "$result"
for mode in control selected off masked sanitize; do
 build="build-k18-final-$mode"
 options=(-DSPIN_BCH_AVX512=ON -DSPIN_K18_DIRECT=ON -DSPIN_K18_EXPERIMENT=auto -DSPIN_BUILD_BENCHMARK=OFF)
 if [[ $mode == selected || $mode == control ]]; then options+=(-DSPIN_BUILD_BENCHMARK=ON); fi
 if [[ $mode == control ]]; then options+=(-DSPIN_K18_DIRECT=OFF); fi
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
python3 workstreams/spin_optimized/isa_check.py build-k18-final-selected | tee "$result/isa.log"
for seed in 1 17; do
 for repeat in 1 2 3; do
  cases=(control selected)
  if [[ $repeat == 2 ]]; then cases=(selected control); fi
  for mode in "${cases[@]}"; do
   for kind in spin_benchmark spin_bidirectional_benchmark; do
    "build-k18-final-$mode/$kind" 18 auto 101 "$seed" > "$result/$mode-$kind-m18-s$seed-r$repeat.json"
   done
  done
 done
done
# Sizes outside the specialized path are controls, not new optimization claims.
for m in 16 20; do
 for mode in control selected; do
  for kind in spin_benchmark spin_bidirectional_benchmark; do
   "build-k18-final-$mode/$kind" "$m" auto 101 1 > "$result/$mode-$kind-m$m-s1-r1.json"
  done
 done
done
find build-k18-final-selected/generated build-k18-final-selected/bidirectional -type f -print0 | sort -z | xargs -0 sha256sum > "$result/sources.sha256"
sha256sum build-k18-final-selected/spin_benchmark build-k18-final-selected/spin_bidirectional_benchmark > "$result/binaries.sha256"
python3 workstreams/spin_optimized/k18_summary.py "$result" | tee "$result/summary.txt"
