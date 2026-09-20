#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements
for mode in control default off masked sanitize; do
 options=(-DSPIN_BCH_AVX512=ON -DSPIN_K16_DIRECT=ON)
 if [[ $mode == control ]]; then options+=(-DSPIN_K16_DIRECT=OFF); fi
 if [[ $mode == off ]]; then options+=(-DSPIN_BCH_AVX512=OFF); fi
 if [[ $mode == masked ]]; then options+=(-DSPIN_TEST_NO_AVX512=ON); fi
 if [[ $mode == sanitize ]]; then options+=(-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined\ -fno-omit-frame-pointer -DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined); fi
 cmake -S workstreams/spin_optimized -B "build-final-$mode" -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_SMALL_EXPERIMENT=auto -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" "${options[@]}"
 cmake --build "build-final-$mode" -j2
 ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "build-final-$mode" --output-on-failure -j1 | tee "$result/production-$mode-tests.log"
done
python3 workstreams/spin_optimized/isa_check.py build-final-default | tee "$result/production-isa.log"
for m in 16 18 20; do
 seeds=(1)
 if [[ $m == 16 ]]; then seeds=(1 17); fi
 for seed in "${seeds[@]}"; do
  for repeat in 1 2 3; do
   modes=(control default)
   if [[ $repeat == 2 ]]; then modes=(default control); fi
   for mode in "${modes[@]}"; do
    for kind in spin_benchmark spin_bidirectional_benchmark; do
     "build-final-$mode/$kind" "$m" auto 101 "$seed" > "$result/production-$mode-$kind-m$m-s$seed-r$repeat.json"
    done
   done
  done
 done
done
find build-final-default/generated build-final-default/bidirectional -type f -print0 | sort -z | xargs -0 sha256sum > "$result/production-sources.sha256"
sha256sum build-final-default/spin_benchmark build-final-default/spin_bidirectional_benchmark > "$result/production-binaries.sha256"
