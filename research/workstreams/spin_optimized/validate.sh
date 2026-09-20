#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements
mkdir -p "$result"
for mode in on off masked; do
  options=(-DSPIN_BCH_AVX512=ON -DSPIN_TEST_NO_AVX512=OFF)
  if [[ $mode == off ]]; then options=(-DSPIN_BCH_AVX512=OFF); fi
  if [[ $mode == masked ]]; then options+=(-DSPIN_TEST_NO_AVX512=ON); fi
  cmake -S workstreams/spin_optimized -B "build-integrated-$mode" -DCMAKE_BUILD_TYPE=Release -DSPIN_BUILD_BENCHMARK=ON "${options[@]}"
  cmake --build "build-integrated-$mode" -j2
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "build-integrated-$mode" --output-on-failure -j1 | tee "$result/test-$mode.log"
done
# No builds or tests run alongside the measured calls. Each benchmark holds both locks.
for m in 16 18 20; do
  for repeat in 1 2 3; do
    backends=(avx2 auto)
    if [[ $repeat == 2 ]]; then backends=(auto avx2); fi
    for backend in "${backends[@]}"; do
      build-integrated-on/spin_benchmark "$m" "$backend" 101 1 > "$result/$backend-m$m-r$repeat.json"
    done
  done
done
find build-integrated-on/generated -type f -print0 | sort -z | xargs -0 sha256sum > "$result/generated.sha256"
sha256sum build-integrated-on/spin_benchmark > "$result/binary.sha256"
