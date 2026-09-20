#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
exec 8>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 8
exec 9>/tmp/bare-spin-benchmark.lock
flock -n 9
cmake -S workstreams/spin_optimized -B build-integrated-sanitize -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BCH_AVX512=ON -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
  -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined'
cmake --build build-integrated-sanitize -j2
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build-integrated-sanitize --output-on-failure -j1 | tee workstreams/spin_optimized/measurements/sanitizer.log
