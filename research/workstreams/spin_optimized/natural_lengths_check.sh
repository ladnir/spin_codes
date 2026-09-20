#!/usr/bin/env bash
# Serial correctness matrix; no timed benchmarks are run by this script.
set -euo pipefail
root=${1:?repository root}
upstream=${2:-$root/hypercat/native/spin}
mode=${3:-release} # release, masked, avx2, sanitize
cd "$root"
result=workstreams/spin_optimized/measurements/natural_lengths
mkdir -p "$result"
options=(-DSPIN_TEST_NO_AVX512=OFF -DSPIN_BUILD_BENCHMARK=OFF)
case "$mode" in
  release) options+=(-DSPIN_BUILD_BENCHMARK=ON) ;;
  masked) options=(-DSPIN_TEST_NO_AVX512=ON) ;;
  avx2) options+=(-DSPIN_BCH_AVX512=OFF -DSPIN_FORWARD_FOUR=OFF
    -DSPIN_FORWARD_DIRECT=OFF -DSPIN_FORWARD_SCHEDULE=global) ;;
  sanitize) options+=('-DCMAKE_CXX_FLAGS=-fsanitize=address,undefined -fno-omit-frame-pointer'
    '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined') ;;
  *) echo "unknown validation mode: $mode" >&2;exit 1 ;;
esac
build="build-natural-$mode"
cmake -S workstreams/spin_optimized -B "$build" \
  -C workstreams/spin_optimized/cmake/ForwardRecommended.cmake \
  -DSPIN_BIDIRECTIONAL_SOURCE="$upstream" -DSPIN_TUNE=znver4 "${options[@]}" > "$result/$mode-configure.log" 2>&1
cmake --build "$build" --target spin_lengths_test spin_forward_lengths_test spin_wide_lengths_test \
  spin_bidirectional_test spin_k16_bidirectional_test spin_wide_test -j2 > "$result/$mode-build.log" 2>&1
if [[ $mode == release ]]; then
  cmake --build "$build" --target spin_single_benchmark spin_single_transpose_benchmark spin_wide_benchmark \
    -j2 >> "$result/$mode-build.log" 2>&1
fi
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "$build" -R '^(aligned_lengths|forward_natural_lengths|wide_.*|bidirectional|k16_bidirectional)$' \
  --output-on-failure -j1 | tee "$result/$mode-tests.log"
