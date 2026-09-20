#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/forward_safety
mkdir -p "$result"
for mode in masked sanitize; do
  options=(-DSPIN_TEST_NO_AVX512=ON)
  if [[ $mode == sanitize ]]; then
    options=(-DSPIN_TEST_NO_AVX512=OFF '-DCMAKE_CXX_FLAGS=-O1 -fsanitize=address,undefined -fno-omit-frame-pointer'
      '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined')
  fi
  build="build-forward-$mode"
  cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_BUILD_WIDE=ON -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_PREFETCH=0 \
    -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" "${options[@]}" > "$result/$mode-configure.log" 2>&1
  cmake --build "$build" --target spin_bidirectional_test spin_k16_bidirectional_test spin_wide_test -j2 > "$result/$mode-build.log" 2>&1
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
    flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "$build" -R '^(bidirectional|k16_bidirectional|wide_forward)$' --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
