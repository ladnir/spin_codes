#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
schedule=${2:-global}
cd "$root"
result=workstreams/spin_optimized/measurements/single_direct
prefix=build-single-direct
if [[ $schedule != global ]]; then
  result="workstreams/spin_optimized/measurements/forward_schedule/safety-$schedule"
  prefix="build-forward-schedule-$schedule"
fi
mkdir -p "$result"
for mode in masked sanitize; do
  options=(-DSPIN_TEST_NO_AVX512=ON)
  if [[ $mode == sanitize ]]; then
    options=(-DSPIN_TEST_NO_AVX512=OFF '-DCMAKE_CXX_FLAGS=-O1 -fsanitize=address,undefined -fno-omit-frame-pointer'
      '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined')
  fi
  build="$prefix-$mode"
  cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_DIRECT=ON -DSPIN_FORWARD_SCHEDULE="$schedule" \
    -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" "${options[@]}" > "$result/$mode-configure.log" 2>&1
  cmake --build "$build" --target spin_bidirectional_test spin_k16_bidirectional_test -j2 > "$result/$mode-build.log" 2>&1
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
    flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "$build" -R '^(bidirectional|k16_bidirectional)$' --output-on-failure -j1 | tee "$result/$mode-tests.log"
done
python3 workstreams/spin_optimized/isa_check.py "$prefix" --optional-only | tee "$result/isa.log"
