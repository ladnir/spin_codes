#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=${2:-workstreams/spin_optimized/measurements/wide}
mkdir -p "$result"
for mode in masked sanitize; do
  options=(-DSPIN_TEST_NO_AVX512=ON)
  if [[ $mode == sanitize ]]; then
    options=(-DSPIN_TEST_NO_AVX512=OFF '-DCMAKE_CXX_FLAGS=-O1 -fsanitize=address,undefined -fno-omit-frame-pointer'
      '-DCMAKE_CXX_FLAGS_RELEASE=-O1 -DNDEBUG' '-DCMAKE_EXE_LINKER_FLAGS=-fsanitize=address,undefined')
  fi
  cmake -S workstreams/spin_optimized -B "build-wide-$mode" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_BUILD_WIDE=ON -DSPIN_BUILD_BENCHMARK=ON \
    -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" "${options[@]}" > "$result/$mode-configure.log" 2>&1
  cmake --build "build-wide-$mode" --target spin_wide_test spin_wide_benchmark -j2 > "$result/$mode-build.log" 2>&1
  ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
    flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    "build-wide-$mode/spin_wide_test" | tee "$result/$mode-tests.log"
  if [[ $mode == sanitize ]]; then
    # A correctness run of the assembly kernels too, not a performance sample.
    ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
      "build-wide-$mode/spin_wide_benchmark" 14 1 > "$result/sanitize-assembly.log"
    ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
      "build-wide-$mode/spin_wide_benchmark" 16 1 6412r2 > "$result/sanitize-k16-assembly.log"
  else
    set +e
    "build-wide-$mode/spin_wide_benchmark" 14 1 > "$result/masked-assembly.log" 2>&1
    status=$?
    set -e
    [[ $status == 77 ]]
  fi
done
