#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
here=workstreams/inner_design/scheduling_20260919
result="$here/measurements"
# Acquire in the common order for all resource-intensive validation.
exec 8>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 8
exec 9>/tmp/bare-spin-benchmark.lock
flock -n 9
cmake -S "$here" -B build-sanitize -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
  -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined' > "$result/sanitize-configure.log"
cmake --build build-sanitize -j3 --target schedule_test pf16_test decode_pf_test > "$result/sanitize-build.log" 2>&1
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
  ctest --test-dir build-sanitize --output-on-failure -j1 -R '^(schedule_tails|pf16|decode_pf)$' | tee "$result/sanitizer.log"
