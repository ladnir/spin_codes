#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
here=workstreams/inner_design/bch_20260919
result="$here/measurements"
exec 8>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 8
exec 9>/tmp/bare-spin-benchmark.lock
flock -n 9
cmake -S "$here" -B build-bch-sanitize -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_CXX_FLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
 -DCMAKE_EXE_LINKER_FLAGS='-fsanitize=address,undefined' > "$result/sanitize-configure.log"
cmake --build build-bch-sanitize -j3 --target packedshare4_test packedshare4_outer_bench packedshare4_tiles > "$result/sanitize-build.log" 2>&1
ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
 ctest --test-dir build-bch-sanitize --output-on-failure -j1 -R '^packedshare4(_basis|_tiles)?$' | tee "$result/sanitizer.log"
