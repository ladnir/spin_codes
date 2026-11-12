#!/bin/bash
# Build, test, and time the NEON and portable Turbo tiers on a 64-bit ARM machine
# (e.g. Apple Silicon with Xcode command-line tools and CMake, or Linux
# AArch64). Run from the spin/ directory: tools/arm_check.sh [build-dir]
# Benchmarks run serially after the tests; close other heavy work first.
set -euo pipefail
B=${1:-build-arm}
cmake -S . -B "$B" -DCMAKE_BUILD_TYPE=Release -DSPIN_BUILD_TESTS=ON -DSPIN_BUILD_BENCHMARKS=ON
cmake --build "$B" -j 8
ctest --test-dir "$B" --output-on-failure -j1
uname -a
for K in 65536 262144 1048576; do
  R=21; [ "$K" = 1048576 ] && R=9
  "$B"/spin_turbo_backends_bench neon "$K" "$R"
  "$B"/spin_turbo_backends_bench portable "$K" "$R"
done
for K in 65536 262144 1048576; do
  R=21; [ "$K" = 1048576 ] && R=9
  "$B"/spin_setup_once_bench turbo "$K" "$R" neon
  "$B"/spin_setup_once_bench turbobank "$K" "$R" neon
done
"$B"/spin_turbo_backends_bench reference 65536 5
