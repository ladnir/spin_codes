#!/usr/bin/env bash
# Host-specific experiment, not an installed-library ISA requirement.
set -euo pipefail
root=$(cd "$(dirname "$0")/../../.." && pwd)
output=${1:-"$root/out/permutation-bank-target"}
cpu=${2:-15}
mkdir -p "$output"
output=$(realpath "$output")
cxx=${CXX:-g++}
for feature in avx2 avx512f avx512vl avx512dq avx512bw avx512_vpopcntdq; do
    grep -qw "$feature" /proc/cpuinfo || { echo "Missing host feature: $feature" >&2; exit 1; }
done
cmake -S "$root/spin" -B "$output/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER="$cxx" -DSPIN_ENABLE_AVX512=ON -DSPIN_TUNE=znver4 \
    -DSPIN_BUILD_TESTS=OFF -DSPIN_BUILD_EXPERIMENTS=OFF
cmake --build "$output/build" --target spin -j2
flags=(-O3 -DNDEBUG -std=c++20 -mavx2 -mavx512f -mavx512vl -mavx512dq
       -mavx512bw -mavx512vpopcntdq -mtune=znver4 -DSPIN_BANK_AFFINE16=1
       -DSPIN_COMPOSED_SHIFT256_UNROLL=4 -DSPIN_WORKSPACE_ROUTING_OPT=1
       -DSPIN_BCH_AVX512=1 -DSPIN_COMPOSED_LOOKAHEAD=1 -DSPIN_COMPOSED_ROW_INDEX=1
       -DSPIN_COMPOSED_BYTE_ROUTE=1 -DSPIN_COMPOSED_ROTATE_KEYS=1
       -DSPIN_BANK_TIMING_REPS=101 -DSPIN_BANK_TIMING_WARMUP=20)
sources="$root/spin/experiments/permutation_bank"
"$cxx" "${flags[@]}" -I"$root/spin/include" "$sources/optimize.cpp" \
    "$output/build/libspin.a" -o "$output/build/spin_bank_target"
"$cxx" "${flags[@]}" -I"$root/spin/include" "$sources/composed_test.cpp" \
    -o "$output/build/spin_bank_target_routing"
"$cxx" "${flags[@]}" -I"$root/spin/include" "$sources/masks_test.cpp" \
    -o "$output/build/spin_bank_target_masks"
"$output/build/spin_bank_target_routing"
"$output/build/spin_bank_target_masks"
"$output/build/spin_bank_target" row-affine1 verify
"$output/build/spin_bank_target" row-rotate1 verify
# Each comparison takes all three shared locks; benchmarks run serially.
SPIN_COMPARE_PASSES=${SPIN_COMPARE_PASSES:-15} \
SPIN_COMPARE_MODES=fixed,row-affine1,row-rotate1 \
    bash "$sources/compare_optimize.sh" "$output" "$cpu" spin_bank_target target
python3 "$sources/summarize_optimize.py" "$output/results-bank-optimize-target"
SPIN_COMPARE_PASSES=1 SPIN_COMPARE_MODES=setup-composed \
    bash "$sources/compare_optimize.sh" "$output" "$cpu" spin_bank_target setup
