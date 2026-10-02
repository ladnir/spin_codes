#!/usr/bin/env bash
set -euo pipefail
cd "${1:?repository root}"
study=workstreams/inner_design/k16_parameters_20260919
result="$study/measurements/subspace_sanitizers"
[[ ! -e "$result" ]] || { echo "Refusing to overwrite sanitizer evidence" >&2; exit 1; }
mkdir -p "$result"
vendor=constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation/vendor
for shared in 0 1; do
 source="build-k16-64-12-${shared}-subspace_v1"
 build="build-k16-subspace-sanitize-${shared}"
 mkdir -p "$build"
 common=(-O2 -g -std=c++20 -march=x86-64-v3 -mtune=znver4 -mpclmul -mvpclmulqdq -fno-lto \
   -fsanitize=address,undefined -fno-omit-frame-pointer \
   -DSPIN_BCH_AVX512=1 -DSPIN_GROUPED_A=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
   -DW5_MASKED=0 -DW5_SHARED_EMISSION=0 -DW5_WRITE_PREFETCH=0 -DW5_PREFETCH=32 \
   -I"$source/generated" -I"$vendor")
 c++ "${common[@]}" -c "$source/generated/Spin.cpp" -o "$build/Spin.o"
 c++ "${common[@]}" -mavx512f -mavx512vl -c "$source/generated/Fast.cpp" -o "$build/Fast.o"
 bch=(build-final-sanitize/CMakeFiles/spin_half_transpose.dir/generated/generated/BchCircuit.cpp.o \
   build-final-sanitize/CMakeFiles/spin_bch512.dir/generated/generated/BchAvx512.cpp.o)
 c++ "${common[@]}" "$source/correctness.cpp" "$build/Spin.o" "$build/Fast.o" "${bch[@]}" -o "$build/test"
 ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=halt_on_error=1 \
   flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
   "$build/test" 2>&1 | tee "$result/${shared}.log"
done
