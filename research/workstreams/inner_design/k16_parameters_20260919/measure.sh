#!/usr/bin/env bash
set -euo pipefail
cd "${1:?repository root}"
study=workstreams/inner_design/k16_parameters_20260919
vendor=constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/implementation/vendor
result="$study/measurements/${3:?supply a fresh batch name}"
maps=${4:-"$study/measurements"}
[[ ! -e "$result" ]] || { echo "Refusing to overwrite batch $result" >&2; exit 1; }
mkdir -p "$result"
# Reuse the selected direct-routing/BCH build and its exact compiler flags.
# Each copy is a new experiment. No measured production binary is overwritten.
specs=${2:-"32_15 64_15 64_19"}
for spec in $specs; do
 t=${spec%_*}; s=${spec#*_}
 for shared in 0 1; do
  build="build-k16-${t}-${s}-${shared}-${3}"
  mkdir -p "$build"
  cp -a build-final-default/generated "$build/"
  python3 "$study/overlay.py" "$build" "$maps/Map${t}S${s}.h" "$t" "$s" "$shared" \
    workstreams/inner_design/asymmetric/bch256/weight5/implementation/correctness.cpp
  common=(-O3 -DNDEBUG -std=c++20 -march=x86-64-v3 -mtune=znver4 -mpclmul -mvpclmulqdq -fno-lto \
    -DSPIN_BCH_AVX512=1 -DSPIN_GROUPED_A=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -DW5_MASKED=0 -DW5_SHARED_EMISSION=0 -DW5_WRITE_PREFETCH=0 -DW5_PREFETCH=32 \
    -I"$build/generated" -I"$vendor")
  c++ "${common[@]}" -c "$build/generated/Spin.cpp" -o "$build/Spin.o"
  c++ "${common[@]}" -mavx512f -mavx512vl -c "$build/generated/Fast.cpp" -o "$build/Fast.o"
  bch=(build-final-default/CMakeFiles/spin_half_transpose.dir/generated/generated/BchCircuit.cpp.o \
    build-final-default/CMakeFiles/spin_bch512.dir/generated/generated/BchAvx512.cpp.o)
  c++ "${common[@]}" "$build/correctness.cpp" "$build/Spin.o" "$build/Fast.o" "${bch[@]}" -o "$build/test"
  c++ "${common[@]}" workstreams/spin_optimized/benchmark.cpp "$build/Spin.o" "$build/Fast.o" "${bch[@]}" -o "$build/bench"
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    "$build/test" | tee "$result/test-${spec}-${shared}.log"
 done
done
# No compilation or correctness jobs overlap these timings. The harness
# obtains both shared benchmark locks and refuses other running benchmarks.
for seed in 1 17; do
 for repeat in 1 2 3; do
  cells=(baseline)
  for spec in $specs; do cells+=("${spec}_0" "${spec}_1"); done
  if [[ $repeat == 2 ]]; then
   reversed=()
   for ((i=${#cells[@]}-1;i>=0;i--)); do reversed+=("${cells[i]}"); done
   cells=("${reversed[@]}")
  fi
  for cell in "${cells[@]}"; do
   binary=build-final-default/spin_benchmark
   if [[ $cell != baseline ]]; then binary="build-k16-${cell//_/-}-${3}/bench"; fi
   "$binary" 16 auto 101 "$seed" > "$result/perf-${cell}-s${seed}-r${repeat}.json"
  done
 done
done
find build-k16-*-${3}/generated -type f -print0 | sort -z | xargs -0 sha256sum > "$result/sources.sha256"
sha256sum build-k16-*-${3}/bench build-final-default/spin_benchmark > "$result/binaries.sha256"
c++ --version > "$result/compiler.txt"
cp "$study/measure.sh" "$result/measure.sh"
