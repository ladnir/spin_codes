#!/usr/bin/env bash
# Research-only exact-map R4 comparison. Every process runs serially.
set -euo pipefail
root=$(realpath "${1:?fresh build directory}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
kernels=${4:-/tmp/spin-packed-mixer-NnskDD}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
if [[ $phase == build ]]; then
  g++ -O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mtune=znver4 \
    -fno-strict-aliasing -fstack-usage -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -I"$reference/spin/src/kernels" -I"$reference/spin/include" \
    "$root/packed_fused_r4.cpp" "$kernels/PackedMixer.o" /tmp/spin-bch-compare-2e2UPS/BchAvx512.o \
    -o "$root/packed_fused_r4"
  sha256sum "$root/FusedR4.h" "$root/packed_fused_r4.cpp" "$root/packed_driver.cpp" \
    "$root/packed_fused_r4" "$kernels/PackedMixer.o"
  nm -S -C "$root/packed_fused_r4" > "$root/packed_fused_r4_symbols.txt"
  objdump -d -C "$root/packed_fused_r4" > "$root/packed_fused_r4_disassembly.txt"
  exit 0
fi
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/fused-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == check ]]; then
  for seed in 1 17; do
    for exponent in 14 20; do
      taskset -c 15 "$root/packed_fused_r4" "$exponent" "$seed" 0 > "$logs/check-$exponent-$seed.txt"
      cat "$logs/check-$exponent-$seed.txt"
    done
  done
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=101; profile=0; repeats=2
if [[ $phase == screen ]]; then calls=31; repeats=1; fi
if [[ $phase == profile ]]; then calls=31; profile=1; repeats=1; fi
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    if ((rep%2)); then order=(2 1 0); else order=(0 1 2); fi
    for mode in "${order[@]}"; do
      taskset -c 15 "$root/packed_fused_r4" 20 "$seed" "$calls" "$mode" "$profile" > "$logs/$seed-$rep-mode$mode.csv"
      cat "$logs/$seed-$rep-mode$mode.csv"
    done
  done
done
