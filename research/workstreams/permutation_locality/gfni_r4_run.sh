#!/usr/bin/env bash
# Research-only exact-map GFNI inner probe. All runs are serial.
set -euo pipefail
root=$(realpath "${1:?fresh build directory}")
phase=${2:-check}
reference=/tmp/spin-joint-9n57TT
prior=/tmp/spin-bch-tune-fN0NiP
objects=/tmp/spin-bch-compare-2e2UPS
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
binary="$root/gfni-r4"
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/gfni-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == build ]]; then
  g++ -O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni -mtune=znver4 \
    -fno-strict-aliasing -fstack-usage -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -I"$reference/spin/src/kernels" -I"$reference/spin/include" -I"$prior" \
    "$root/gfni_r4_probe.cpp" "$prior/PackedCoeff1.o" \
    /tmp/spin-packed-mixer-NnskDD/PackedMixer.o "$objects/BchAvx512.o" \
    -o "$binary" 2>&1 | tee "$logs/build.txt"
  sha256sum "$root/FusedR4Gfni.h" "$root/gfni_r4_probe.cpp" "$binary" \
    "$prior/FusedR4.h" "$prior/packed_driver.cpp" "$prior/PackedCoeff1.o" \
    /tmp/spin-packed-mixer-NnskDD/PackedMixer.o "$objects/BchAvx512.o" > "$logs/hashes.txt"
  nm -S -C "$binary" > "$root/symbols.txt"
  cat "$logs/hashes.txt"
  exit 0
fi
if [[ $phase == check ]]; then
  for seed in 1 17; do
    taskset -c 15 "$binary" core "$seed" > "$logs/core-$seed.txt"
    cat "$logs/core-$seed.txt"
    for exponent in 14 20; do
      taskset -c 15 "$binary" "$exponent" "$seed" 0 0 > "$logs/full-$exponent-$seed.txt"
      cat "$logs/full-$exponent-$seed.txt"
    done
  done
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=101; repeats=2; profile=0
if [[ $phase == screen ]]; then calls=31; repeats=1; fi
if [[ $phase == profile ]]; then calls=31; repeats=1; profile=1; fi
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    order=(0 1 2); if ((rep)); then order=(2 1 0); fi
    for mode in "${order[@]}"; do
      taskset -c 15 "$binary" 20 "$seed" "$calls" "$mode" "$profile" > "$logs/$seed-$rep-mode$mode.csv"
      cat "$logs/$seed-$rep-mode$mode.csv"
    done
  done
done
