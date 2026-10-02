#!/usr/bin/env bash
# Complete S16 research encoder: serial build, validation and comparison.
set -euo pipefail
root=$(realpath "${1:?fresh build directory}")
phase=${2:-check}
reference=/tmp/spin-joint-9n57TT
prior=/tmp/spin-bch-tune-fN0NiP
control=/tmp/spin-bch-tune-Q1wko3/gfni-r4
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
binary="$root/gfni-s16"
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/s16-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == build ]]; then
  g++ -O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni -mtune=znver4 \
    -fno-strict-aliasing -fstack-usage -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -I"$reference/spin/src/kernels" -I"$reference/spin/include" -I"$prior" \
    "$root/gfni_s16_probe.cpp" "$prior/PackedCoeff1.o" \
    /tmp/spin-packed-mixer-NnskDD/PackedMixer.o /tmp/spin-bch-compare-2e2UPS/BchAvx512.o \
    -o "$binary" 2>&1 | tee "$logs/build.txt"
  sha256sum "$root/GfniS16.h" "$root/FusedR4Gfni.h" "$root/gfni_s16_probe.cpp" "$binary" \
    "$prior/packed_driver.cpp" "$prior/PackedCoeff1.o" "$control" \
    /tmp/spin-packed-mixer-NnskDD/PackedMixer.o /tmp/spin-bch-compare-2e2UPS/BchAvx512.o > "$logs/hashes.txt"
  nm -S -C "$binary" > "$root/symbols.txt"
  cat "$logs/hashes.txt"
  exit 0
fi
if [[ $phase == check ]]; then
  for seed in 1 17; do
    for feedback in 0 1; do
      for refresh in 0 8 16; do
        taskset -c 15 "$binary" 14 "$seed" 0 "$refresh" "$feedback" > "$logs/full-14-$seed-$refresh-$feedback.txt"
        cat "$logs/full-14-$seed-$refresh-$feedback.txt"
      done
    done
  done
  taskset -c 15 "$binary" 20 1 0 0 1 > "$logs/full-20-1-0-1.txt"
  cat "$logs/full-20-1-0-1.txt"
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=31; repeats=1; profile=0
if [[ $phase == confirm ]]; then calls=101; repeats=2; fi
if [[ $phase == profile ]]; then profile=1; fi
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    order=(control bch16 weight5); if ((rep)); then order=(weight5 bch16 control); fi
    for kind in "${order[@]}"; do
      if [[ $kind == control ]]; then
        taskset -c 15 "$control" 20 "$seed" "$calls" 1 "$profile" > "$logs/$seed-$rep-control.csv"
      else
        feedback=1; if [[ $kind == weight5 ]]; then feedback=0; fi
        taskset -c 15 "$binary" 20 "$seed" "$calls" 0 "$feedback" "$profile" > "$logs/$seed-$rep-$kind.csv"
      fi
      cat "$logs/$seed-$rep-$kind.csv"
    done
  done
done
