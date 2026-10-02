#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?source directory}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
kernels=${4:-/tmp/spin-packed-mixer-NnskDD}
if [[ $phase == build ]]; then
  # Match retained joint driver ISA flags; the separately compiled BCH object has GFNI.
  g++ -O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mtune=znver4 \
    -fno-strict-aliasing -fstack-usage -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -I"$reference/spin/src/kernels" -I"$reference/spin/include" \
    "$root/packed_driver.cpp" "$kernels/PackedMixer.o" /tmp/spin-bch-compare-2e2UPS/BchAvx512.o -o "$root/packed_driver"
  sha256sum "$root/packed_driver.cpp" "$root/packed_driver" "$kernels/PackedMixer.o"
  nm -S -C "$root/packed_driver" > "$root/packed_driver_symbols.txt"
  objdump -d -C "$root/packed_driver" > "$root/packed_driver_disassembly.txt"
  exit 0
fi
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
mkdir -p "$root/measurements"; logs=$(mktemp -d "$root/measurements/driver-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == check ]]; then
  for seed in 1 17; do
    for mode in 0 1 2; do
      for exponent in 14 20; do
        taskset -c 15 "$root/packed_driver" "$exponent" "$mode" "$seed" 0 > "$logs/check-$exponent-$mode-$seed.txt"
        cat "$logs/check-$exponent-$mode-$seed.txt"
      done
    done
  done
  exit 0
fi
calls=101; profile=0; repeats=2
if [[ $phase == screen ]]; then calls=31; profile=1; repeats=1; fi
if [[ $phase == profile ]]; then calls=31; profile=1; repeats=1; fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    for slot in 0 1 2 3; do
      mode=$slot;if ((rep%2==1)); then mode=$((3-slot)); fi
      if ((mode==3)); then
        taskset -c 15 /tmp/spin-fused-bch-oHiXdl/joint 20 4 2 9 "$seed" "$calls" "$profile" > "$logs/$seed-$rep-retained.csv"
        cat "$logs/$seed-$rep-retained.csv"
      else
        taskset -c 15 "$root/packed_driver" 20 "$mode" "$seed" "$calls" "$profile" > "$logs/$seed-$rep-$mode.csv"
        cat "$logs/$seed-$rep-$mode.csv"
      fi
    done
  done
done
