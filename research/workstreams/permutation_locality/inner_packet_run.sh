#!/usr/bin/env bash
# Isolated exact-map experiments. No simultaneous builds or benchmarks.
set -euo pipefail
root=$(realpath "${1:?fresh directory}")
phase=${2:?phase}
shift 2
prior=/tmp/spin-bch-tune-fN0NiP
ref=/tmp/spin-joint-9n57TT
base=/tmp/spin-bch-tune-vKm7EZ
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -n 7 || exit 75
logs=$(mktemp -d "$root/$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/inner-packet"
if [[ $phase == build ]]; then
  flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni -mprfchw -mtune=znver4
    -fno-strict-aliasing -fstack-usage -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
    -I"$root" -I"$ref/spin/src/kernels" -I"$ref/spin/include" -I"$prior" -I"$base")
  g++ "${flags[@]}" "$root/inner_packet_probe.cpp" "$prior/PackedCoeff1.o" \
    /tmp/spin-packed-mixer-NnskDD/PackedMixer.o /tmp/spin-bch-compare-2e2UPS/BchAvx512.o \
    -o "$binary" 2>&1 | tee "$logs/build.txt"
  sha256sum "$root/inner_packet_probe.cpp" "$root/InnerPacketMaps.h" "$root/InnerPacketFeedback.h" "$root/InnerWordUpdate.h" \
    "$root/T64Reference.h" "$binary" "$prior/PackedCoeff1.o" > "$logs/hashes.txt"
  nm -S -C "$binary" > "$logs/symbols.txt"
  exit 0
fi
if [[ $phase == check ]]; then
  for mode in "$@"; do
    for seed in 1 17; do
      for exponent in 14 20; do
        taskset -c 15 "$binary" "$mode" "$exponent" "$seed" 0 > "$logs/check-$mode-$exponent-$seed.txt"
        cat "$logs/check-$mode-$exponent-$seed.txt"
      done
    done
  done
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=31; repeats=1; profile=0
if [[ $phase == confirm ]]; then calls=101; repeats=2; fi
if [[ $phase == profile ]]; then profile=1; fi
modes=("$@")
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    for ((slot=0;slot<${#modes[@]};++slot)); do
      index=$slot; if ((rep)); then index=$((${#modes[@]}-1-slot)); fi
      mode=${modes[$index]}
      taskset -c 15 "$binary" "$mode" 20 "$seed" "$calls" "$profile" > "$logs/full-$mode-$seed-$rep.csv"
      cat "$logs/full-$mode-$seed-$rep.csv"
    done
  done
done
