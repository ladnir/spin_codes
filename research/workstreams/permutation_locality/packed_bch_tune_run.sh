#!/usr/bin/env bash
# Exact-map packed GL32/BCH tuning; all compilation and timing is serial.
set -euo pipefail
root=$(realpath "${1:?fresh build directory}")
phase=${2:-check}
modes=${3:-"0 1 2 3 4"}
reference=/tmp/spin-joint-9n57TT
objects=/tmp/spin-bch-compare-2e2UPS
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
common=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mtune=znver4
  -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$reference/spin/src/kernels" -I"$reference/spin/include")
if [[ $phase == build ]]; then
  python3 "$root/packed_bch_tune_codegen.py" "$objects/BchCircuit.h" > "$root/PackedTune.cpp"
  g++ "${common[@]}" -fstack-usage -c "$root/packed_fused_r4.cpp" -o "$root/FusedDriver.o"
  g++ "${common[@]}" -mgfni -c "$root/packed_mixer_perf.cpp" -o "$root/MixerChecks.o"
  for mode in $modes; do
    g++ "${common[@]}" -mgfni -fstack-usage -DSPIN_PACKED_TUNE="$mode" \
      -c "$root/PackedTune.cpp" -o "$root/PackedTune$mode.o"
    g++ "$root/FusedDriver.o" "$root/PackedTune$mode.o" "$objects/BchAvx512.o" -o "$root/full-$mode"
    g++ "$root/MixerChecks.o" "$root/PackedTune$mode.o" "$objects/BchAvx512.o" "$objects/BchRestrict.o" -o "$root/check-$mode"
    sha256sum "$root/full-$mode" "$root/PackedTune$mode.o"
    nm -S -C "$root/full-$mode" > "$root/symbols-$mode.txt"
  done
  sha256sum "$root/packed_bch_tune_codegen.py" "$root/PackedTune.cpp" "$root/FusedDriver.o"
  exit 0
fi
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/bch-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == check ]]; then
  for mode in $modes; do
    for seed in 1 17; do
      taskset -c 15 "$root/check-$mode" check 14 4 "$seed" 0 > "$logs/basis-$mode-$seed.txt"
      cat "$logs/basis-$mode-$seed.txt"
      for exponent in 14 20; do
        taskset -c 15 "$root/full-$mode" "$exponent" "$seed" 0 1 > "$logs/full-$mode-$exponent-$seed.txt"
        cat "$logs/full-$mode-$exponent-$seed.txt"
      done
    done
  done
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=101; profile=0; repeats=2
if [[ $phase == screen ]]; then calls=31; repeats=1; fi
if [[ $phase == profile ]]; then calls=31; profile=1; repeats=1; fi
read -ra selected <<< "$modes"
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    for ((slot=0;slot<${#selected[@]};++slot)); do
      index=$slot
      if ((rep%2)); then index=$((${#selected[@]}-1-slot)); fi
      mode=${selected[$index]}
      taskset -c 15 "$root/full-$mode" 20 "$seed" "$calls" 1 "$profile" > "$logs/$seed-$rep-mode$mode.csv"
      cat "$logs/$seed-$rep-mode$mode.csv"
    done
  done
done
