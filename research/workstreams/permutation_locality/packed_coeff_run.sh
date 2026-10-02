#!/usr/bin/env bash
# Exact coefficient-layout comparison; never runs benchmarks concurrently.
set -euo pipefail
root=$(realpath "${1:?build directory}")
phase=${2:-check}
tile=${3:-1}
selected=${4:-"3 0 2"}
reference=/tmp/spin-joint-9n57TT
objects=/tmp/spin-bch-compare-2e2UPS
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
binary="$root/coeff-$tile"
common=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mtune=znver4
  -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$reference/spin/src/kernels" -I"$reference/spin/include")
if [[ $phase == build ]]; then
  python3 "$root/packed_coeff_codegen.py" "$objects/BchCircuit.h" --tile-mode "$tile" > "$root/PackedCoeff$tile.cpp"
  g++ "${common[@]}" -mgfni -fstack-usage -c "$root/PackedCoeff$tile.cpp" -o "$root/PackedCoeff$tile.o"
  g++ "${common[@]}" -fstack-usage "$root/packed_coeff_r4.cpp" "$root/PackedCoeff$tile.o" \
    /tmp/spin-packed-mixer-NnskDD/PackedMixer.o "$objects/BchAvx512.o" -o "$binary"
  sha256sum "$root/packed_coeff_r4.cpp" "$root/packed_coeff_codegen.py" "$root/PackedCoeff$tile.cpp" "$binary"
  nm -S -C "$binary" > "$root/coeff-symbols-$tile.txt"
  exit 0
fi
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/coeff-$tile-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == tile-compare ]]; then
  for seed in 1 17; do
    for rep in 0 1; do
      order=(0 1 3); if ((rep)); then order=(3 1 0); fi
      for candidate in "${order[@]}"; do
        candidateTile=$candidate; mode=2
        if ((candidate==0)); then candidateTile=1; mode=3; fi
        taskset -c 15 "$root/coeff-$candidateTile" 20 "$seed" 101 "$mode" 0 > "$logs/$seed-$rep-candidate$candidate.csv"
        cat "$logs/$seed-$rep-candidate$candidate.csv"
      done
    done
  done
  exit 0
fi
if [[ $phase == check ]]; then
  for seed in 1 17; do
    taskset -c 15 "$binary" basis "$seed" > "$logs/basis-$seed.txt"
    cat "$logs/basis-$seed.txt"
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
    read -ra order <<< "$selected"
    for ((slot=0;slot<${#order[@]};++slot)); do
      index=$slot; if ((rep%2)); then index=$((${#order[@]}-1-slot)); fi
      mode=${order[$index]}
      taskset -c 15 "$binary" 20 "$seed" "$calls" "$mode" "$profile" > "$logs/$seed-$rep-mode$mode.csv"
      cat "$logs/$seed-$rep-mode$mode.csv"
    done
  done
done
