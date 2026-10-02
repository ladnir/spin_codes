#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?source directory}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
objects=${4:-/tmp/spin-bch-compare-2e2UPS}
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni -mtune=znver4
 -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
 -I"$reference/spin/src/kernels" -I"$reference/spin/include")
if [[ $phase == build || $phase == relink ]]; then
  if [[ $phase == build ]]; then
    python3 "$root/packed_mixer_codegen.py" "$objects/BchCircuit.h" > "$root/PackedMixer.cpp"
    g++ "${flags[@]}" -fstack-usage -c "$root/PackedMixer.cpp" -o "$root/PackedMixer.o"
  fi
  g++ "${flags[@]}" "$root/packed_mixer_perf.cpp" "$objects/BchAvx512.o" "$objects/BchRestrict.o" "$root/PackedMixer.o" -o "$root/packed_mixer_perf"
  sha256sum "$root"/packed_mixer* "$root/PackedMixer.cpp" "$root/PackedMixer.o"
  exit 0
fi
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w 0 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/packed-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/packed_mixer_perf"
if [[ $phase == check ]]; then
  for seed in 1 17; do
    for mode in 0 1 2 3 4 5; do
      taskset -c 15 "$binary" check 14 "$mode" "$seed" 0 | tee "$logs/basis-$seed-$mode.txt"
      taskset -c 15 "$binary" full 14 "$mode" "$seed" 0 | tee "$logs/full-$seed-$mode.txt"
    done
  done
  exit 0
fi
if [[ $phase == check-relocated ]]; then
  for seed in 1 17; do
    for mode in 6 7; do
      taskset -c 15 "$binary" check 14 "$mode" "$seed" 0 | tee "$logs/basis-$seed-$mode.txt"
      taskset -c 15 "$binary" full 14 "$mode" "$seed" 0 | tee "$logs/full14-$seed-$mode.txt"
      taskset -c 15 "$binary" full 20 "$mode" "$seed" 0 | tee "$logs/full20-$seed-$mode.txt"
    done
  done
  exit 0
fi
case $phase in
 screen) kinds=(hot bulk full profile); modes=(0 1 2 3 4 5); seeds=(1); calls=31; repeats=1;;
 confirm) kinds=(full); modes=(0 2 3 4); seeds=(1 17); calls=101; repeats=2;;
 relocated) kinds=(full); modes=(0 6 7 4); seeds=(1 17); calls=101; repeats=2;;
 relocated-profile) kinds=(profile); modes=(0 6 7 4); seeds=(1 17); calls=31; repeats=1;;
 *) exit 2;;
esac
for kind in "${kinds[@]}"; do
  for seed in "${seeds[@]}"; do
    for ((rep=0;rep<repeats;++rep)); do
      for ((slot=0;slot<${#modes[@]};++slot)); do
        idx=$slot
        if ((rep%2==1)); then idx=$((${#modes[@]}-1-slot)); fi
        mode=${modes[$idx]}
        taskset -c 15 "$binary" "$kind" 20 "$mode" "$seed" "$calls" > "$logs/$kind-$seed-$rep-$mode.csv"
        cat "$logs/$kind-$seed-$rep-$mode.csv"
      done
    done
  done
done
