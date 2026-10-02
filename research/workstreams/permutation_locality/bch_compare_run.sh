#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?comparison source directory}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mtune=znver4
       -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
       -I"$reference/spin/src/kernels" -I"$reference/spin/include")
if [[ $phase == build ]]; then
  python3 "$root/schedule_bch.py" "$root/BchAvx512.cpp" "$root/BchCircuit.h" Restrict > "$root/BchRestrict.cpp"
  python3 "$root/gfni_bch.py" "$root/BchCircuit.h" > "$root/BchGfni.cpp"
  for source in BchAvx512 BchRestrict BchGfni; do
    g++ "${flags[@]}" -mgfni -fstack-usage -c "$root/$source.cpp" -o "$root/$source.o"
  done
  g++ "${flags[@]}" "$root/bch_compare.cpp" "$root/BchAvx512.o" "$root/BchRestrict.o" "$root/BchGfni.o" -o "$root/bch_compare"
  sha256sum "$root"/*.cpp "$root"/*.o "$root/bch_compare"
  nm -S -C --size-sort "$root"/Bch*.o > "$root/symbol-sizes.txt"
  objdump -d -C "$root"/Bch*.o > "$root/disassembly.txt"
  exit 0
fi
if [[ $phase == relink ]]; then
  g++ "${flags[@]}" "$root/bch_compare.cpp" "$root/BchAvx512.o" "$root/BchRestrict.o" "$root/BchGfni.o" -o "$root/bch_compare"
  sha256sum "$root/bch_compare.cpp" "$root/bch_compare"
  exit 0
fi
if [[ $phase == build-algebraic ]]; then
  python3 "$root/bch_algebraic.py" --header "$root/BchCircuit.h" --leaf 32 --output "$root/BchAlgebraic.cpp"
  g++ "${flags[@]}" -mgfni -fstack-usage -c "$root/BchAlgebraic.cpp" -o "$root/BchAlgebraic.o"
  g++ "${flags[@]}" "$root/bch_algebraic_perf.cpp" "$root/BchAvx512.o" "$root/BchRestrict.o" "$root/BchGfni.o" "$root/BchAlgebraic.o" -o "$root/bch_algebraic_perf"
  sha256sum "$root/BchAlgebraic.cpp" "$root/BchAlgebraic.o" "$root/bch_algebraic_perf"
  nm -S -C --size-sort "$root/BchAlgebraic.o" > "$root/algebraic-symbol-sizes.txt"
  exit 0
fi
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/bch-compare-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/bch_compare"
if [[ $phase == algebraic ]]; then
  "$root/bch_algebraic_perf" check | tee "$logs/basis.txt"
  "$root/bch_algebraic_perf" full exact 17 0 > "$logs/full-check.txt"
  for seed in 1 17; do
    for kind in hot bulk full; do
      size=20; if [[ $kind == hot ]]; then size=1; fi
      for backend in 0 1 2; do
        taskset -c 15 "$binary" "$kind" "$size" "$backend" "$seed" 31 > "$logs/$kind-$seed-$backend.csv"
        cat "$logs/$kind-$seed-$backend.csv"
      done
      taskset -c 15 "$root/bch_algebraic_perf" "$kind" exact "$seed" 31 > "$logs/$kind-$seed-exact.csv"
      cat "$logs/$kind-$seed-exact.csv"
      if [[ $kind != full ]]; then
        taskset -c 15 "$root/bch_algebraic_perf" "$kind" raw "$seed" 31 > "$logs/$kind-$seed-raw.csv"
        cat "$logs/$kind-$seed-raw.csv"
      fi
    done
  done
  exit 0
fi
if [[ $phase == legacy-control ]]; then
  for seed in 1 17; do
    taskset -c 15 /tmp/spin-fused-bch-oHiXdl/joint 20 4 2 9 "$seed" 101 1 > "$logs/$seed-legacy.csv"
    cat "$logs/$seed-legacy.csv"
    taskset -c 15 "$binary" full-profile 20 3 "$seed" 101 > "$logs/$seed-fresh.csv"
    cat "$logs/$seed-fresh.csv"
  done
  exit 0
fi
if [[ $phase == control ]]; then
  for seed in 1 17; do
    for slot in 0 1 2 3; do
      backend=$((2+(slot/2==slot%2?0:1)))
      taskset -c 15 "$binary" full-profile 20 "$backend" "$seed" 101 > "$logs/$seed-$slot-$backend.csv"
      cat "$logs/$seed-$slot-$backend.csv"
    done
  done
  exit 0
fi
if [[ $phase == check ]]; then
  "$binary" check | tee "$logs/basis.txt"
  for backend in 0 1 2; do
    for seed in 1 17; do
      "$binary" full 14 "$backend" "$seed" 0 > "$logs/full-14-$backend-$seed.txt"
    done
    "$binary" full 20 "$backend" 17 0 > "$logs/full-20-$backend.txt"
  done
  printf 'Full-code/suffix checks pass for all backends.\n'
  exit 0
fi
case $phase in
  hot) kinds=(hot); sizes=(1 4 16); seeds=(1 17); calls=101;;
  bulk) kinds=(bulk); sizes=(16 18 20 22); seeds=(1 17); calls=101;;
  full) kinds=(full); sizes=(16 18 20 22); seeds=(1 17); calls=101;;
  profile) kinds=(full-profile); sizes=(16 18 20 22); seeds=(1 17); calls=31;;
  *) exit 2;;
esac
# Six balanced process orders at each size, using two seeds and three rotations.
for kind in "${kinds[@]}"; do
  for size in "${sizes[@]}"; do
    for seed in "${seeds[@]}"; do
      for repetition in 0 1 2; do
        for slot in 0 1 2; do
          backend=$(((repetition+slot)%3))
          if [[ $seed == 17 ]]; then backend=$((2-backend)); fi
          taskset -c 15 "$binary" "$kind" "$size" "$backend" "$seed" "$calls" > "$logs/$kind-$size-$seed-$repetition-$backend.csv"
          cat "$logs/$kind-$size-$seed-$repetition-$backend.csv"
        done
      done
    done
  done
done
sha256sum "$binary"
