#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?directory containing joint.cpp and headers}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
if [[ $phase == build || $phase == build-sanitize ]]; then
  flags=(-O3 -DNDEBUG)
  binary=joint
  if [[ $phase == build-sanitize ]]; then
    flags=(-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer)
    binary=joint-sanitize
  fi
  g++ "${flags[@]}" -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw \
    -mtune=znver4 -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -I"$reference/spin/src/kernels" -I"$reference/spin/include" \
    "$root/joint.cpp" "$reference/build/CMakeFiles/locality_bch.dir/BchGfni.cpp.o" \
    "$reference/build/spin/libspin.a" -o "$root/$binary"
  sha256sum "$root/$binary"
  exit 0
fi
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/shared-gf16-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/joint"
if [[ $phase == check || $phase == sanitize ]]; then
  if [[ $phase == sanitize ]]; then binary="$root/joint-sanitize"; fi
  for seed in 1 17; do
    for updates in 2 3 4; do
      for mode in 9 10 11; do
        "$binary" 14 4 "$updates" "$mode" "$seed" 0 > "$logs/$seed-$updates-$mode.txt"
      done
    done
  done
  "$binary" 20 4 4 9 17 0 > "$logs/shared-full-size.txt"
  "$binary" 18 4 4 7 17 0 > "$logs/independent-control.txt"
  "$binary" 18 2 2 5 17 0 > "$logs/two-bit-control.txt"
  printf 'Passed shared GF/lane routes and unchanged controls, three input patterns each.\n'
  exit 0
fi
[[ $phase == confirm || $phase == profile ]] || exit 2
calls=101;profile=0
if [[ $phase == profile ]]; then calls=31;profile=1;fi
# Compare lane-only R2, shared GF R2/R4, and independent GF R4 in
# balanced serial order. Setup and reference checks are outside timing.
for seed in 1 17; do
  slot=0
  for spec in 2:11 2:9 4:9 4:7 4:7 4:9 2:9 2:11; do
    updates=${spec%:*};mode=${spec#*:}
    taskset -c 15 "$binary" 20 4 "$updates" "$mode" "$seed" "$calls" "$profile" > "$logs/$seed-$slot.csv"
    cat "$logs/$seed-$slot.csv"
    slot=$((slot+1))
  done
done
sha256sum "$binary"
