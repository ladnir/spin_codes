#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?directory containing joint.cpp and headers}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
if [[ $phase == build || $phase == build-sanitize || $phase == build-mixer || $phase == build-mixer-sanitize ]]; then
  flags=(-O3 -DNDEBUG); binary=joint
  source=joint.cpp
  if [[ $phase == build-mixer ]]; then binary=pair_mixer; source=pair_mixer_perf.cpp; fi
  if [[ $phase == build-mixer-sanitize ]]; then
    flags=(-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer)
    binary=pair_mixer-sanitize; source=pair_mixer_perf.cpp
  fi
  if [[ $phase == build-sanitize ]]; then
    flags=(-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer)
    binary=joint-sanitize
  fi
  g++ "${flags[@]}" -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw \
    -mtune=znver4 -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
    -I"$reference/spin/src/kernels" -I"$reference/spin/include" \
    "$root/$source" "$reference/build/CMakeFiles/locality_bch.dir/BchGfni.cpp.o" \
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
logs=$(mktemp -d "$root/measurements/fused-bch-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/joint"
if [[ $phase == mixer-r1-screen ]]; then
  "$root/pair_mixer" 14 1 17 0 > "$logs/check.txt"
  taskset -c 15 "$root/pair_mixer" 20 1 17 31 1 > "$logs/screen.csv"
  cat "$logs/screen.csv"
  exit 0
fi
if [[ $phase == sanitize-mixer ]]; then
  for seed in 1 17; do
    for updates in 1 2 4; do
      for local_copy in 0 1; do
        "$root/pair_mixer-sanitize" 14 "$updates" "$seed" 0 0 "$local_copy" > "$logs/$seed-$updates-$local_copy.txt"
      done
    done
  done
  printf 'ASan/UBSan pair mixer checks pass. Linked prebuilt BCH/library objects are not instrumented.\n'
  exit 0
fi
if [[ $phase == check-mixer ]]; then
  for seed in 1 17; do
    for updates in 1 2 4; do
      "$root/pair_mixer" 14 "$updates" "$seed" 0 > "$logs/$seed-$updates.txt"
    done
  done
  "$root/pair_mixer" 20 2 17 0 > "$logs/mixer-full-size.txt"
  printf 'Passed scalar full encoder, local adjoint, sparse/dense replay and suffix.\n'
  exit 0
fi
if [[ $phase == mixer || $phase == profile-mixer ]]; then
  calls=101; profile=0
  if [[ $phase == profile-mixer ]]; then calls=31; profile=1; fi
  for seed in 1 17; do
    slot=0
    for spec in baseline:2 mixer:2 baseline:4 mixer:4 mixer:4 baseline:4 mixer:2 baseline:2; do
      kind=${spec%:*}; updates=${spec#*:}
      if [[ $kind == mixer ]]; then
        taskset -c 15 "$root/pair_mixer" 20 "$updates" "$seed" "$calls" "$profile" > "$logs/$seed-$slot-$kind-$updates.csv"
      else
        taskset -c 15 "$binary" 20 4 "$updates" 9 "$seed" "$calls" "$profile" > "$logs/$seed-$slot-$kind-$updates.csv"
      fi
      cat "$logs/$seed-$slot-$kind-$updates.csv"
      slot=$((slot+1))
    done
  done
  sha256sum "$root/pair_mixer" "$binary"
  exit 0
fi
if [[ $phase == check || $phase == sanitize ]]; then
  if [[ $phase == sanitize ]]; then binary="$root/joint-sanitize"; fi
  for seed in 1 17; do
    for updates in 2 4; do
      for mode in 7 16 17; do
        "$binary" 14 4 "$updates" "$mode" "$seed" 0 > "$logs/$seed-$updates-$mode.txt"
      done
    done
  done
  "$binary" 20 4 4 16 17 0 > "$logs/fused-full-size.txt"
  "$binary" 20 4 4 17 17 0 > "$logs/direct-full-size.txt"
  printf 'Passed full encoder/suffix, dense inner, adjoint, GF16, route and sparse/dense replay.\n'
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=101; profile=0; seeds=(1 17); modes=(7 16 17 17 16 7)
if [[ $phase == screen ]]; then calls=31; seeds=(17); modes=(7 16 17); fi
if [[ $phase == profile ]]; then calls=31; profile=1; fi
for seed in "${seeds[@]}"; do
  slot=0
  for mode in "${modes[@]}"; do
    taskset -c 15 "$binary" 20 4 4 "$mode" "$seed" "$calls" "$profile" > "$logs/$seed-$slot-$mode.csv"
    cat "$logs/$seed-$slot-$mode.csv"
    slot=$((slot+1))
  done
done
sha256sum "$binary"
