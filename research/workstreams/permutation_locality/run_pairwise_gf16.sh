#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?directory containing joint.cpp and headers}")
phase=${2:-check}
reference=${3:-/tmp/spin-joint-9n57TT}
if [[ $phase == build || $phase == build-sanitize ]]; then
  exec bash "$root/run_shared_gf16.sh" "$root" "$phase" "$reference"
fi
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/pairwise-gf16-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/joint"
if [[ $phase == check || $phase == sanitize || $phase == direct-check || $phase == direct-sanitize ]]; then
  if [[ $phase == sanitize || $phase == direct-sanitize ]]; then binary="$root/joint-sanitize"; fi
  modes=(12 13)
  if [[ $phase == direct-check || $phase == direct-sanitize ]]; then modes=(12 13 14 15);fi
  for seed in 1 17; do
    for updates in 2 3 4; do
      for mode in "${modes[@]}"; do
        "$binary" 14 4 "$updates" "$mode" "$seed" 0 > "$logs/$seed-$updates-$mode.txt"
      done
    done
  done
  "$binary" 20 4 4 12 17 0 > "$logs/pairwise-full-size.txt"
  if [[ $phase == direct-check || $phase == direct-sanitize ]]; then
    "$binary" 20 4 4 14 17 0 > "$logs/direct-full-size.txt"
  fi
  "$binary" 18 4 4 7 17 0 > "$logs/independent-control.txt"
  "$binary" 18 4 4 9 17 0 > "$logs/shared-control.txt"
  "$binary" 18 2 2 5 17 0 > "$logs/two-bit-control.txt"
  printf 'Passed pairwise routes and unchanged controls, three input patterns each.\n'
  exit 0
fi
[[ $phase == confirm || $phase == profile || $phase == direct-confirm || $phase == direct-profile ]] || exit 2
calls=101;profile=0
if [[ $phase == profile || $phase == direct-profile ]]; then calls=31;profile=1;fi
modes=(7 12 13 9 9 13 12 7)
if [[ $phase == direct-confirm || $phase == direct-profile ]]; then modes=(7 12 14 9 9 14 12 7);fi
# All are GF16, four updates, same full-line router. Balanced serial order.
for seed in 1 17; do
  slot=0
  for mode in "${modes[@]}"; do
    taskset -c 15 "$binary" 20 4 4 "$mode" "$seed" "$calls" "$profile" > "$logs/$seed-$slot.csv"
    cat "$logs/$seed-$slot.csv"
    slot=$((slot+1))
  done
done
sha256sum "$binary" "$root/joint.cpp"
