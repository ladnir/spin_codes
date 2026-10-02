#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing joint and joint-sanitize}")
phase=${2:-confirm}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/gf16-updates-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/joint"
if [[ $phase == check || $phase == sanitize ]]; then
  if [[ $phase == sanitize ]]; then binary="$root/joint-sanitize"; fi
  for seed in 1 17; do
    for updates in 2 3 4; do
      for mode in 7 8; do
        "$binary" 14 4 "$updates" "$mode" "$seed" 0 > "$logs/$seed-$updates-$mode.txt"
      done
    done
  done
  "$binary" 18 4 4 7 17 0 > "$logs/r4-larger.txt"
  "$binary" 18 4 2 1 17 0 > "$logs/lane-control.txt"
  "$binary" 18 2 2 5 17 0 > "$logs/two-bit-control.txt"
  printf 'Passed GF16 update-count checks, larger R4 case, and unchanged controls.\n'
  exit 0
fi
[[ $phase == confirm ]] || exit 2
# Balanced serial order. Setup, allocation, and correctness checks are
# outside the 101 timed precomputed encoder calls. No phase timers.
order=(2 3 4 4 3 2)
for seed in 1 17; do
  for slot in "${!order[@]}"; do
    updates=${order[$slot]}
    taskset -c 15 "$binary" 20 4 "$updates" 7 "$seed" 101 0 > "$logs/$seed-$slot-r$updates.csv"
    cat "$logs/$seed-$slot-r$updates.csv"
  done
done
sha256sum "$binary"
