#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing joint}")
phase=${2:-check}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/gf16-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/joint"
if [[ $phase == check || $phase == sanitize ]]; then
  if [[ $phase == sanitize ]]; then binary="$root/joint-sanitize"; fi
  for seed in 1 17; do
    for updates in 2 3; do
      for mode in 7 8; do
        for exponent in 14 18 20; do
          "$binary" "$exponent" 4 "$updates" "$mode" "$seed" 0 > "$logs/$exponent-$seed-$updates-$mode.txt"
        done
      done
    done
  done
  "$binary" 18 4 2 1 17 0 > "$logs/lane-control.txt"
  "$binary" 18 2 2 5 17 0 > "$logs/two-bit-control.txt"
  printf 'Passed 24 GF16 cases and two unchanged controls; three input patterns each.\n'
  exit 0
fi
[[ $phase == screen || $phase == confirm ]] || exit 2
calls=31;reps=1;profile=1
if [[ $phase == confirm ]]; then calls=101;reps=3;profile=0;fi
for seed in 1 17; do
  for ((rep=1;rep<=reps;rep++)); do
    modes=(1 7 8 5)
    if ((rep==2));then modes=(5 8 7 1);fi
    for mode in "${modes[@]}"; do
      packet=4;if [[ $mode == 5 ]];then packet=2;fi
      taskset -c 15 "$binary" 20 "$packet" 2 "$mode" "$seed" "$calls" "$profile" > "$logs/$seed-$rep-$mode.csv"
      cat "$logs/$seed-$rep-$mode.csv"
    done
    expected=$(awk -F, '/^20,4,128,19,2,/ {print $NF}' "$logs/$seed-$rep-7.csv")
    actual=$(awk -F, '/^20,4,128,19,2,/ {print $NF}' "$logs/$seed-$rep-8.csv")
    [[ -n $expected && $actual == "$expected" ]] || { printf 'GF16 mode checksum mismatch\n' >&2;exit 1; }
  done
done
sha256sum "$binary"
