#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build and build-sanitize}")
phase=${2:-screen}
tile=${3:-16}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/two-bit-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/build/joint"
if [[ $phase == sanitize ]]; then binary="$root/build-sanitize/joint"; fi
if [[ $phase == check || $phase == sanitize ]]; then
  for seed in 1 17; do
    for updates in 1 2 3; do
      for mode in 0 1 2 3 4 5; do
        for rows in 4 64 2048; do
          "$binary" 14 2 "$updates" "$mode" "$seed" 0 0 "$rows" > "$logs/$seed-$updates-$mode-$rows.txt"
        done
      done
    done
  done
  for rows in 64 256 2048; do
    "$binary" 18 2 2 3 17 0 0 "$rows" > "$logs/k18-$rows.txt"
  done
  for exponent in 16 18 20; do
    "$binary" "$exponent" 2 2 5 17 0 0 16 > "$logs/padded-k$exponent.txt"
  done
  printf 'Passed 108 small cases plus 6 larger cases through K=2^20 (three input patterns each).\n'
  exit 0
fi
[[ $phase == screen || $phase == padded || $phase == confirm ]] || exit 2
for seed in 1 17; do
  if [[ $phase == padded ]]; then
    for rows in 8 16 32 64 128 256 512 1024 2048; do
      taskset -c 15 "$binary" 20 2 2 5 "$seed" 31 1 "$rows" > "$logs/$seed-padded-$rows.csv"
      cat "$logs/$seed-padded-$rows.csv"
    done
  elif [[ $phase == screen ]]; then
    for rows in 64 128 256 512 1024 2048 4096; do
      taskset -c 15 "$binary" 20 2 2 3 "$seed" 31 1 "$rows" > "$logs/$seed-tile-$rows.csv"
      cat "$logs/$seed-tile-$rows.csv"
    done
    for mode in 0 2 4; do
      taskset -c 15 "$binary" 20 2 2 "$mode" "$seed" 31 1 "$tile" > "$logs/$seed-control-$mode.csv"
      cat "$logs/$seed-control-$mode.csv"
    done
  else
    for rep in 1 2 3; do
      modes=(baseline 0 2 3 4 5 one-update)
      if ((rep==2)); then modes=(one-update 5 4 3 2 0 baseline); fi
      for mode in "${modes[@]}"; do
        if [[ $mode == baseline ]]; then
          taskset -c 15 "$root/build/locality" 20 1 0 tiled "$seed" 101 > "$logs/$seed-$rep-$mode.csv"
        elif [[ $mode == one-update ]]; then
          taskset -c 15 "$binary" 20 2 1 5 "$seed" 101 0 "$tile" > "$logs/$seed-$rep-$mode.csv"
        else
          taskset -c 15 "$binary" 20 2 2 "$mode" "$seed" 101 0 "$tile" > "$logs/$seed-$rep-$mode.csv"
        fi
        cat "$logs/$seed-$rep-$mode.csv"
      done
      expected=$(awk -F, '/^20,2,128,19,2,/ {print $NF}' "$logs/$seed-$rep-0.csv")
      for mode in 2 3 4 5; do
        actual=$(awk -F, '/^20,2,128,19,2,/ {print $NF}' "$logs/$seed-$rep-$mode.csv")
        [[ -n $expected && $actual == "$expected" ]] || { printf 'Repeated-encoding checksum mismatch\n' >&2; exit 1; }
      done
    done
  fi
done
sha256sum "$binary" "$root/build/spin/libspin.a"
