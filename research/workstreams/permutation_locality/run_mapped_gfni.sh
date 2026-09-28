#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/mapped-gfni"
for seed in 1 17; do
  taskset -c 15 "$root/build/locality" 14 4 1 block-gfni-check "$seed" 1 4
done
for mode in block-gfni-blend random-gfni random-gfni-repack random-gfni-prepared block-gfni-profile random-gfni-profile random-gfni-repack-profile random-gfni-prepared-profile; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 4 > "$root/measurements/mapped-gfni/$mode.csv"
  cat "$root/measurements/mapped-gfni/$mode.csv"
done
