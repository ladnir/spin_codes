#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/physical-columns"
for mode in random-gfni-vector-mapped random-gfni-physical random-gfni-physical-copy random-gfni-physical-prefetch random-gfni-vector-current-prefetch; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 4 > "$root/measurements/physical-columns/$mode.csv"
  cat "$root/measurements/physical-columns/$mode.csv"
done
for mode in random-gfni-physical-profile random-gfni-physical-copy-profile random-gfni-vector-current-prefetch-profile; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 4
done
