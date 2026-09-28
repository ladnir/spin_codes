#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
for seed in 1 17; do
  taskset -c 15 "$root/build/locality" 14 4 1 block-gfni-check "$seed" 1 4
done
mkdir -p "$root/measurements/copy-packed"
for mode in random-gfni-copy-mapped random-gfni-manual-copy random-gfni-copy-mapped-prefetch random-gfni-copy-packed random-gfni-manual-packed random-gfni-vector-packed random-gfni-vector-mapped; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 4 > "$root/measurements/copy-packed/$mode.csv"
  cat "$root/measurements/copy-packed/$mode.csv"
done
