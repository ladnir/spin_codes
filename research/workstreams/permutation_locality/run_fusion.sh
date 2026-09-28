#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/fusion"
for mode in tiled sequential-control fused-1k fused-4k fused-16k; do
  taskset -c 15 "$root/build/locality" 20 1 0 "$mode" 17 51 > "$root/measurements/fusion/$mode.csv"
  cat "$root/measurements/fusion/$mode.csv"
done
