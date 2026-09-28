#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/floor"
for pass in 1 2 3; do
  modes='tiled sequential-control'
  if ((pass==2)); then modes='sequential-control tiled'; fi
  for mode in $modes; do
    taskset -c 15 "$root/build/locality" 20 1 0 "$mode" 17 101 > "$root/measurements/floor/$mode-$pass.csv"
    cat "$root/measurements/floor/$mode-$pass.csv"
  done
done
