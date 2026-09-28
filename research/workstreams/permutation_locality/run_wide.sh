#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/wide"
taskset -c 15 "$root/build/locality" 20 1 0 tiled 17 51 > "$root/measurements/wide/control.csv"
cat "$root/measurements/wide/control.csv"
for g in 8 16 32; do
  for mode in tiled wide wide-stream; do
    taskset -c 15 "$root/build/locality" 20 "$g" 1 "$mode" 17 51 > "$root/measurements/wide/$mode-g$g.csv"
    cat "$root/measurements/wide/$mode-g$g.csv"
  done
done
