#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/state-gather"
taskset -c 15 "$root/build/locality" 20 1 0 tiled 17 51 > "$root/measurements/state-gather/control.csv"
cat "$root/measurements/state-gather/control.csv"
for c in 1 4 8; do
  for mode in tiled state-gather; do
    taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 51 "$c" > "$root/measurements/state-gather/$mode-c$c.csv"
    cat "$root/measurements/state-gather/$mode-c$c.csv"
  done
done
