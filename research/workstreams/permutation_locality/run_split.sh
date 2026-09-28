#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/split"
taskset -c 15 "$root/build/locality" 20 1 0 tiled 17 51 > "$root/measurements/split/control.csv"
cat "$root/measurements/split/control.csv"
taskset -c 15 "$root/build/locality" 20 16 1 wide-stream 17 51 > "$root/measurements/split/wide16.csv"
cat "$root/measurements/split/wide16.csv"
for g in 32 64 128; do
  for mode in split-tiled split-stream; do
    taskset -c 15 "$root/build/locality" 20 "$g" 1 "$mode" 17 51 > "$root/measurements/split/$mode-g$g.csv"
    cat "$root/measurements/split/$mode-g$g.csv"
  done
done
