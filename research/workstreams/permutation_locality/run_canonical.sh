#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/canonical"
taskset -c 15 "$root/build/locality" 20 1 0 tiled 17 101 > "$root/measurements/canonical/control.csv"
cat "$root/measurements/canonical/control.csv"
taskset -c 15 "$root/build/locality" 20 16 1 wide-stream 17 101 > "$root/measurements/canonical/prior-wide.csv"
cat "$root/measurements/canonical/prior-wide.csv"
for columns in 2 4; do
  for mode in block-tiled block-cached block-stream block-profile; do
    taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 "$columns" > "$root/measurements/canonical/$mode-c$columns.csv"
    cat "$root/measurements/canonical/$mode-c$columns.csv"
  done
done
