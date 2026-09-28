#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
cpu=${2:-15}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/mapped"
taskset -c "$cpu" "$root/build/locality" 20 1 0 tiled 17 51 > "$root/measurements/mapped/control.csv"
cat "$root/measurements/mapped/control.csv"
for c in 1 2 4 8; do
  for mode in tiled mapped mapped-stream; do
    taskset -c "$cpu" "$root/build/locality" 20 4 1 "$mode" 17 51 "$c" > "$root/measurements/mapped/$c-$mode.csv"
    cat "$root/measurements/mapped/$c-$mode.csv"
  done
done
