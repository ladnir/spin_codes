#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/bundled"
for mode in plain-bch repack-bch; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 51 8 > "$root/measurements/bundled/$mode.csv"
  cat "$root/measurements/bundled/$mode.csv"
done
taskset -c 15 "$root/build/locality" 20 1 0 tiled 17 51 > "$root/measurements/bundled/control.csv"
cat "$root/measurements/bundled/control.csv"
for c in 1 2 4 8; do
  taskset -c 15 "$root/build/locality" 20 4 1 bundled 17 51 "$c" > "$root/measurements/bundled/$c.csv"
  cat "$root/measurements/bundled/$c.csv"
done
