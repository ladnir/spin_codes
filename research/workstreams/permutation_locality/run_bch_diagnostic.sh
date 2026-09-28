#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/bch-diagnostic"
for mode in plain-bch mapped-bch; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 8 > "$root/measurements/bch-diagnostic/$mode.csv"
  cat "$root/measurements/bch-diagnostic/$mode.csv"
done
