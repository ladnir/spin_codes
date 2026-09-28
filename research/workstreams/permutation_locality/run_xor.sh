#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/xor"
for mode in tiled wide-stream wide-profile xor-tiled xor-wide-stream xor-wide-profile; do
  if [[ $mode == tiled ]]; then g=1; shared=0; else g=16; shared=1; fi
  taskset -c 15 "$root/build/locality" 20 "$g" "$shared" "$mode" 17 101 > "$root/measurements/xor/$mode.csv"
  cat "$root/measurements/xor/$mode.csv"
done
