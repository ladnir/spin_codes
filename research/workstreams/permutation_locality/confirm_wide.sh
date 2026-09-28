#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/wide-confirm"
for seed in 1 17; do
  for pass in 1 2 3; do
    modes='tiled wide-stream'
    if ((pass==2)); then modes='wide-stream tiled'; fi
    for mode in $modes; do
      if [[ $mode == tiled ]]; then g=1; shared=0; else g=16; shared=1; fi
      taskset -c 15 "$root/build/locality" 20 "$g" "$shared" "$mode" "$seed" 101 > "$root/measurements/wide-confirm/$mode-s$seed-p$pass.csv"
      cat "$root/measurements/wide-confirm/$mode-s$seed-p$pass.csv"
    done
  done
done
