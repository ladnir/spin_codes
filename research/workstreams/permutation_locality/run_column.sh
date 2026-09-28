#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/column"
for pass in 1 2 3; do
  modes='tiled column'
  if ((pass==2)); then modes='column tiled'; fi
  for mode in $modes; do
    if [[ $mode == tiled ]]; then g=1; shared=0; else g=4; shared=1; fi
    taskset -c 15 "$root/build/locality" 20 "$g" "$shared" "$mode" 17 101 > "$root/measurements/column/$mode-$pass.csv"
    cat "$root/measurements/column/$mode-$pass.csv"
  done
done
