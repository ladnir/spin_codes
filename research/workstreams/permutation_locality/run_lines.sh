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
mkdir -p "$root/measurements"
for pass in 1 2 3; do
  modes='control tiled direct line stream'
  if ((pass==2)); then modes='stream line direct tiled control'; fi
  for mode in $modes; do
    if [[ $mode == control ]]; then g=1; shared=0; kernel=tiled;
    else g=4; shared=1; kernel=$mode; fi
    taskset -c "$cpu" "$root/build/locality" 20 "$g" "$shared" "$kernel" 17 101 > "$root/measurements/lines-$mode-$pass.csv"
    cat "$root/measurements/lines-$mode-$pass.csv"
  done
done
