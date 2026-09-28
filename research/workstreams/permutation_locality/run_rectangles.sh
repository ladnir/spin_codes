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
mkdir -p "$root/measurements/rectangles"
for g in 1 4; do
  if ((g==1)); then shared=0; else shared=1; fi
  for c in 1 2 4 8 16 32; do
    # Kernel-rank obstructions reject g*c >=64 when g*c divides 128.
    if ((g*c>=64)); then continue; fi
    "$root/build/locality" 20 "$g" "$shared" rank 17 1 "$c" > "$root/measurements/rectangles/rank-$g-$c.csv"
    cat "$root/measurements/rectangles/rank-$g-$c.csv"
    modes='tiled direct'
    if ((g==4)); then modes='tiled direct stream'; fi
    for mode in $modes; do
      taskset -c "$cpu" "$root/build/locality" 20 "$g" "$shared" "$mode" 17 31 "$c" > "$root/measurements/rectangles/$g-$c-$mode.csv"
      cat "$root/measurements/rectangles/$g-$c-$mode.csv"
    done
  done
done
