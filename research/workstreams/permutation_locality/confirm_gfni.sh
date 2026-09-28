#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/gfni-confirm"
for seed in 1 17; do
  for rep in 1 2 3; do
    modes=(tiled block-restrict block-gfni-ternary block-gfni-blend)
    if [[ $rep == 2 ]]; then modes=(block-gfni-blend block-gfni-ternary block-restrict tiled); fi
    for mode in "${modes[@]}"; do
      group=4; shared=1; columns=4
      if [[ $mode == tiled ]]; then group=1; shared=0; columns=1; fi
      taskset -c 15 "$root/build/locality" 20 "$group" "$shared" "$mode" "$seed" 101 "$columns" > "$root/measurements/gfni-confirm/$seed-$rep-$mode.csv"
      cat "$root/measurements/gfni-confirm/$seed-$rep-$mode.csv"
    done
  done
done
