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
for g in 1 2 4 8 16; do
  for shared in 0 1; do
    if ((g==1 && shared==1)); then continue; fi
    "$root/build/locality" 20 "$g" "$shared" rank 1 > "$root/measurements/rank-$g-$shared.csv"
    cat "$root/measurements/rank-$g-$shared.csv"
    for mode in tiled direct; do
      taskset -c "$cpu" "$root/build/locality" 20 "$g" "$shared" "$mode" 1 31 > "$root/measurements/screen-$g-$shared-$mode.csv"
      cat "$root/measurements/screen-$g-$shared-$mode.csv"
    done
  done
done
