#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
for seed in 1 17; do
  "$root/build-sanitize/locality" 14 4 1 r2-physical "$seed" 1 1
done
mkdir -p "$root/measurements/two-updates"
for seed in 1 17; do
  for rep in 1 2 3; do
    modes=(random-gfni-physical r2-physical)
    if [[ $rep == 2 ]]; then modes=(r2-physical random-gfni-physical); fi
    for mode in "${modes[@]}"; do
      taskset -c 15 "$root/build/locality" 20 4 1 "$mode" "$seed" 101 1 > "$root/measurements/two-updates/$seed-$rep-$mode.csv"
      cat "$root/measurements/two-updates/$seed-$rep-$mode.csv"
    done
  done
done
sha256sum "$root/build/locality" "$root/build/spin/libspin.a"
