#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/one-column"
for seed in 1 17; do
  for rep in 1 2 3; do
    modes=(baseline c2 c1)
    if [[ $rep == 2 ]]; then modes=(c1 c2 baseline); fi
    for selection in "${modes[@]}"; do
      mode=random-gfni-physical; group=4; shared=1
      case $selection in
        baseline) mode=tiled; group=1; shared=0; columns=1;;
        c2) columns=2;;
        c1) columns=1;;
      esac
      taskset -c 15 "$root/build/locality" 20 "$group" "$shared" "$mode" "$seed" 101 "$columns" > "$root/measurements/one-column/$seed-$rep-$selection.csv"
      cat "$root/measurements/one-column/$seed-$rep-$selection.csv"
    done
  done
done
sha256sum "$root/build/locality" "$root/build/spin/libspin.a"
