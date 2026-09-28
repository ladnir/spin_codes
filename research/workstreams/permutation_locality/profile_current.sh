#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/current-profile"
for mode in random-gfni-copy-mapped-profile block-gfni-profile; do
  taskset -c 15 "$root/build/locality" 20 4 1 "$mode" 17 101 4
done
perf record -q -F 499 -o "$root/measurements/current-profile/perf.data" -- \
  taskset -c 15 "$root/build/locality" 20 4 1 random-gfni-copy-mapped 17 1001 4
perf report --stdio --no-children --sort=symbol -i "$root/measurements/current-profile/perf.data" \
  > "$root/measurements/current-profile/report.txt"
head -n 65 "$root/measurements/current-profile/report.txt"
