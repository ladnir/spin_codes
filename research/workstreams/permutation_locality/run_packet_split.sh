#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
phase=${2:-screen}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
logs=$(mktemp -d "$root/measurements/packet-split-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == sanitize ]]; then
  for seed in 1 17; do
    for mode in r2-repack r2-fused r2-scatter r2-repack-copy r2-repack-prefetch r2-fused-copy; do
      "$root/build-sanitize/locality" 14 4 0 "$mode" "$seed" 1 1
    done
  done
  exit 0
fi
calls=31
repetitions=1
if [[ $phase == confirm ]]; then calls=101; repetitions=3; fi
for seed in 1 17; do
  for ((rep=1; rep<=repetitions; ++rep)); do
    modes=(r2-physical r2-repack r2-fused r2-scatter)
    if (( rep==2 )); then modes=(r2-scatter r2-fused r2-repack r2-physical); fi
    if [[ $phase == warm ]]; then modes=(r2-physical r2-repack-copy r2-repack-prefetch r2-fused); fi
    if [[ $phase == final || $phase == confirm ]]; then
      modes=(r2-physical r2-repack-copy r2-fused-copy)
      if (( rep==2 )); then modes=(r2-fused-copy r2-repack-copy r2-physical); fi
    fi
    for mode in "${modes[@]}"; do
      shared=0
      if [[ $mode == r2-physical ]]; then shared=1; fi
      taskset -c 15 "$root/build/locality" 20 4 "$shared" "$mode-profile" "$seed" "$calls" 1 > "$logs/$seed-$rep-$mode.csv"
      cat "$logs/$seed-$rep-$mode.csv"
    done
  done
done
if [[ $phase == screen ]]; then
  # Same-map controls isolate the cost of splitting/repacking locally.
  for mode in r2-repack r2-fused; do
    taskset -c 15 "$root/build/locality" 20 4 1 "$mode-profile" 1 31 1 > "$logs/shared-$mode.csv"
    cat "$logs/shared-$mode.csv"
  done
fi
sha256sum "$root/build/locality" "$root/build/spin/libspin.a"
