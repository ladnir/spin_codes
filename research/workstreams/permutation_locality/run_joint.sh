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
if [[ $phase == check || $phase == sanitize ]]; then
  build=build
  if [[ $phase == sanitize ]]; then build=build-sanitize; fi
  for seed in 1 17; do
    for packet in 2 4; do
      for updates in 2 3; do
        for stream in 0 1; do
          "$root/$build/joint" 14 "$packet" "$updates" "$stream" "$seed" 0
        done
      done
    done
  done
  exit 0
fi
[[ $phase == screen || $phase == confirm ]] || exit 2
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/joint-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
calls=31;repetitions=1
if [[ $phase == confirm ]]; then calls=101;repetitions=3;fi
for seed in 1 17; do
  for ((rep=1;rep<=repetitions;++rep)); do
    modes=(2:2:0 2:3:0 4:2:1 4:3:1 2:2:1 2:3:1 4:2:0 4:3:0)
    if [[ $phase == confirm ]];then modes=(2:2:0 2:3:0 4:2:1 4:3:1);fi
    if ((rep==2));then modes=(4:3:0 4:2:0 2:3:1 2:2:1 4:3:1 4:2:1 2:3:0 2:2:0);fi
    if [[ $phase == confirm ]] && ((rep==2));then modes=(4:3:1 4:2:1 2:3:0 2:2:0);fi
    for mode in "${modes[@]}"; do
      IFS=: read -r packet updates stream <<< "$mode"
      taskset -c 15 "$root/build/joint" 20 "$packet" "$updates" "$stream" "$seed" "$calls" 1 > "$logs/$seed-$rep-$mode.csv"
      cat "$logs/$seed-$rep-$mode.csv"
    done
  done
done
sha256sum "$root/build/joint" "$root/build/spin/libspin.a"
