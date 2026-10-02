#!/usr/bin/env bash
# Serial R2/R3/R4 follow-on; retain the original R2 benchmark driver/script.
set -euo pipefail
root=$(realpath "${1:?fresh build directory}")
phase=${2:-check}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -w 0 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -w 0 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -w 0 7 || exit 75
if [[ $phase == build ]]; then
  bash "$root/packed_driver_run.sh" "$root" build
  exit 0
fi
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/updates-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
if [[ $phase == check ]]; then
  for updates in 2 3 4; do
    for seed in 1 17; do
      for exponent in 14 20; do
        taskset -c 15 "$root/packed_driver" "$exponent" 2 "$seed" 0 0 "$updates" > "$logs/check-$exponent-$updates-$seed.txt"
        cat "$logs/check-$exponent-$updates-$seed.txt"
      done
    done
  done
  exit 0
fi
[[ $phase == confirm || $phase == profile ]] || exit 2
calls=101; profile=0; repeats=2
if [[ $phase == profile ]]; then calls=31; profile=1; repeats=1; fi
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    if ((rep%2)); then order=(4 2); else order=(2 4); fi
    for updates in "${order[@]}"; do
      taskset -c 15 "$root/packed_driver" 20 2 "$seed" "$calls" "$profile" "$updates" > "$logs/$seed-$rep-r$updates.csv"
      cat "$logs/$seed-$rep-r$updates.csv"
    done
  done
done
