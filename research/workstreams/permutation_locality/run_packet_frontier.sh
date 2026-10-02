#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
phase=${2:-screen}
tile=${3:-16}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 75
mkdir -p "$root/measurements"
logs=$(mktemp -d "$root/measurements/packet-frontier-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
binary="$root/build/joint"
if [[ $phase == sanitize ]]; then binary="$root/build-sanitize/joint"; fi
if [[ $phase == check || $phase == sanitize ]]; then
  for exponent in 14 18 20; do
    for seed in 1 17; do
      for rows in 4 8 16 32 64 128; do
        "$binary" "$exponent" 2 2 6 "$seed" 0 0 "$rows" > "$logs/$exponent-$seed-$rows.txt"
      done
    done
  done
  printf 'Passed 36 gather cases, with three input patterns and layout rejection tests.\n'
  exit 0
fi
[[ $phase == screen || $phase == confirm ]] || exit 2
for seed in 1 17; do
  if [[ $phase == screen ]]; then
    for rows in 4 8 16 32 64 128 256; do
      taskset -c 15 "$binary" 20 2 2 6 "$seed" 31 1 "$rows" > "$logs/$seed-gather-$rows.csv"
      cat "$logs/$seed-gather-$rows.csv"
    done
    taskset -c 15 "$binary" 20 2 2 5 "$seed" 31 1 16 > "$logs/$seed-scatter.csv"
    cat "$logs/$seed-scatter.csv"
  else
    for rep in 1 2 3; do
      modes=(2:2:5 2:2:6 4:2:1 4:3:1)
      if ((rep==2)); then modes=(4:3:1 4:2:1 2:2:6 2:2:5); fi
      for mode in "${modes[@]}"; do
        IFS=: read -r packet updates method <<< "$mode"
        rows=16; if [[ $method == 6 ]]; then rows=$tile; fi
        taskset -c 15 "$binary" 20 "$packet" "$updates" "$method" "$seed" 101 0 "$rows" > "$logs/$seed-$rep-$mode.csv"
        cat "$logs/$seed-$rep-$mode.csv"
      done
      expected=$(awk -F, '/^20,2,/ {print $NF}' "$logs/$seed-$rep-2:2:5.csv")
      actual=$(awk -F, '/^20,2,/ {print $NF}' "$logs/$seed-$rep-2:2:6.csv")
      [[ -n $expected && $expected == "$actual" ]] || exit 1
    done
  fi
done
sha256sum "$binary" "$root/build/spin/libspin.a"
