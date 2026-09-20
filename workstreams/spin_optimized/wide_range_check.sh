#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
build=${2:-build-wide-schedule-dfs}
result=workstreams/spin_optimized/measurements/wide_range
mkdir -p "$result"
# Use the validated DFS build. Do not build anything while timing.
grep -qx 'SPIN_WIDE_SCHEDULE:STRING=dfs' "$build/CMakeCache.txt"
grep -qx 'SPIN_FORWARD_SCHEDULE:STRING=dfs' "$build/CMakeCache.txt"
grep -qx 'SPIN_FORWARD_DIRECT:BOOL=ON' "$build/CMakeCache.txt"
for repeat in 1 2 3; do
  tiles=(256 512)
  [[ $repeat != 2 ]] || tiles=(512 256)
  for m in 17 19; do
    for tile in "${tiles[@]}"; do
      for lanes in 1 2 4; do
        "$build/spin_wide_benchmark" "$m" 5 12819 "$tile" "$lanes" \
          > "$result/dfs-m$m-t$tile-w$lanes-r$repeat.csv"
      done
    done
  done
done
