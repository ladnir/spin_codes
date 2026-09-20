#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/forward_schedule
cmake --build build-forward-schedule-dfs --target spin_single_transpose_benchmark -j2 > "$result/transpose-build.log" 2>&1
for repeat in 1 2 3; do
  directions=(forward transpose inplace)
  [[ $repeat != 2 ]] || directions=(inplace transpose forward)
  for m in 16 18 20; do
    for direction in "${directions[@]}"; do
      exe=spin_single_benchmark; tile=256
      [[ $direction == forward ]] || { exe=spin_single_transpose_benchmark; tile=0; }
      "build-forward-schedule-dfs/$exe" "$m" "$direction" 101 "$tile" > "$result/confirm-m$m-$direction-r$repeat.json"
    done
  done
done
