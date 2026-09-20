#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/wide_tuning
# Run only after wide_tune.sh has completed. The executable holds both locks.
for repeat in 1 2 3; do
  tiles=(2048 256 512)
  if [[ $repeat == 2 ]]; then tiles=(512 256 2048); fi
  if [[ $repeat == 3 ]]; then tiles=(256 2048 512); fi
  for tile in "${tiles[@]}"; do
    build-wide/spin_wide_benchmark 20 7 12819 "$tile" > "$result/confirm-t$tile-r$repeat.csv"
  done
done
