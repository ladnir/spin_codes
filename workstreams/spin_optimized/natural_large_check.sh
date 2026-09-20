#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/natural_lengths
mkdir -p "$result"
# Requires the release build. Cross the packed-index boundary by one natural unit.
for target in spin_lengths_test spin_forward_lengths_test spin_wide_lengths_test; do
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    "build-natural-release/$target" 8404992 > "$result/large-$target.log" 2>&1
done
