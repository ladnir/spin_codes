#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/natural_lengths
mkdir -p "$result"
# All builds and correctness checks must finish before this serial timing pass.
python3 workstreams/spin_optimized/isa_check.py build-natural-release > "$result/isa.log"
for repeat in 1 2 3; do
  sizes=(16 18 20)
  [[ $repeat != 2 ]] || sizes=(20 18 16)
  for m in "${sizes[@]}"; do
    build-natural-release/spin_single_benchmark "$m" forward 31 256 > "$result/forward-m$m-r$repeat.json"
    build-natural-release/spin_single_transpose_benchmark "$m" transpose 31 0 > "$result/transpose-m$m-r$repeat.json"
    cfg=12819;tile=256
    [[ $m != 16 ]] || cfg=6412r2
    [[ $m != 20 ]] || tile=512
    build-natural-release/spin_wide_benchmark "$m" 3 "$cfg" "$tile" 0 > "$result/wide-m$m-r$repeat.csv"
  done
done
