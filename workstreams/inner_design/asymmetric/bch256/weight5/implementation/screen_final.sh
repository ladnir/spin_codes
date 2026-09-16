#!/usr/bin/env bash
set -euo pipefail
cd "${1:?pass isolated remote root}"
for variant in baseline baseline_tuned sparse sparse_pages sparse_tuned masked_tuned shared_tuned grouped_tuned sparse_pf128; do
    flock -w 45 /tmp/bare-spin-benchmark.lock true
    taskset -c 15 "build/${variant}_bench" 51 2048 0 20 0 2 1 | tee "measurements/screen-final-${variant}.jsonl"
done
