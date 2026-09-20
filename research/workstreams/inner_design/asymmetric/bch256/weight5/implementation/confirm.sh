#!/usr/bin/env bash
set -euo pipefail
cd "${1:?pass isolated remote root}"
for repeat in 1 2 3; do
    variants=(baseline baseline_tuned sparse sparse_pages sparse_tuned masked_tuned shared_tuned grouped_tuned sparse_pf128)
    if [[ $repeat == 2 ]]; then variants=(sparse_pf128 grouped_tuned shared_tuned masked_tuned sparse_tuned sparse_pages sparse baseline_tuned baseline); fi
    if [[ $repeat == 3 ]]; then variants=(masked_tuned baseline sparse_tuned baseline_tuned grouped_tuned sparse shared_tuned sparse_pages sparse_pf128); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        taskset -c 15 "build/${variant}_bench" 101 2048 0 20 0 2 1 | tee "measurements/confirm-${variant}-${repeat}.jsonl"
    done
done
