#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated build root}
cd "$root"
for repeat in 1 2 3; do
    variants=(control bucketed)
    if [[ $repeat == 2 ]]; then variants=(bucketed control); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        if [[ $variant == control ]]; then
            taskset -c 15 build/basis_asymmetric_greedy3_2_sparse_bench 101 0 0 20 1 2 1 | tee "redesign-accumulator-bucket-control-repeat${repeat}.jsonl"
        else
            taskset -c 15 redesign-build/accumulator_bench 2 1 | tee "redesign-accumulator-bucket-2-repeat${repeat}.jsonl"
        fi
    done
done
sha256sum build/basis_asymmetric_greedy3_2_sparse_bench redesign-build/accumulator_bench | tee redesign-accumulator-bucket-binaries.sha256
