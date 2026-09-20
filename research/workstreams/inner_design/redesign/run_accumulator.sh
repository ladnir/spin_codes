#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated build root}
cd "$root"
for repeat in 1 2 3; do
    variants=(control 2 3)
    if [[ $repeat == 2 ]]; then variants=(3 2 control); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        if [[ $variant == control ]]; then
            taskset -c 15 build/basis_asymmetric_greedy3_2_sparse_bench 101 0 0 20 1 2 1 | tee "redesign-accumulator-control-repeat${repeat}.jsonl"
        else
            taskset -c 15 redesign-build/accumulator_bench "$variant" | tee "redesign-accumulator-${variant}-repeat${repeat}.jsonl"
        fi
    done
done
sha256sum build/basis_asymmetric_greedy3_2_sparse_bench redesign-build/accumulator_bench | tee redesign-accumulator-binaries.sha256
