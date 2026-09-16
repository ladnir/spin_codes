#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated build root}
cd "$root"
ctest --test-dir build -R 'basis_redesign_gather.*_test' -V -j1 | tee redesign-routing-correctness.log
for repeat in 1 2 3; do
    variants=(asymmetric_greedy3_2_sparse redesign_gather_pf0 redesign_gather_pf64)
    if [[ $repeat == 2 ]]; then variants=(redesign_gather_pf64 redesign_gather_pf0 asymmetric_greedy3_2_sparse); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        taskset -c 15 "build/basis_${variant}_bench" 101 0 0 20 1 2 1 | tee "redesign-routing-${variant}-repeat${repeat}.jsonl"
    done
done
sha256sum build/basis_asymmetric_greedy3_2_sparse_bench build/basis_redesign_gather_pf{0,64}_{bench,test} | tee redesign-routing-binaries.sha256
