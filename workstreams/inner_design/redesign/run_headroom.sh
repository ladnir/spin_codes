#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated build root}
cd "$root"
ctest --test-dir build -R 'basis_redesign_(no_inner|no_inner_unrolled|outer_only)_test' -V -j1 | tee redesign-headroom-correctness.log
for repeat in 1 2 3; do
    variants=(asymmetric_greedy3_2_sparse redesign_no_inner redesign_no_inner_unrolled redesign_outer_only)
    if [[ $repeat == 2 ]]; then variants=(redesign_outer_only redesign_no_inner_unrolled redesign_no_inner asymmetric_greedy3_2_sparse); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        taskset -c 15 "build/basis_${variant}_bench" 101 0 0 20 1 2 1 | tee "redesign-headroom-${variant}-repeat${repeat}.jsonl"
    done
done
sha256sum build/basis_asymmetric_greedy3_2_sparse_bench build/basis_{redesign_no_inner,redesign_no_inner_unrolled,redesign_outer_only}_{bench,test} | tee redesign-headroom-binaries.sha256
