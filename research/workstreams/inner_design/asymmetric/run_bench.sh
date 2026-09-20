#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass isolated build root}
cd "$root"
for run in 1 2 3; do
    variants=(balanced mixed2_3_sparse mixed2_3_masked greedy3_2_sparse greedy3_2_masked)
    if [[ $run == 2 ]]; then variants=(greedy3_2_masked greedy3_2_sparse mixed2_3_masked mixed2_3_sparse balanced); fi
    if [[ $run == 3 ]]; then variants=(balanced greedy3_2_masked greedy3_2_sparse mixed2_3_masked mixed2_3_sparse); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        binary="build/basis_asymmetric_${variant}_bench"
        if [[ $variant == balanced ]]; then binary=build/basis_balanced_masked_bench; fi
        # The executable acquires the lock; no outer lock is held here.
        taskset -c 15 "$binary" 101 0 0 20 1 2 1 | tee "asymmetric-${variant}-repeat${run}.jsonl"
    done
done
ctest --test-dir build -R asymmetric -V -j1 | tee asymmetric-correctness.log
sha256sum build/basis_balanced_masked_bench build/basis_asymmetric_{mixed2_3,greedy3_2}_{sparse,masked}_{bench,test} | tee asymmetric-binaries.sha256
