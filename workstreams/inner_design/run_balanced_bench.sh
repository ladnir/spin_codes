#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass the isolated remote build root}
cd "$root"
# The executable owns the shared lock. Do not hold a second lock around it.
# The preliminary wait yields to another task; the executable rechecks races.
for run in 1 2 3; do
    variants=(baseline sparse masked)
    if [[ $run == 2 ]]; then variants=(masked sparse baseline); fi
    if [[ $run == 3 ]]; then variants=(baseline masked sparse); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        binary="build/basis_balanced_${variant}_bench"
        if [[ $variant == baseline ]]; then binary=baseline/spin_benchmark; fi
        taskset -c 15 "$binary" 101 0 0 20 1 2 1 | tee "balanced-${variant}-repeat${run}.jsonl"
    done
done
ctest --test-dir build -R balanced -V -j1 | tee balanced-correctness.log
sha256sum baseline/spin_benchmark build/basis_balanced_{sparse,masked}_{bench,test} | tee balanced-binaries.sha256
