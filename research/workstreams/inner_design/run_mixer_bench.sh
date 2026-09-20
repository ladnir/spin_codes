#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass the isolated remote build root}
cd "$root"
# Each executable acquires the shared lock and rejects other benchmark
# processes. Never wrap this script in a second lock on that same file.
for run in 1 2 3; do
    taskset -c 15 baseline/spin_benchmark 101 0 0 20 1 2 1 | tee "mixer-baseline-repeat${run}.jsonl"
    for variant in rows sparse masked; do
        taskset -c 15 "build/basis_mixer_${variant}_bench" 101 0 0 20 1 2 1 | tee "mixer-${variant}-repeat${run}.jsonl"
    done
done
