#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass isolated build root}
cd "$root"
for repeat in 1 2 3; do
    variants=(balanced greedy3_2 mixed2_3)
    if [[ $repeat == 2 ]]; then variants=(mixed2_3 greedy3_2 balanced); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        taskset -c 15 "inner-build/inner_${variant}_bench" | tee "asymmetric-kernel-${variant}-repeat${repeat}.jsonl"
    done
done
sha256sum inner-build/inner_{balanced,greedy3_2,mixed2_3}_bench | tee asymmetric-kernel-binaries.sha256
