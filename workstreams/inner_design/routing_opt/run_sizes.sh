#!/usr/bin/env bash
set -euo pipefail
cd "${1:?isolated build root}"
for m in 16 18; do
    for r in 1 2 3; do
        variants=(asymmetric_greedy3_2_sparse routeopt_pages_bothwrite)
        if [[ $r == 2 ]]; then variants=(routeopt_pages_bothwrite asymmetric_greedy3_2_sparse); fi
        for v in "${variants[@]}"; do
            flock -w 45 /tmp/bare-spin-benchmark.lock true
            "build/basis_${v}_bench" 101 0 0 "$m" 1 2 1 \
                2>"routing-opt/sizes-${m}-${v}-${r}.log" | tee "routing-opt/sizes-${m}-${v}-${r}.jsonl"
        done
    done
done
