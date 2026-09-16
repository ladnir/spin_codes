#!/usr/bin/env bash
set -euo pipefail
cd "${1:?isolated build root}"
ctest --test-dir build -R 'basis_routeopt_pages_bothwrite_test' -V -j1 >routing-opt/correctness-combined.log
for r in 1 2 3; do
    for v in asymmetric_greedy3_2_sparse routeopt_pages_both routeopt_pages_bothwrite; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        "build/basis_${v}_bench" 101 4096 0 20 1 2 1 \
            2>"routing-opt/combined-${v}-${r}.log" | tee "routing-opt/combined-${v}-${r}.jsonl"
    done
done
sha256sum build/basis_routeopt_pages_bothwrite_{bench,test} >routing-opt/combined-binaries.sha256
sha256sum workstreams/inner_design/generated/routeopt_pages_bothwrite/* >routing-opt/combined-sources.sha256
