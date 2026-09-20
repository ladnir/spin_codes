#!/usr/bin/env bash
set -euo pipefail
cd "${1:?isolated build root}"
mkdir -p routing-opt
ctest --test-dir build -R '^basis_routeopt_.*_test$' -V -j1 | tee routing-opt/correctness.log
for v in asymmetric_greedy3_2_sparse routeopt_scatter0 routeopt_scatter32 routeopt_scatter128 routeopt_write32 routeopt_write128 routeopt_gather0 routeopt_gather64; do
    flock -w 45 /tmp/bare-spin-benchmark.lock true
    "build/basis_${v}_bench" 31 4096 0 20 1 2 1 | tee "routing-opt/screen-${v}.jsonl"
done
sha256sum build/basis_routeopt*_{bench,test} build/basis_asymmetric_greedy3_2_sparse_bench > routing-opt/binaries.sha256
