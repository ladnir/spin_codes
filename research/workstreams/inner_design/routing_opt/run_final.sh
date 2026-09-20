#!/usr/bin/env bash
set -euo pipefail
cd "${1:?isolated build root}"
for repeat in 1 2 3; do
    variants=(control huge2048 write32 huge4096)
    if [[ $repeat == 2 ]]; then variants=(huge4096 write32 huge2048 control); fi
    for name in "${variants[@]}"; do
        tile=4096
        case $name in
            control) variant=asymmetric_greedy3_2_sparse;;
            write32) variant=routeopt_scalarwrite32;;
            huge2048) variant=routeopt_pages_both;tile=2048;;
            huge4096) variant=routeopt_pages_both;;
        esac
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        "build/basis_${variant}_bench" 101 "$tile" 0 20 1 2 1 \
            2>"routing-opt/final-${name}-${repeat}.log" | tee "routing-opt/final-${name}-${repeat}.jsonl"
    done
done
for name in control huge2048; do
    variant=asymmetric_greedy3_2_sparse;tile=4096
    if [[ $name == huge2048 ]]; then variant=routeopt_pages_both;tile=2048; fi
    flock -w 45 /tmp/bare-spin-benchmark.lock true
    perf stat -e cycles,instructions,dTLB-load-misses,cache-misses -o "routing-opt/perf-${name}.log" \
        "build/basis_${variant}_bench" 501 "$tile" 0 20 1 2 1 \
        >"routing-opt/perf-${name}.jsonl" 2>"routing-opt/perf-${name}-pages.log"
done
sha256sum build/basis_routeopt*_{bench,test} build/basis_asymmetric_greedy3_2_sparse_bench >routing-opt/binaries.sha256
sha256sum workstreams/inner_design/generated/routeopt*/*.cpp workstreams/inner_design/generated/routeopt*/*.h \
    workstreams/inner_design/generated/asymmetric_greedy3_2_sparse/* \
    workstreams/inner_design/asymmetric/AsymmetricInner.h workstreams/bare_bch_rm2sub/Spin.h \
    workstreams/bare_bch_rm2sub/Inner.h workstreams/bare_bch_rm2sub/benchmark.cpp \
    workstreams/rate_quarter_bch/implementation/generated/QuarterCircuit.* \
    workstreams/inner_design/CMakeLists.txt >routing-opt/compiled-sources.sha256
