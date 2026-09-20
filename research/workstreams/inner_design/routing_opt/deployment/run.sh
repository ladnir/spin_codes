#!/usr/bin/env bash
set -euo pipefail
cd "${1:?isolated repository root}"
mkdir -p routing-opt/deployment
ctest --test-dir deployment-build -V -j1 >routing-opt/deployment/correctness.log
for m in 16 18 20; do
    for repeat in 1 2 3; do
        variants=(supported_off supported_on asymmetric_off asymmetric_on)
        if [[ $repeat == 2 ]]; then variants=(asymmetric_on asymmetric_off supported_on supported_off); fi
        for variant in "${variants[@]}"; do
            flock -w 45 /tmp/bare-spin-benchmark.lock true
            "deployment-build/${variant}_bench" 101 0 0 "$m" 1 2 1 \
                | tee "routing-opt/deployment/${variant}-m${m}-r${repeat}.jsonl"
        done
    done
done
sha256sum deployment-build/*_bench deployment-build/*_test >routing-opt/deployment/binaries.sha256
sha256sum workstreams/bare_bch_rm2sub/{Spin.cpp,Spin.h,WorkspaceRouting.h,Inner.h,benchmark.cpp,CMakeLists.txt,correctness.cpp} \
    workstreams/inner_design/routing_opt/deployment/{AsymmetricSpin.cpp,AsymmetricMap.h,CMakeLists.txt,page_test.cpp} \
    workstreams/inner_design/asymmetric/AsymmetricInner.h \
    workstreams/rate_quarter_bch/implementation/{correctness.cpp,generated/QuarterCircuit.cpp,generated/QuarterCircuit.h} \
    >routing-opt/deployment/sources.sha256
