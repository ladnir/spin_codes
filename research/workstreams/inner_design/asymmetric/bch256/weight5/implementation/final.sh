#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass isolated remote root}
cd "$root"
for repeat in 1 2 3; do
    variants=(baseline baseline_tuned sparse sparse_pages sparse_tuned grouped_tuned)
    if [[ $repeat == 2 ]]; then variants=(grouped_tuned sparse_tuned sparse_pages sparse baseline_tuned baseline); fi
    if [[ $repeat == 3 ]]; then variants=(sparse_pages baseline sparse_tuned baseline_tuned grouped_tuned sparse); fi
    for variant in "${variants[@]}"; do
        flock -w 45 /tmp/bare-spin-benchmark.lock true
        taskset -c 15 "build/${variant}_bench" 101 2048 0 20 0 2 1 | tee "measurements/final-${variant}-${repeat}.jsonl"
    done
done
ctest --test-dir build -V -j1 > measurements/correctness-verbose.log
sha256sum build/*_bench build/*_test > measurements/binaries-final.sha256
for variant in baseline baseline_tuned sparse sparse_pages sparse_tuned masked_tuned shared_tuned grouped_tuned sparse_pf128; do
    cp "build/CMakeFiles/${variant}.dir/flags.make" "measurements/flags-${variant}.txt"
done
find workstreams constructions -type f \( -name '*.cpp' -o -name '*.h' -o -name 'CMakeLists.txt' \) -exec sha256sum '{}' \; > measurements/sources.sha256
