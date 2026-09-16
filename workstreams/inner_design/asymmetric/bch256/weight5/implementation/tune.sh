#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass isolated remote root}
cd "$root"
ctest --test-dir build --output-on-failure -j1 | tee measurements/correctness-final.log
sha256sum build/*_bench build/*_test > measurements/binaries-final.sha256
taskset -c 15 build/grouped_tuned_bench 51 0 0 20 0 2 1 | tee measurements/screen-grouped_tuned.jsonl
for variant in baseline_tuned sparse_pages sparse_tuned grouped_tuned; do
    for layout in 0 1; do
        for tile in 256 512 1024 2048 4096; do
            flock -w 45 /tmp/bare-spin-benchmark.lock true
            taskset -c 15 "build/${variant}_bench" 31 "$tile" "$layout" 20 0 2 1 | tee "measurements/tune-${variant}-${tile}-${layout}.jsonl"
        done
    done
done
