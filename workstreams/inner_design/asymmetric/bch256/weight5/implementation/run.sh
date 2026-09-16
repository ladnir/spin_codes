#!/usr/bin/env bash
set -euo pipefail
root=${1:?pass isolated remote root}
cd "$root"
mkdir -p measurements
ctest --test-dir build --output-on-failure -j1 | tee measurements/correctness.log
sha256sum build/*_bench build/*_test > measurements/binaries.sha256
{ uname -a; lscpu; c++ --version; } > measurements/environment.txt
for variant in baseline baseline_tuned sparse sparse_pages sparse_tuned masked_tuned shared_tuned sparse_pf128; do
    # The executable checks for other benchmark processes and takes the shared lock.
    flock -w 45 /tmp/bare-spin-benchmark.lock true
    taskset -c 15 "build/${variant}_bench" 51 0 0 20 0 2 1 | tee "measurements/screen-${variant}.jsonl"
done
