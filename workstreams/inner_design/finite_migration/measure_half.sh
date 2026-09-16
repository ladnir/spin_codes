#!/usr/bin/env bash
set -euo pipefail
root=${1:?existing isolated build root required}
cd "$root"
sha256sum -c measurements/sources.sha256
sha256sum -c measurements/binaries.sha256
result=$(mktemp -d /tmp/imt-finite-half-XXXXXX)
printf 'Measurement directory: %s\n' "$result"
cp measurements/sources.sha256 "$result/sources.sha256"
cp measurements/binaries.sha256 "$result/binaries.sha256"
cp "$0" "$result/measure_half.sh"
{ uname -a; lscpu; c++ --version; } > "$result/environment.txt"
ctest --test-dir build --output-on-failure -j1 > "$result/correctness.log"
for variant in baseline baseline_tuned sparse_pages; do
    cp "build/CMakeFiles/${variant}.dir/flags.make" "$result/flags-${variant}.txt"
done
for m in 16 18 20; do
    for repeat in 1 2 3; do
        variants=(baseline baseline_tuned sparse_pages)
        if [[ $repeat == 2 ]]; then variants=(sparse_pages baseline_tuned baseline); fi
        for variant in "${variants[@]}"; do
            # The executable itself holds this lock throughout the process
            # and rejects other benchmark executables before measuring.
            flock -w 45 /tmp/bare-spin-benchmark.lock true
            taskset -c 15 "build/${variant}_bench" 101 0 0 "$m" 0 2 1 \
                | tee "$result/${variant}-m${m}-r${repeat}.jsonl"
        done
    done
done
sha256sum -c "$result/sources.sha256"
sha256sum -c "$result/binaries.sha256"
printf 'Completed measurements: %s\n' "$result"
