#!/usr/bin/env bash
# Baseline: compile bench/setup.cpp without SPIN_BENCH_REUSE against the old lib.
# Candidate: cmake -DSPIN_BUILD_BENCHMARKS=ON, with matching Release/tuning flags.
set -euo pipefail
root=$(realpath "${1:?comparison directory containing baseline and candidate builds}")
cpu=${2:-15}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || { echo 'Another encoder benchmark holds the lock' >&2; exit 1; }
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || { echo 'Another SPIN benchmark holds the lock' >&2; exit 1; }
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || { echo 'Another Hypercat benchmark holds the lock' >&2; exit 1; }
mkdir -p "$root/results"
for pass in 1 2 3; do
    order='baseline fresh reuse'
    if [[ $pass == 2 ]]; then order='reuse fresh baseline'; fi
    for mode in $order; do
        case $mode in
            baseline) command=("$root/baseline/bench");;
            fresh) command=("$root/candidate/build/spin_setup_bench");;
            reuse) command=("$root/candidate/build/spin_setup_bench" --reuse);;
        esac
        taskset -c "$cpu" "${command[@]}" > "$root/results/$mode-$pass.csv"
        printf '%s pass %s\n' "$mode" "$pass"
        cat "$root/results/$mode-$pass.csv"
    done
done
