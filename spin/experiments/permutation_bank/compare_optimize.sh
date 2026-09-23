#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?directory containing build}")
cpu=${2:-15}
binary=${3:-spin_bank_optimize}
label=${4:-baseline}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 1
mkdir -p "$root/results-bank-optimize-$label"
for ((pass=1;pass<=${SPIN_COMPARE_PASSES:-3};++pass));do
    modes=${SPIN_COMPARE_MODES:-'fixed composed composed4 composed2'}
    modes=${modes//,/ }
    if (( pass % 2 == 0 ));then modes=$(echo "$modes" | xargs -n1 | tac | xargs);fi
    for mode in $modes;do
        file="$root/results-bank-optimize-$label/$mode-$pass.csv"
        taskset -c "$cpu" "$root/build/$binary" "$mode" > "$file"
        cat "$file"
    done
done
