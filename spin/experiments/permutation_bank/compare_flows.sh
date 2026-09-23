#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?directory containing build}")
cpu=${2:-15}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 1
mkdir -p "$root/results-bank-flows"
for pass in 1 2 3 4 5;do
    modes='fixed-exact fixed16 fixed64 fresh16 fresh64'
    if (( pass % 2 == 0 ));then modes='fresh64 fresh16 fixed64 fixed16 fixed-exact';fi
    for mode in $modes;do
        file="$root/results-bank-flows/$mode-$pass.csv"
        taskset -c "$cpu" "$root/build/spin_bank_flows" "$mode" > "$file"
        cat "$file"
    done
done
