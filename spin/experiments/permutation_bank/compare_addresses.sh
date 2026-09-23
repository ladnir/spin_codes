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
mkdir -p "$root/results-bank-addresses"
for pass in 1 2 3;do
    families='bank16 bank64 feistel6';if [[ $pass == 2 ]];then families='feistel6 bank64 bank16';fi
    for family in $families;do
        for mode in scalar batch;do
            file="$root/results-bank-addresses/$family-$mode-$pass.csv"
            taskset -c "$cpu" "$root/build/spin_bank_addresses" "$family" "$mode" > "$file"
            cat "$file"
        done
    done
done
