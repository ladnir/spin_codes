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
mkdir -p "$root/results-phases"
for pass in 1 2 3;do
    rounds='0 4 6 8';if [[ $pass == 2 ]];then rounds='8 6 4 0';fi
    for r in $rounds;do
        for mode in materialized online;do
            if [[ $r == 0 && $mode == online ]];then continue;fi
            file="$root/results-phases/$r-$mode-$pass.csv"
            taskset -c "$cpu" "$root/build/spin_feistel_phases" "$r" "$mode" > "$file"
            cat "$file"
        done
    done
done
