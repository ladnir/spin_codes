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
mkdir -p "$root/results"
for pass in 1 2 3; do
    order='0 4 6 8'
    if [[ $pass == 2 ]]; then order='8 6 4 0'; fi
    for rounds in $order; do
        for reuse in fresh reuse; do
            args=();if [[ $reuse == reuse ]]; then args+=(--reuse);fi
            name="full-$rounds-$reuse-$pass"
            taskset -c "$cpu" "$root/build/spin_feistel_bench" "$rounds" "${args[@]}" > "$root/results/$name.csv"
            printf '%s\n' "$name";cat "$root/results/$name.csv"
        done
    done
    for rounds in $order; do
        for mode in materialized online; do
            if [[ $rounds == 0 && $mode == online ]]; then continue;fi
            for reuse in fresh reuse; do
                args=();if [[ $reuse == reuse ]]; then args+=(--reuse);fi
                name="direct-$rounds-$mode-$reuse-$pass"
                taskset -c "$cpu" "$root/build/spin_feistel_direct" "$rounds" "$mode" "${args[@]}" > "$root/results/$name.csv"
                printf '%s\n' "$name";cat "$root/results/$name.csv"
            done
        done
    done
done
