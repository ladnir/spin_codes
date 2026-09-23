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
mkdir -p "$root/results-batched"
for pass in 1 2 3; do
    methods='scalar batched materialized'
    if [[ $pass == 2 ]]; then methods='materialized batched scalar';fi
    for reuse in fresh reuse; do
        args=();if [[ $reuse == reuse ]];then args+=(--reuse);fi
        taskset -c "$cpu" "$root/build/spin_feistel_direct" 0 materialized "${args[@]}" > "$root/results-batched/exact-$reuse-$pass.csv"
        for rounds in 4 6 8;do
            for method in $methods;do
                binary=spin_feistel_direct;mode=online
                if [[ $method == scalar ]];then binary=spin_feistel_scalar;fi
                if [[ $method == materialized ]];then mode=materialized;fi
                file="$root/results-batched/$method-$rounds-$reuse-$pass.csv"
                taskset -c "$cpu" "$root/build/$binary" "$rounds" "$mode" "${args[@]}" > "$file"
                printf '%s %s %s pass %s\n' "$method" "$rounds" "$reuse" "$pass";cat "$file"
            done
        done
    done
done
