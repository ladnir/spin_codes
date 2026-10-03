#!/usr/bin/env bash
# Hold all shared host locks for the entire campaign. Never run concurrently.
set -euo pipefail
if (( $# < 3 )); then
    echo 'usage: run_ab.sh baseline candidate fresh-output [calls=1001] [cpu=15] [seeds...]' >&2
    exit 2
fi
baseline=$1
candidate=$2
output=$3
calls=${4:-1001}
cpu=${5:-15}
candidate_mode=${PACKET8_MODE:-0}
shift "$(( $# < 5 ? $# : 5 ))"
seeds=("$@")
if (( ${#seeds[@]} == 0 )); then seeds=(1 17 41 113); fi
if [[ -e "$output" ]]; then echo 'fresh output required' >&2; exit 2; fi
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
exec 8>/tmp/bare-spin-benchmark.lock
exec 7>/tmp/hypercat-benchmark.lock
exec 6>/tmp/hypercat-global-benchmark.lock
flock 9
flock 8
flock 7
flock 6
{
    date -u
    sha256sum "$baseline" "$candidate"
    for seed in "${seeds[@]}"; do
        # Reverse ordering within each seed: no clocks between encoder stages.
        for kind in baseline candidate candidate baseline; do
            echo "sample $kind seed=$seed calls=$calls cpu=$cpu"
            if [[ "$kind" == baseline ]]; then
                taskset -c "$cpu" "$baseline" 52 65536 "$seed" "$calls" 0 0
            else
                taskset -c "$cpu" "$candidate" 65536 "$seed" "$calls" "$candidate_mode" 0
            fi
        done
    done
} | tee "$output"
