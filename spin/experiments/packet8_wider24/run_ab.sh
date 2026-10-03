#!/usr/bin/env bash
# Three-way matched campaign. Locks cover tests, warmups, and every timed call.
set -euo pipefail
if (( $# < 4 )); then
    echo 'usage: run_ab.sh nibble-control byte-control candidate fresh-output [calls=2001] [cpu=15] [seeds...]' >&2
    exit 2
fi
nibble=$1
byte=$2
candidate=$3
output=$4
calls=${5:-2001}
cpu=${6:-15}
shift "$(( $# < 6 ? $# : 6 ))"
seeds=("$@")
if (( ${#seeds[@]} == 0 )); then seeds=(41 113 257 997); fi
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
    sha256sum "$nibble" "$byte" "$candidate"
    for seed in "${seeds[@]}"; do
        for kind in nibble byte candidate candidate byte nibble; do
            echo "sample $kind seed=$seed calls=$calls cpu=$cpu"
            case "$kind" in
                nibble) taskset -c "$cpu" "$nibble" 52 65536 "$seed" "$calls" 0 0;;
                byte) taskset -c "$cpu" "$byte" 65536 "$seed" "$calls" 3 0;;
                candidate) taskset -c "$cpu" "$candidate" 65536 "$seed" "$calls" 0 0;;
            esac
        done
    done
} | tee "$output"
