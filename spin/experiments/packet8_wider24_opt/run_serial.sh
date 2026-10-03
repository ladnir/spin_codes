#!/usr/bin/env bash
set -euo pipefail
binary=$1; output=$2; calls=${3:-1001}; cpu=${4:-15}; seed=${5:-41}
shift 5
[[ ! -e "$output" ]] || { echo 'fresh output required' >&2; exit 1; }
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
exec 8>/tmp/bare-spin-benchmark.lock
exec 7>/tmp/hypercat-benchmark.lock
exec 6>/tmp/hypercat-global-benchmark.lock
flock -x 9; flock -x 8; flock -x 7; flock -x 6
{
  date -u
  sha256sum "$binary"
  for mode in "$@"; do
    taskset -c "$cpu" "$binary" 65536 "$seed" "$calls" "$mode" 0
  done
} | tee "$output"
