#!/usr/bin/env bash
set -euo pipefail
binary=$1; output=$2; k=$3; seed=$4; calls=$5; cpu=$6; memory=$7; phases=$8
shift 8
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
    taskset -c "$cpu" "$binary" "$k" "$seed" "$calls" "$mode" "$memory" "$phases"
  done
} | tee "$output"
