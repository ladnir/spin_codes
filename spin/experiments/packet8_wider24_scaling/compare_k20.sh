#!/usr/bin/env bash
set -euo pipefail
new=$1; old=$2; output=$3; calls=${4:-301}; cpu=${5:-15}
[[ ! -e "$output" ]] || { echo 'fresh output required' >&2; exit 1; }
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
exec 8>/tmp/bare-spin-benchmark.lock
exec 7>/tmp/hypercat-benchmark.lock
exec 6>/tmp/hypercat-global-benchmark.lock
flock -x 9; flock -x 8; flock -x 7; flock -x 6
{
  date -u
  sha256sum "$new" "$old"
  taskset -c "$cpu" "$old" check 1048576 11 20
  for seed in 271 419 557 863; do
    order='new old old new'
    [[ "$seed" != 419 && "$seed" != 863 ]] || order='old new new old'
    for candidate in $order; do
      if [[ "$candidate" == new ]]; then
        taskset -c "$cpu" "$new" 1048576 "$seed" "$calls" 4 huge 0
      else
        taskset -c "$cpu" "$old" fused-nt 1048576 20 "$seed" "$calls" huge
      fi
    done
  done
} | tee "$output"
