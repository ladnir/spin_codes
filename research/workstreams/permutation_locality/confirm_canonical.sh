#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/canonical-confirm"
for seed in 1 17; do
  for repetition in 1 2 3; do
    modes=(control prior-wide block-stream block-restrict)
    if [[ "$repetition" == 2 ]]; then modes=(block-restrict block-stream prior-wide control); fi
    for name in "${modes[@]}"; do
      case "$name" in
        control) args=(1 0 tiled 1);;
        prior-wide) args=(16 1 wide-stream 1);;
        *) args=(4 1 "$name" 4);;
      esac
      receipt="$root/measurements/canonical-confirm/$name-s$seed-r$repetition.csv"
      taskset -c 15 "$root/build/locality" 20 "${args[0]}" "${args[1]}" "${args[2]}" "$seed" 101 "${args[3]}" > "$receipt"
      cat "$receipt"
    done
  done
done
