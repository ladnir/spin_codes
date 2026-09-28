#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build-sanitize}")
cmake --build "$root/build-sanitize" --target locality -j2
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
for seed in 1 17; do
  "$root/build-sanitize/locality" 14 4 1 verify "$seed" 1 1
  "$root/build-sanitize/locality" 14 4 1 random-gfni-physical "$seed" 1 1
done
ctest --test-dir "$root/build-sanitize/spin" --output-on-failure
