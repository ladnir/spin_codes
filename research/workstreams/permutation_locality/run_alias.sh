#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing build}")
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
mkdir -p "$root/measurements/alias"
for pass in 1 2 3; do
  binaries='locality-before-restrict locality'
  if ((pass==2)); then binaries='locality locality-before-restrict'; fi
  for binary in $binaries; do
    taskset -c 15 "$root/build/$binary" 20 4 1 column-profile 17 101 > "$root/measurements/alias/$binary-$pass.csv"
    cat "$root/measurements/alias/$binary-$pass.csv"
  done
done
