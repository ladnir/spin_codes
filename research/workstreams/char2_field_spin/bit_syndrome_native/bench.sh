#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec 9>/tmp/hypercat-global-benchmark.lock
flock -x 9
exec 8>/tmp/spin-char2-benchmark.lock
flock -x 8
for round in 1 2 3; do
  if (( round%2 )); then order=(g16 bit); else order=(bit g16); fi
  for binary in "${order[@]}"; do
    echo "round,$round,case,$binary"
    taskset -c 2 "./$binary" --bench 501 full 20 nt 16 4
  done
done
