#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec 9>/tmp/hypercat-global-benchmark.lock
flock -x 9
exec 8>/tmp/spin-char2-benchmark.lock
flock -x 8
for round in 1 2 3; do
  case "$round" in
    1|3) order=(bit lambda);;
    2) order=(lambda bit);;
  esac
  for binary in "${order[@]}"; do
    echo "round,$round,case,$binary"
    taskset -c 2 "./$binary" --bench 501 full 20 nt 16 4
  done
done
for binary in bit lambda; do
  echo "outer-only,case,$binary"
  taskset -c 2 "./$binary" --bench 501 outer 20 nt 16 4
done
