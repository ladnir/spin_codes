#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec 9>/tmp/hypercat-global-benchmark.lock
flock -x 9
exec 8>/tmp/spin-char2-benchmark.lock
flock -x 8
for round in 1 2 3; do
  case "$round" in
    1) order=(bit cross1 cross2);;
    2) order=(cross2 bit cross1);;
    3) order=(cross1 cross2 bit);;
  esac
  for binary in "${order[@]}"; do
    echo "round,$round,case,$binary"
    taskset -c 2 "./$binary" --bench 501 full 20 nt 16 4
  done
done
for binary in bit cross1 cross2; do
  echo "outer-only,case,$binary"
  taskset -c 2 "./$binary" --bench 501 outer 20 nt 16 4
done
