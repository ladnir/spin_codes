#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
exec 9>/tmp/hypercat-global-benchmark.lock
flock -x 9
exec 8>/tmp/spin-char2-benchmark.lock
flock -x 8
case "${1:-phases}" in
  phases)
    for phase in outer inner; do
      for binary in baseline decoded source both; do
        echo "case,$binary,phase,$phase"
        taskset -c 2 "./$binary" --bench 501 "$phase" 20 nt 4 4
      done
    done
    ;;
  matched)
    for round in 1 2 3; do
      if (( round%2 )); then order=(baseline decoded source both); else order=(both source decoded baseline); fi
      for binary in "${order[@]}"; do
        echo "round,$round,case,$binary"
        taskset -c 2 "./$binary" --bench 501 full 20 nt 4 4
      done
    done
    ;;
  address)
    for round in 1 2 3; do
      if (( round%2 )); then order=(decoded address); else order=(address decoded); fi
      for binary in "${order[@]}"; do
        echo "round,$round,case,$binary"
        taskset -c 2 "./$binary" --bench 501 full 20 nt 4 4
      done
    done
    ;;
  batched)
    for round in 1 2 3; do
      if (( round%2 )); then order=(decoded batched-decoded); else order=(batched-decoded decoded); fi
      for binary in "${order[@]}"; do
        echo "round,$round,case,$binary"
        taskset -c 2 "./$binary" --bench 501 full 20 nt 4 4
      done
    done
    ;;
esac
