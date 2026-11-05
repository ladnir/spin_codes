#!/usr/bin/env bash
set -euo pipefail
cd /tmp/hd-compact7-20261009-uot619
exec 9>/tmp/hypercat-global-benchmark.lock
exec 8>/tmp/spin-char2-benchmark.lock
flock -n 9
flock -n 8
for round in 1 2 3; do
    date -u
    if (( round % 2 )); then
        taskset -c 2 ./control --bench 501 full 20 nt 4 4
        taskset -c 2 ./compact7 --bench 501 full 20 nt 4 4
    else
        taskset -c 2 ./compact7 --bench 501 full 20 nt 4 4
        taskset -c 2 ./control --bench 501 full 20 nt 4 4
    fi
done
