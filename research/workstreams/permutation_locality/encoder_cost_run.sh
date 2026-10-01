#!/usr/bin/env bash
# Research diagnostics only; reuse the immutable measured encoder/object files.
set -euo pipefail
root=$(realpath "${1:?diagnostic directory containing encoder_cost_probe.cpp}")
phase=${2:-footprint}
prior=/tmp/spin-bch-tune-fN0NiP
reference=/tmp/spin-joint-9n57TT
t64=/tmp/spin-bch-tune-vKm7EZ/gfni-t64
s19=/tmp/spin-bch-tune-Q1wko3/gfni-r4
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -n 7 || exit 75
logs=$(mktemp -d "$root/cost-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
case "$phase" in
  build)
    g++ -O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni -mtune=znver4 \
      -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1 \
      -I"$reference/spin/src/kernels" -I"$reference/spin/include" -I"$prior" \
      "$root/encoder_cost_probe.cpp" "$prior/PackedCoeff1.o" \
      /tmp/spin-packed-mixer-NnskDD/PackedMixer.o /tmp/spin-bch-compare-2e2UPS/BchAvx512.o \
      -o "$root/cost-probe" 2>&1 | tee "$logs/build.txt"
    sha256sum "$root/encoder_cost_probe.cpp" "$root/cost-probe" "$prior/PackedCoeff1.o" "$t64" "$s19" > "$logs/hashes.txt"
    ;;
  coarse)
    for seed in 1 17; do
      taskset -c 15 "$s19" 20 "$seed" 31 1 1 > "$logs/phase-$seed-control.csv"
      taskset -c 15 "$t64" 20 "$seed" 31 1 > "$logs/phase-$seed-t64.csv"
    done
    ;;
  footprint)
    for seed in 1 17; do
      taskset -c 15 "$root/cost-probe" "$seed" 101 > "$logs/footprint-$seed.csv"
    done
    ;;
  sample)
    perf record -o "$logs/cycles.data" -e cycles:u -F 1999 -- \
      taskset -c 15 "$t64" 20 17 5001 0 > "$logs/run.txt" 2> "$logs/perf.txt"
    perf report -i "$logs/cycles.data" --stdio --no-children --percent-limit 0.2 --sort dso,symbol > "$logs/report.txt"
    perf annotate -i "$logs/cycles.data" --stdio --no-source --percent-type global-period > "$logs/annotate.txt"
    ;;
  counters)
    for footprint in 1 2048; do
      perf stat -x, -e cycles:u,instructions:u,cache-references:u,cache-misses:u \
        -o "$logs/counters-$footprint.csv" -- taskset -c 15 "$root/cost-probe" 1 1001 "$footprint" \
        > "$logs/run-$footprint.txt"
    done
    ;;
  *) printf 'Unknown phase: %s\n' "$phase" >&2; exit 2;;
esac
