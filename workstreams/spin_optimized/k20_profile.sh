#!/usr/bin/env bash
set -euo pipefail
root=${1:?root}
cd "$root"
result=workstreams/spin_optimized/measurements/k20-profile
mkdir -p "$result"
build=build-k20-profile
cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
 -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
 -DSPIN_STAGE_PROFILE=ON -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/build.log" 2>&1
cmake --build "$build" -j2 >> "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/tests.log"
for seed in 1 17; do
 for repeat in 1 2; do
  for kind in spin_benchmark spin_bidirectional_benchmark; do
   for tile in 512 1024 2048; do
    SPIN_PROFILE=1 "$build/$kind" 20 auto 101 "$seed" "$tile" \
     > "$result/$kind-s$seed-r$repeat-t$tile.json" 2> "$result/$kind-s$seed-r$repeat-t$tile-stages.json"
   done
   build-k20-final/"$kind" 20 auto 101 "$seed" > "$result/control-$kind-s$seed-r$repeat.json"
  done
 done
done
# Sample the unchanged binary separately; setup samples are retained and labeled.
perf record -e cycles:u -F 999 -o "$result/cycles.data" -- \
 build-k20-final/spin_bidirectional_benchmark 20 auto 501 1 > "$result/sampled.json" 2> "$result/perf-record.log"
perf report --stdio -i "$result/cycles.data" --no-children --percent-limit 0.5 > "$result/perf-report.txt"
perf annotate --stdio --no-source -i "$result/cycles.data" > "$result/annotate.txt"
sha256sum "$build/generated/Fast.cpp" "$build/bidirectional/Fast.cpp" \
 "$build/spin_benchmark" "$build/spin_bidirectional_benchmark" \
 build-k20-final/spin_bidirectional_benchmark > "$result/sources-binaries.sha256"
python3 workstreams/spin_optimized/profile_summary.py "$result" | tee "$result/summary.txt"
python3 workstreams/spin_optimized/annotate_summary.py "$result/annotate.txt" | tee "$result/assembly-summary.txt"
