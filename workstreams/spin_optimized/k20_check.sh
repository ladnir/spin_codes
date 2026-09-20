#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/k20-final
mkdir -p "$result"
build=build-k20-final
cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
 -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
 -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" -DSPIN_K20_EXPERIMENT=0 > "$result/build.log" 2>&1
cmake --build "$build" -j2 >> "$result/build.log" 2>&1
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir "$build" --output-on-failure -j1 | tee "$result/test.log"
cp "$build/Testing/Temporary/LastTest.log" "$result/test-details.log"
python3 workstreams/spin_optimized/isa_check.py "$build" | tee "$result/isa.log"
for seed in 1 17; do
 for repeat in 1 2 3; do
  kinds=(spin_benchmark spin_bidirectional_benchmark)
  if [[ $repeat == 2 ]]; then kinds=(spin_bidirectional_benchmark spin_benchmark); fi
  for kind in "${kinds[@]}"; do
   "$build/$kind" 20 auto 101 "$seed" > "$result/control-$kind-m20-s$seed-r$repeat.json"
  done
 done
done
find "$build/generated" "$build/bidirectional" -type f -print0 | sort -z | xargs -0 sha256sum > "$result/sources.sha256"
sha256sum "$build/spin_benchmark" "$build/spin_bidirectional_benchmark" > "$result/binaries.sha256"
python3 workstreams/spin_optimized/k18_summary.py "$result" | tee "$result/summary.txt"
