#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
result=workstreams/inner_design/scheduling_20260919/measurements
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-scheduling --output-on-failure -j1 -R '^(pf8|pf16|pf64|group4p32)$' | tee "$result/refine-correctness.log"
for name in control pf8 pf16 pf64 group4p32; do
  "build-scheduling/${name}_bench" 31 2048 0 20 0 2 1 | tee "$result/refine-${name}.jsonl"
done
for tile in 1024 2048 4096; do
  for layout in 0 1; do
    build-scheduling/control_bench 31 "$tile" "$layout" 20 0 2 1 | tee "$result/layout-t${tile}-l${layout}.jsonl"
  done
done
build-scheduling/profile_bench 31 2048 0 20 0 2 1 > "$result/profile.jsonl" 2> "$result/profile.log"
cat "$result/profile.log"
sha256sum build-scheduling/*_bench build-scheduling/*_test > "$result/binaries-refine.sha256"
