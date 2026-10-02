#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
result=workstreams/inner_design/bch_20260919/measurements
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir build-bch --output-on-failure -j1 -R '^(share3|share4|share8|four)$' | tee "$result/refine-correctness.log"
for name in control dfs share3 share4 share8 four control; do
  suffix=refine
  if [[ $name == control && -f "$result/refine-control.jsonl" ]]; then suffix=refine-end; fi
  "build-bch/${name}_bench" 31 2048 0 20 0 2 1 > "$result/${suffix}-${name}.jsonl"
  python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[1],d["median_ms"],d["output_hash"])' "$result/${suffix}-${name}.jsonl"
done
find build-bch -name '*.su' -exec cat {} + > "$result/stack-usage-refine.txt"
