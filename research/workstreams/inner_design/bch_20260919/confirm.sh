#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
result=workstreams/inner_design/bch_20260919/measurements
find workstreams constructions -type f \( -name '*.h' -o -name '*.cpp' -o -name '*.py' -o -name '*.sh' -o -name CMakeLists.txt \) -print0 | sort -z | xargs -0 sha256sum > "$result/sources-final.sha256"
sha256sum build-bch/*_bench build-bch/*_test build-bch/packedshare4_tiles > "$result/binaries-final.sha256"
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir build-bch --output-on-failure -j1 | tee "$result/final-correctness.log"
for seed in 1 17; do
 for repeat in 1 2 3; do
  names=(control packedfour packedshare3 packedshare4)
  if [[ $repeat == 2 ]]; then names=(packedshare4 packedshare3 packedfour control); fi
  for name in "${names[@]}"; do
   file="$result/final-${name}-m20-s${seed}-r${repeat}.jsonl"
   "build-bch/${name}_bench" 101 2048 0 20 0 2 1 "$seed" > "$file"
   python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[1],d["median_ms"],d["output_hash"])' "$file"
  done
 done
done
for m in 16 18; do
 for repeat in 1 2 3; do
  names=(control packedshare4)
  if [[ $repeat == 2 ]]; then names=(packedshare4 control); fi
  for name in "${names[@]}"; do
   file="$result/final-${name}-m${m}-s1-r${repeat}.jsonl"
   "build-bch/${name}_bench" 101 0 0 "$m" 0 2 1 > "$file"
   python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[1],d["median_ms"])' "$file"
  done
 done
done
for mode in hot stream; do
 for repeat in 1 2 3; do
  names=(control packedshare4)
  if [[ $repeat == 2 ]]; then names=(packedshare4 control); fi
  for name in "${names[@]}"; do
   "build-bch/${name}_outer_bench" "$mode" > "$result/outer-${name}-${mode}-r${repeat}.jsonl"
  done
 done
done
find build-bch -name '*.su' -exec cat {} + > "$result/stack-usage-final.txt"
sha256sum -c "$result/sources-final.sha256" > "$result/source-check-final.log"
sha256sum -c "$result/binaries-final.sha256" > "$result/binary-check-final.log"
