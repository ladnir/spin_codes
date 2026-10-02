#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
here=workstreams/inner_design/bch_20260919
result="$here/measurements"
mkdir -p "$result"
{ uname -a; lscpu; c++ --version; cat /sys/devices/system/cpu/cpu15/cpufreq/scaling_governor; cat /sys/devices/system/cpu/cpufreq/boost; } > "$result/environment.txt"
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
 ctest --test-dir build-bch --output-on-failure -j1 | tee "$result/correctness.log"
sha256sum build-bch/*_bench build-bch/*_test > "$result/binaries.sha256"
find workstreams constructions -type f \( -name '*.cpp' -o -name '*.h' \) -print0 | sort -z | xargs -0 sha256sum > "$result/sources.sha256"
for name in control dfs chunk32 synth32 synth16 single control; do
  suffix=pilot
  if [[ $name == control && -f "$result/pilot-control.jsonl" ]]; then suffix=pilot-end; fi
  "build-bch/${name}_bench" 31 2048 0 20 0 2 1 > "$result/${suffix}-${name}.jsonl"
  python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[1],d["median_ms"],d["output_hash"])' "$result/${suffix}-${name}.jsonl"
done
find build-bch -name '*.su' -exec cat {} + > "$result/stack-usage.txt"
