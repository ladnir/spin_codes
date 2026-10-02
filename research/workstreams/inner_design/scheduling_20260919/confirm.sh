#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
here=workstreams/inner_design/scheduling_20260919
result="$here/measurements"
# Final confirmation is bound to its own sources and binaries. The preliminary
# screens used earlier generated files and remain selection evidence only.
find workstreams constructions -type f \( -name '*.h' -o -name '*.cpp' -o -name '*.py' -o -name '*.sh' -o -name CMakeLists.txt \) -print0 | sort -z | xargs -0 sha256sum > "$result/sources-final.sha256"
sha256sum build-scheduling/*_bench build-scheduling/*_test > "$result/binaries-final.sha256"
for seed in 1 17; do
  for repeat in 1 2 3; do
    names=(control pf16 tile1024)
    if [[ $repeat == 2 ]]; then names=(tile1024 pf16 control); fi
    for name in "${names[@]}"; do
      exe=$name;tile=2048
      if [[ $name == tile1024 ]]; then exe=control;tile=1024; fi
      "build-scheduling/${exe}_bench" 101 "$tile" 0 20 0 2 1 "$seed" > "$result/final-${name}-s${seed}-r${repeat}.jsonl"
      python3 -c 'import json,sys; d=json.load(open(sys.argv[1])); print(sys.argv[1],d["median_ms"],d["output_hash"])' "$result/final-${name}-s${seed}-r${repeat}.jsonl"
    done
  done
done
for tile in 1024 2048 4096; do
  build-scheduling/profile_bench 31 "$tile" 0 20 0 2 1 > "$result/profile-t${tile}.jsonl" 2> "$result/profile-t${tile}.log"
  cat "$result/profile-t${tile}.log"
done
sha256sum -c "$result/sources-final.sha256" > "$result/source-check-final.log"
sha256sum -c "$result/binaries-final.sha256" > "$result/binary-check-final.log"
