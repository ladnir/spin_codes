#!/usr/bin/env bash
# Run against a clean archive of tracked source, not the experiment checkout.
set -euo pipefail
root=${1:?clean repository root}
upstream=${2:-}
result=${3:-$root/workstreams/spin_optimized/measurements/cleanup}
mkdir -p "$result"
cmake -S "$root/workstreams/spin_optimized" -B "$root/build-clean-transpose" \
  -DCMAKE_BUILD_TYPE=Release -DSPIN_BCH_AVX512=ON -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_TUNE=znver4 > "$result/configure.log" 2>&1
cmake --build "$root/build-clean-transpose" -j2 > "$result/build.log" 2>&1
test ! -d "$root/build-clean-transpose/bidirectional"
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir "$root/build-clean-transpose" --output-on-failure -j1 | tee "$result/tests.log"
python3 "$root/workstreams/spin_optimized/isa_check.py" "$root/build-clean-transpose" > "$result/isa.log"
if [[ -n $upstream ]]; then
  cmake -S "$root/workstreams/spin_optimized" -B "$root/build-clean-optional" \
    -DCMAKE_BUILD_TYPE=Release -DSPIN_BCH_AVX512=ON -DSPIN_BUILD_BENCHMARK=ON \
    -DSPIN_TUNE=znver4 -DSPIN_BIDIRECTIONAL_SOURCE="$upstream" > "$result/optional-configure.log" 2>&1
  cmake --build "$root/build-clean-optional" -j2 --target spin_bidirectional_test \
    spin_k16_bidirectional_test spin_bidirectional_benchmark > "$result/optional-build.log" 2>&1
  flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
    ctest --test-dir "$root/build-clean-optional" -R bidirectional --output-on-failure -j1 | tee "$result/optional-tests.log"
  "$root/build-clean-optional/spin_bidirectional_benchmark" 20 auto 3 1 0 0 12819 81920 > "$result/optional-length.json"
  python3 -c 'import json,sys;v=json.load(open(sys.argv[1]));assert v["K"]==81920 and v["m"] is None' "$result/optional-length.json"
fi
if cmake -S "$root/workstreams/spin_optimized" -B "$root/build-clean-invalid" \
  -DSPIN_BUILD_WIDE=ON > "$result/reject-forward.log" 2>&1; then
  echo 'Forward option without source was accepted' >&2;exit 1
fi
grep -q 'Forward/wide options require' "$result/reject-forward.log"
# Smoke checks of JSON metadata, not a performance comparison. Never concurrent.
for k in 65536 81920; do
  "$root/build-clean-transpose/spin_benchmark" 20 auto 3 1 0 0 12819 "$k" > "$result/metadata-$k.json"
done
python3 -c 'import json,sys;from pathlib import Path;p=Path(sys.argv[1]);a=json.loads((p/"metadata-65536.json").read_text());b=json.loads((p/"metadata-81920.json").read_text());assert (a["K"],a["m"])==(65536,16);assert (b["K"],b["m"])==(81920,None);print("benchmark metadata PASS")' "$result"
