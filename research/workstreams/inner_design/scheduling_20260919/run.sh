#!/usr/bin/env bash
set -euo pipefail
root=${1:?isolated root}
cd "$root"
here=workstreams/inner_design/scheduling_20260919
mkdir -p "$here/measurements"
result="$here/measurements"
{ uname -a; lscpu; c++ --version; cat /sys/devices/system/cpu/cpu15/cpufreq/scaling_cur_freq;
  cat /sys/devices/system/cpu/cpu15/cpufreq/scaling_governor; cat /sys/devices/system/cpu/cpufreq/boost; } > "$result/environment.txt"
find workstreams constructions -type f \( -name '*.h' -o -name '*.cpp' -o -name '*.py' -o -name '*.sh' -o -name CMakeLists.txt \) -print0 | sort -z | xargs -0 sha256sum > "$result/sources.sha256"
sha256sum build-scheduling/*_bench build-scheduling/*_test > "$result/binaries.sha256"
# Hold both locks during validation too; timed executables acquire the same
# locks themselves, in this same order, and scan /proc for other benchmarks.
flock -n /tmp/prindal-addition-encoder-benchmark.lock flock -n /tmp/bare-spin-benchmark.lock \
  ctest --test-dir build-scheduling --output-on-failure -j1 | tee "$result/correctness.log"
for name in control pf0 pf128 pf256 group4 group8 decode decode_pf emit_pf control; do
  suffix=pilot
  if [[ $name == control && -f "$result/pilot-control.jsonl" ]]; then suffix=pilot-end; fi
  "build-scheduling/${name}_bench" 31 2048 0 20 0 2 1 | tee "$result/${suffix}-${name}.jsonl"
done
