#!/usr/bin/env bash
# Build-only exact wide-feedback combinations.
set -euo pipefail
root=$(realpath -e -- "${1:?root}")
previous=$(realpath -e -- "${2:?previous build}")
control=$(realpath -e -- "${3:?retained control}")
repo=$(realpath -e -- "${4:?source repo}")
ws="$repo/research/workstreams/permutation_locality"
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -n 7 || exit 75
logs=$(mktemp -d "$root/feedback-build-XXXXXX")
printf 'Logs: %s\n' "$logs"
export PYTHONPATH="$root:$previous:$ws" PYTHONDONTWRITEBYTECODE=1
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
 -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
 -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
 -I"$root" -I"$previous" -I"$control" -I"$ws" -I"$repo/spin/src/kernels" -I"$repo/spin/include")
objects=("$previous/Transpose.o" "$previous/Gather.o" "$previous/RoutePacked.o" "$control/PackedMixer.o" "$control/BchAvx512.o")
build() {
 local name=$1 source=$2
 shift 2
 g++ "${flags[@]}" "$@" -c "$source" -o "$root/$name.o"
 g++ "$root/$name.o" "${objects[@]}" -o "$root/inner-size-$name"
}
python3 -B "$root/further_feedback_codegen.py" --output "$root/further_feedbackWide.h" 2> "$logs/generate.txt"
build feedback "$root/further_feedback_probe.cpp"
python3 -B "$root/further_control_codegen.py" "$root/further_feedback_probe.cpp" --kind align > "$root/FeedbackAligned.cpp"
build feedbackaligned "$root/FeedbackAligned.cpp" -DSPIN_INPUT_OFFSET=0
python3 -B "$root/further_memory_codegen.py" "$root/further_feedback_probe.cpp" > "$root/FeedbackMemory.cpp"
build feedbackhuge "$root/FeedbackMemory.cpp" -DSPIN_INPUT_OFFSET=0 -DSPIN_MEMORY_ALIGNMENT=2097152
sha256sum "$root"/inner-size* "$root"/*.py "$root"/*.h "$0" > "$logs/hashes.txt"
printf 'Built feedback candidates; no benchmark run.\n'
