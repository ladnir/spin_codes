#!/usr/bin/env bash
# Build-only follow-up. Root contains overlay generators; PREVIOUS is the
# measured creative winner's build, CONTROL supplies original reference objects.
set -euo pipefail
root=$(realpath -e -- "${1:?root}")
previous=$(realpath -e -- "${2:?previous build}")
control=$(realpath -e -- "${3:?retained control}")
repo=$(realpath -e -- "${4:?source repo}")
ws="$repo/research/workstreams/permutation_locality"
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -n 7 || exit 75
logs=$(mktemp -d "$root/control-build-XXXXXX")
printf 'Logs: %s\n' "$logs"
export PYTHONPATH="$root:$previous:$ws" PYTHONDONTWRITEBYTECODE=1
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
 -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing
 -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
 -I"$root" -I"$previous" -I"$control" -I"$ws" -I"$repo/spin/src/kernels" -I"$repo/spin/include")
objects=("$previous/Gather.o" "$previous/RoutePacked.o" "$control/PackedMixer.o" "$control/BchAvx512.o")
python3 -B "$root/further_control_codegen.py" "$previous/creative_probe.cpp" --kind align > "$root/AlignedDriver.cpp"
for offset in 0 16 32 48; do
 g++ "${flags[@]}" -DSPIN_INPUT_OFFSET=$offset -c "$root/AlignedDriver.cpp" -o "$root/Aligned$offset.o"
 g++ "$root/Aligned$offset.o" "$previous/Transpose.o" "${objects[@]}" -o "$root/inner-size-align$offset"
done
for variant in reverse unroll; do
 python3 -B "$root/further_control_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" \
  --kind "$variant" > "$root/$variant.cpp" 2> "$logs/generate-$variant.txt"
 g++ "${flags[@]}" -c "$root/$variant.cpp" -o "$root/$variant.o"
 g++ "$previous/CreativeDriver.o" "$root/$variant.o" "${objects[@]}" -o "$root/inner-size-$variant"
done
sha256sum "$root"/inner-size* "$root/further_control_codegen.py" "$0" > "$logs/hashes.txt"
printf 'Build complete; no benchmark was run.\n'
