#!/usr/bin/env bash
# Build only; candidates reuse one already checked creative driver object.
# Usage: creative_more_run.sh ROOT CONTROL_BUILD SOURCE_REPO
set -euo pipefail
root=$(realpath -e -- "${1:?output and overlay directory}")
control=$(realpath -e -- "${2:?retained final directory}")
repo=$(realpath -e -- "${3:?matching source repository}")
workstream="$repo/research/workstreams/permutation_locality"
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75
logs=$(mktemp -d "$root/more-build-XXXXXX")
printf 'Logs: %s\n' "$logs"
export PYTHONPATH="$root:$workstream"
export PYTHONDONTWRITEBYTECODE=1
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
  -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing
  -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$root" -I"$control" -I"$repo/spin/src/kernels" -I"$repo/spin/include" -I"$workstream")
cp "$root/creative" "$root/inner-size-creative"
for lead in 1 2 4; do
  python3 -B "$root/outer_gather_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" \
    --prefetch "$lead" > "$root/GatherP$lead.cpp" 2> "$logs/generate-gather$lead.txt"
  g++ "${flags[@]}" -c "$root/GatherP$lead.cpp" -o "$root/GatherP$lead.o"
  g++ "$root/CreativeDriver.o" "$root/GatherP$lead.o" "$root/RoutePacked.o" \
    "$control/tile5layout2.o" "$control/PackedMixer.o" "$control/BchAvx512.o" \
    -o "$root/inner-size-gather$lead"
done
if [[ -f $root/outer_unpack_codegen.py ]]; then
  python3 -B "$root/outer_unpack_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" \
    > "$root/Unpack.cpp" 2> "$logs/generate-unpack.txt"
  g++ "${flags[@]}" -c "$root/Unpack.cpp" -o "$root/Unpack.o"
  g++ "$root/CreativeDriver.o" "$root/Gather.o" "$root/RoutePacked.o" \
    "$root/Unpack.o" "$control/PackedMixer.o" "$control/BchAvx512.o" \
    -o "$root/inner-size-unpack"
fi
if [[ -f $root/outer_transpose_codegen.py ]]; then
  python3 -B "$root/outer_transpose_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" \
    > "$root/Transpose.cpp" 2> "$logs/generate-transpose.txt"
  g++ "${flags[@]}" -c "$root/Transpose.cpp" -o "$root/Transpose.o"
  g++ "$root/CreativeDriver.o" "$root/Gather.o" "$root/RoutePacked.o" \
    "$root/Transpose.o" "$control/PackedMixer.o" "$control/BchAvx512.o" \
    -o "$root/inner-size-transpose"
  if [[ -f $root/creative_basis.cpp ]]; then
    g++ "${flags[@]}" -c "$root/creative_basis.cpp" -o "$root/CreativeBasis.o"
    g++ "$root/CreativeBasis.o" "$root/Transpose.o" \
      "$control/PackedMixer.o" "$control/BchAvx512.o" -o "$root/transpose-basis"
  fi
fi
sha256sum "$0" "$root/CreativeDriver.o" "$root/RoutePacked.o" "$root/Gather.o" \
  "$control/tile5layout2.o" "$root"/GatherP* "$root"/inner-size-* > "$logs/hashes.txt"
printf 'Built further candidates; no benchmarks run.\n'
