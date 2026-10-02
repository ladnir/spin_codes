#!/usr/bin/env bash
# Isolated serial builds against the authenticated small-size control campaign.
# Usage: bash creative_outer_run.sh ROOT CONTROL_BUILD SOURCE_REPO
set -euo pipefail
root=$(realpath -e -- "${1:?fresh output directory}")
control=$(realpath -e -- "${2:?retained final build directory}")
repo=$(realpath -e -- "${3:?matching source repository}")
workstream="$repo/research/workstreams/permutation_locality"
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75
logs=$(mktemp -d "$root/build-XXXXXX")
printf 'Logs: %s\n' "$logs"
export PYTHONPATH="$workstream"
export PYTHONDONTWRITEBYTECODE=1
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
  -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing
  -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$root" -I"$control" -I"$repo/spin/src/kernels" -I"$repo/spin/include" -I"$workstream")
cp "$control/inner-size-tile5layout2" "$root/inner-size"
for planes in 1 2 4 8; do
  python3 -B "$root/outer_plane_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" \
    --planes "$planes" > "$root/Plane$planes.cpp" 2> "$logs/generate-$planes.txt"
  g++ "${flags[@]}" -c "$root/Plane$planes.cpp" -o "$root/Plane$planes.o"
  g++ "${flags[@]}" "$control/SizeDriver.o" "$root/Plane$planes.o" \
    "$control/PackedMixer.o" "$control/BchAvx512.o" -o "$root/inner-size-plane$planes"
done
sha256sum "$root/outer_plane_codegen.py" "$0" "$root"/Plane*.cpp "$root"/Plane*.o \
  "$control/SizeDriver.o" "$root"/inner-size* > "$logs/hashes.txt"
printf 'Built plane candidates, no benchmarks run.\n'
