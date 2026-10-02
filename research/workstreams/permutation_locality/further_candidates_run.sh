#!/usr/bin/env bash
# Build-only structural follow-up. Same four directories as further_control_run.
set -euo pipefail
root=$(realpath -e -- "${1:?root}")
previous=$(realpath -e -- "${2:?previous build}")
control=$(realpath -e -- "${3:?retained control}")
repo=$(realpath -e -- "${4:?source repo}")
ws="$repo/research/workstreams/permutation_locality"
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -n 7 || exit 75
logs=$(mktemp -d "$root/candidates-build-XXXXXX")
printf 'Logs: %s\n' "$logs"
export PYTHONPATH="$root:$previous:$ws" PYTHONDONTWRITEBYTECODE=1
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
 -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
 -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
 -I"$root" -I"$previous" -I"$control" -I"$ws" -I"$repo/spin/src/kernels" -I"$repo/spin/include")
retained=("$control/PackedMixer.o" "$control/BchAvx512.o")
extras=("$previous/Gather.o" "$previous/RoutePacked.o")
for planes in 4 2; do
 python3 -B "$root/further_outer_winograd_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" \
  --planes "$planes" > "$root/Winograd$planes.cpp" 2> "$logs/generate-winograd$planes.txt"
 g++ "${flags[@]}" -c "$root/Winograd$planes.cpp" -o "$root/Winograd$planes.o"
 g++ "$previous/CreativeDriver.o" "$root/Winograd$planes.o" "${extras[@]}" "${retained[@]}" -o "$root/inner-size-winograd$planes"
done
python3 -B "$root/further_map_wide_codegen.py" --output "$root/further_mapWide.h" 2> "$logs/generate-wide.txt"
python3 -B "$root/further_control_codegen.py" "$previous/creative_probe.cpp" --kind wide-driver > "$root/WideDriver.cpp"
g++ "${flags[@]}" -c "$root/WideDriver.cpp" -o "$root/WideDriver.o"
g++ "$root/WideDriver.o" "$previous/Transpose.o" "${extras[@]}" "${retained[@]}" -o "$root/inner-size-wide"
python3 -B "$root/outer_clmul_codegen.py" "$repo/spin/src/kernels/generated/BchCircuit.h" > "$root/Clmul.cpp" 2> "$logs/generate-clmul.txt"
g++ "${flags[@]}" -mvpclmulqdq -c "$root/Clmul.cpp" -o "$root/Clmul.o"
python3 -B "$root/further_control_codegen.py" "$previous/raw_outer_probe.cpp" --kind raw-driver > "$root/RawComposedDriver.cpp"
g++ "${flags[@]}" -c "$root/RawComposedDriver.cpp" -o "$root/RawComposedDriver.o"
g++ "$root/RawComposedDriver.o" "$root/Clmul.o" "${retained[@]}" -o "$root/inner-size-clmul"
sha256sum "$root"/inner-size* "$root"/*.py "$root"/*.h "$0" > "$logs/hashes.txt"
printf 'Built structural candidates; no benchmark run.\n'
