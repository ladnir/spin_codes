#!/usr/bin/env bash
# Optional final-neighborhood builds after inner_size_prepare.sh.
# Usage: bash inner_size_extra.sh REPO ROOT
set -euo pipefail
repo=$(realpath -e -- "${1:?source repository}")
root=$(realpath -e -- "${2:?prepared build directory}")
workstream="$repo/research/workstreams/permutation_locality"
kernels="$repo/spin/src/kernels"
[[ -f $root/SizeDriver.o && -f $root/PackedMixer.o && -f $root/BchAvx512.o ]] || exit 2
for name in inner-size-layout2 inner-size-tile5layout2 inner-size-fallback; do
  [[ ! -e $root/$name ]] || { printf 'Refusing to replace %s\n' "$root/$name" >&2; exit 2; }
done
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75
logs=$(mktemp -d "$root/extra-XXXXXX")
printf 'Logs: %s\n' "$logs"
export PYTHONDONTWRITEBYTECODE=1
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
  -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
  -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$root" -I"$kernels" -I"$repo/spin/include" -I"$workstream")
run() {
  printf '%q ' "$@" >> "$logs/commands.txt"
  printf '\n' >> "$logs/commands.txt"
  "$@"
}
for tile in 1 5; do
  variant=layout2
  [[ $tile == 1 ]] || variant=tile5layout2
  run python3 -B "$workstream/outer_layout_codegen.py" "$kernels/generated/BchCircuit.h" \
    --mode 2 --tile-mode "$tile" > "$root/$variant.cpp" 2> "$logs/generate-$variant.txt"
  run g++ "${flags[@]}" -mavx512vbmi -c "$root/$variant.cpp" -o "$root/$variant.o"
  run g++ "${flags[@]}" -mavx512vbmi "$root/SizeDriver.o" "$root/$variant.o" \
    "$root/PackedMixer.o" "$root/BchAvx512.o" -o "$root/inner-size-$variant"
done
# Compile every TU without VBMI; do not infer fallback support merely from
# selecting a legacy recipe in a translation unit allowed to emit VBMI.
run g++ "${flags[@]}" -mno-avx512vbmi -c "$workstream/inner_size_probe.cpp" -o "$root/FallbackDriver.o"
run g++ "${flags[@]}" -mno-avx512vbmi -c "$root/PackedCoeff1.cpp" -o "$root/FallbackCoeff.o"
run g++ "${flags[@]}" -mno-avx512vbmi -c "$root/PackedMixer.cpp" -o "$root/FallbackMixer.o"
run g++ "${flags[@]}" -mno-avx512vbmi -c "$kernels/generated/BchAvx512.cpp" -o "$root/FallbackBch.o"
run g++ "${flags[@]}" -mno-avx512vbmi "$root/FallbackDriver.o" "$root/FallbackCoeff.o" \
  "$root/FallbackMixer.o" "$root/FallbackBch.o" -o "$root/inner-size-fallback"
sha256sum "$0" "$workstream/outer_layout_codegen.py" "$workstream/inner_size_probe.cpp" \
  "$workstream/PacketInnerKernel.h" "$workstream/PacketInnerBoundaryTests.h" \
  "$root/SizeDriver.o" "$root"/Fallback*.o "$root/layout2.o" "$root/tile5layout2.o" \
  "$root/inner-size-layout2" "$root/inner-size-tile5layout2" "$root/inner-size-fallback" > "$logs/hashes.txt"
printf 'Built factored layouts and no-VBMI fallback; no benchmarks run.\n'
