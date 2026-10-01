#!/usr/bin/env bash
# Usage: bash inner_size_prepare.sh REPO ROOT
# REPO is a clean checkout or archive of committed sources. ROOT must already
# exist and be empty. Builds inner-size, inner-size-tile3, inner-size-tile5.
# This script does not benchmark, delete files, or use retained build objects.
set -euo pipefail

if (( $# != 2 )); then
  printf 'Usage: inner_size_prepare.sh REPO ROOT\n' >&2
  exit 2
fi
repo=$(realpath -e -- "$1")
root=$(realpath -e -- "$2")
if [[ ! -d $repo || ! -d $root ]]; then
  printf 'REPO and ROOT must be existing directories\n' >&2
  exit 2
fi
workstream="$repo/research/workstreams/permutation_locality"
kernels="$repo/spin/src/kernels"
public_include="$repo/spin/include"
bch_header="$kernels/generated/BchCircuit.h"
bch_source="$kernels/generated/BchAvx512.cpp"
map_source="$repo/research/workstreams/rate_quarter_bch/inner_calibration/maps/t64_s16_selected.json"

for source in "$workstream/inner_size_probe.cpp" "$workstream/PacketInnerKernel.h" \
  "$workstream/InnerPacketMaps.h" "$workstream/inner_packet_codegen.py" \
  "$workstream/inner_fusion_prepare.py" "$workstream/packed_mixer_codegen.py" \
  "$workstream/packed_coeff_codegen.py" "$bch_header" "$bch_source" "$map_source"; do
  if [[ ! -f $source ]]; then
    printf 'Missing repository source: %s\n' "$source" >&2
    exit 2
  fi
done
if [[ ! -d $public_include ]]; then
  printf 'Missing include directory: %s\n' "$public_include" >&2
  exit 2
fi

exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75

if [[ -n $(find "$root" -mindepth 1 -maxdepth 1 -print -quit) ]]; then
  printf 'ROOT must be empty; nothing will be deleted: %s\n' "$root" >&2
  exit 2
fi
logs=$(mktemp -d "$root/prepare-XXXXXX")
printf 'Logs: %s\n' "$logs"
printf 'repo=%s\nroot=%s\n' "$repo" "$root" > "$logs/context.txt"
if [[ -e $repo/.git ]]; then
  git -C "$repo" rev-parse HEAD >> "$logs/context.txt"
  git -C "$repo" status --short > "$logs/git-status.txt"
else
  printf 'provenance=source archive; see source-hashes.txt\n' >> "$logs/context.txt"
fi
g++ --version > "$logs/compiler.txt"
python3 --version > "$logs/python.txt" 2>&1

# Prevent imports in the generators AND their subprocesses from writing
# __pycache__ into the source archive. Every generated output goes to ROOT.
export PYTHONDONTWRITEBYTECODE=1
find "$workstream" "$repo/spin/src" "$public_include" -type f \
  \( -name '*.py' -o -name '*.h' -o -name '*.cpp' -o -name '*.sh' \) -print0 \
  | sort -z | xargs -0 -r sha256sum > "$logs/source-hashes.txt"
sha256sum "$map_source" "$0" >> "$logs/source-hashes.txt"

run_logged() {
  local label=$1
  shift
  printf '%q ' "$@" >> "$logs/commands.txt"
  printf '\n' >> "$logs/commands.txt"
  "$@" 2>&1 | tee "$logs/$label.txt"
}

generate_stdout() {
  local label=$1
  local destination=$2
  shift 2
  printf '%q ' "$@" >> "$logs/commands.txt"
  printf '> %q\n' "$destination" >> "$logs/commands.txt"
  if "$@" > "$destination" 2> "$logs/$label.txt"; then
    cat "$logs/$label.txt"
  else
    local status=$?
    cat "$logs/$label.txt" >&2
    return "$status"
  fi
}

# Do not shadow the committed InnerPacketMaps.h in ROOT: PacketInnerKernel.h
# includes the same header from its own source directory. A second physical
# copy with only #pragma once can produce duplicate definitions.
run_logged packet-reference python3 -B "$workstream/inner_packet_codegen.py" \
  --output "$root/GeneratedPacketMaps.h" --reference-output "$root/T64Reference.h"
# Git archives can retain CRLF for historical generated headers. Ignore only
# that representation difference; do not ignore whitespace or code changes.
if ! diff --strip-trailing-cr -q "$root/GeneratedPacketMaps.h" "$workstream/InnerPacketMaps.h"; then
  printf 'Generated packet maps differ from committed InnerPacketMaps.h\n' >&2
  exit 1
fi
run_logged prior-reference python3 -B "$workstream/inner_fusion_prepare.py" \
  --output "$root/PacketReference.h"
generate_stdout packed-mixer "$root/PackedMixer.cpp" \
  python3 -B "$workstream/packed_mixer_codegen.py" "$bch_header"
for tile in 1 3 5; do
  generate_stdout "packed-coeff-$tile" "$root/PackedCoeff$tile.cpp" \
    python3 -B "$workstream/packed_coeff_codegen.py" "$bch_header" --tile-mode "$tile"
done

flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
  -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
  -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$root" -I"$kernels" -I"$public_include" -I"$workstream")
run_logged compile-bch g++ "${flags[@]}" -MMD -MF "$root/BchAvx512.d" \
  -c "$bch_source" -o "$root/BchAvx512.o"
run_logged compile-mixer g++ "${flags[@]}" -MMD -MF "$root/PackedMixer.d" \
  -c "$root/PackedMixer.cpp" -o "$root/PackedMixer.o"
for tile in 1 3 5; do
  run_logged "compile-coeff-$tile" g++ "${flags[@]}" \
    -MMD -MF "$root/PackedCoeff$tile.d" -c "$root/PackedCoeff$tile.cpp" \
    -o "$root/PackedCoeff$tile.o"
done
run_logged compile-driver g++ "${flags[@]}" -MMD -MF "$root/SizeDriver.d" \
  -c "$workstream/inner_size_probe.cpp" -o "$root/SizeDriver.o"

for tile in 1 3 5; do
  suffix="-tile$tile"
  if [[ $tile == 1 ]]; then
    suffix=
  fi
  binary="$root/inner-size$suffix"
  run_logged "link-$tile" g++ "${flags[@]}" "$root/SizeDriver.o" \
    "$root/PackedCoeff$tile.o" "$root/PackedMixer.o" "$root/BchAvx512.o" -o "$binary"
  nm -S -C "$binary" > "$logs/symbols-$tile.txt"
  objdump -d -C "$binary" > "$logs/disassembly-$tile.txt"
done

sha256sum "$root"/*.h "$root"/*.cpp "$root"/*.o \
  "$root/inner-size" "$root/inner-size-tile3" "$root/inner-size-tile5" \
  > "$logs/build-hashes.txt"
sha256sum -c "$logs/source-hashes.txt" > "$logs/source-recheck.txt"
printf 'Prepared base, tile3, and tile5 from repository sources; no benchmarks run.\n'
