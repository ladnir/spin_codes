#!/usr/bin/env bash
# Usage: bash raw_outer_run.sh REPO ROOT build|check|basis [LEAF ...]
# Default leaf2; optional leaves4/8. No timing phase is provided here.
# Outputs ROOT/inner-size-rawLEAF so inner_size_compare.sh can select raw2:0
# (exact control) and raw2:1 (different raw polynomial message coordinates).
# SPIN_RAW_SOURCE selects a new-source overlay (defaults to this script's dir).
# Configurable retained inputs: SPIN_RAW_REFERENCE_DIR, SPIN_RAW_MIXER_OBJECT,
# SPIN_RAW_BCH_OBJECT, SPIN_RAW_KERNELS, SPIN_RAW_PUBLIC_INCLUDE.
set -euo pipefail
if (( $# < 3 )); then
  printf 'Usage: raw_outer_run.sh REPO ROOT build|check|basis [LEAF ...]\n' >&2
  exit 2
fi
repo=$(realpath -e -- "$1")
root=$(realpath -e -- "$2")
phase=$3
shift 3
case "$phase" in build|check|basis) ;; *) exit 2 ;; esac
leaves=("$@")
if (( ${#leaves[@]} == 0 )); then leaves=(2); fi
declare -A seen=()
for leaf in "${leaves[@]}"; do
  case "$leaf" in 1|2|4|8|16) ;; *) printf 'Invalid leaf: %s\n' "$leaf" >&2; exit 2 ;; esac
  [[ -z ${seen[$leaf]+present} ]] || exit 2
  seen[$leaf]=1
done
[[ -d $repo && -d $root ]] || exit 2
workstream="$repo/research/workstreams/permutation_locality"
overlay=$(realpath -e -- "${SPIN_RAW_SOURCE:-$(dirname -- "${BASH_SOURCE[0]}")}")
references=$(realpath -e -- "${SPIN_RAW_REFERENCE_DIR:-$root}")
kernels=$(realpath -e -- "${SPIN_RAW_KERNELS:-$repo/spin/src/kernels}")
public_include=$(realpath -e -- "${SPIN_RAW_PUBLIC_INCLUDE:-$repo/spin/include}")
mixer=${SPIN_RAW_MIXER_OBJECT:-$references/PackedMixer.o}
bch=${SPIN_RAW_BCH_OBJECT:-$references/BchAvx512.o}
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75
logs=$(mktemp -d "$root/raw-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
printf 'repo=%s\nroot=%s\noverlay=%s\nreferences=%s\nphase=%s\nleaves=%s\n' \
  "$repo" "$root" "$overlay" "$references" "$phase" "${leaves[*]}" > "$logs/context.txt"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$overlay:$workstream${PYTHONPATH:+:$PYTHONPATH}"
run() {
  printf '%q ' "$@" >> "$logs/commands.txt"
  printf '\n' >> "$logs/commands.txt"
  "$@"
}
if [[ $phase == build ]]; then
  [[ -f $references/PacketReference.h && -f $references/T64Reference.h && -f $mixer && -f $bch ]] || {
    printf 'Need prepared reference headers and retained mixer/BCH objects; set SPIN_RAW_* paths.\n' >&2
    exit 2
  }
  [[ ! -e $root/RawSizeDriver.o ]] || { printf 'Refusing to replace existing RawSizeDriver.o\n' >&2; exit 2; }
  for leaf in "${leaves[@]}"; do
    [[ ! -e $root/inner-size-raw$leaf ]] || exit 2
  done
  flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
    -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
    -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
    -I"$root" -I"$overlay" -I"$references" -I"$kernels" -I"$public_include" -I"$workstream")
  g++ --version > "$logs/compiler.txt"
  python3 --version > "$logs/python.txt" 2>&1
  run g++ "${flags[@]}" -MMD -MF "$root/RawSizeDriver.d" \
    -c "$overlay/raw_outer_probe.cpp" -o "$root/RawSizeDriver.o" 2>&1 | tee "$logs/driver.txt"
  for leaf in "${leaves[@]}"; do
    run python3 -B "$overlay/outer_composed_codegen.py" "$kernels/generated/BchCircuit.h" \
      --basis raw --leaf-bytes "$leaf" > "$root/RawOuter$leaf.cpp" 2> "$logs/generator-$leaf.txt"
    cat "$logs/generator-$leaf.txt"
    run g++ "${flags[@]}" -MMD -MF "$root/RawOuter$leaf.d" \
      -c "$root/RawOuter$leaf.cpp" -o "$root/RawOuter$leaf.o" 2>&1 | tee "$logs/outer-$leaf.txt"
    run g++ "${flags[@]}" "$root/RawSizeDriver.o" "$root/RawOuter$leaf.o" \
      "$mixer" "$bch" -o "$root/inner-size-raw$leaf" 2>&1 | tee "$logs/link-$leaf.txt"
    nm -S -C "$root/inner-size-raw$leaf" > "$logs/symbols-$leaf.txt"
    objdump -d -C "$root/inner-size-raw$leaf" > "$logs/disassembly-$leaf.txt"
    sha256sum "$root/RawOuter$leaf.cpp" "$root/RawOuter$leaf.o" "$root/inner-size-raw$leaf" >> "$logs/hashes.txt"
  done
  sha256sum "$0" "$overlay/raw_outer_probe.cpp" "$overlay/outer_composed_codegen.py" \
    "$workstream/outer_layout_codegen.py" "$workstream/packed_coeff_codegen.py" \
    "$workstream/packed_bch_tune_codegen.py" "$workstream/packed_mixer_codegen.py" "$workstream/gfni_bch.py" \
    "$workstream/inner_size_probe.cpp" "$workstream/PacketInnerKernel.h" \
    "$references/PacketReference.h" "$references/T64Reference.h" "$kernels/generated/BchCircuit.h" \
    "$root/RawSizeDriver.o" "$mixer" "$bch" >> "$logs/hashes.txt"
  printf 'Built raw-basis candidates and exact controls; no checks or timings run.\n'
  exit 0
fi
for leaf in "${leaves[@]}"; do
  binary="$root/inner-size-raw$leaf"
  [[ -x $binary ]] || exit 2
  sha256sum "$binary" >> "$logs/hashes.txt"
  if [[ $phase == basis ]]; then
    for seed in 1 17; do
      run taskset -c 15 "$binary" basis "$seed" > "$logs/basis-$leaf-$seed.txt"
      cat "$logs/basis-$leaf-$seed.txt"
    done
  else
    for exponent in 16 18; do
      for seed in 1 17; do
        for mode in 0 1; do
          run taskset -c 15 "$binary" "$mode" "$exponent" "$seed" 0 > "$logs/check-$leaf-$mode-$exponent-$seed.txt"
          cat "$logs/check-$leaf-$mode-$exponent-$seed.txt"
        done
      done
    done
  fi
done
