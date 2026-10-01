#!/usr/bin/env bash
# Usage: bash inner_size_run.sh ROOT PHASE EXPONENT [MODE ...]
# PHASE: build, link, check, screen, confirm, profile, or carry.
# Optional SPIN_SIZE_COEFF_OBJECT replaces the outer object; SPIN_SIZE_VARIANT
# gives its binary a suffix. link reuses SizeDriver.o to hold inner code fixed.
# All phases hold the same three locks; no background builds or benchmarks.
# carry screens the preserved 27-mode inner-fusion binary with 61 calls.
set -euo pipefail

root=$(realpath -- "${1:?fresh directory}")
phase=${2:?phase}
exponent=${3:?message exponent}
shift 3

case "$phase" in
  build|link|check|screen|confirm|profile|carry) ;;
  *) printf 'Unknown phase: %s\n' "$phase" >&2; exit 2 ;;
esac
case "$exponent" in
  14|15|16|17|18|19|20) ;;
  *) printf 'Exponent must be in 14..20\n' >&2; exit 2 ;;
esac
if [[ $phase == build || $phase == link ]]; then
  if (( $# )); then
    printf 'Build accepts no mode arguments\n' >&2
    exit 2
  fi
else
  if (( $# == 0 )); then
    printf 'At least one explicit mode is required\n' >&2
    exit 2
  fi
  for mode in "$@"; do
    if [[ ! $mode =~ ^[0-9]+$ ]]; then
      printf 'Mode must be a nonnegative integer: %s\n' "$mode" >&2
      exit 2
    fi
  done
fi

prior=/tmp/spin-bch-tune-fN0NiP
ref=/tmp/spin-joint-9n57TT
base=/tmp/spin-bch-tune-vKm7EZ
mixer=/tmp/spin-packed-mixer-NnskDD/PackedMixer.o
bch=/tmp/spin-bch-compare-2e2UPS/BchAvx512.o

exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75

logs=$(mktemp -d "$root/$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
variant=${SPIN_SIZE_VARIANT:-}
[[ $variant =~ ^[A-Za-z0-9-]*$ ]] || exit 2
binary="$root/inner-size${variant:+-$variant}"
coeff_object=${SPIN_SIZE_COEFF_OBJECT:-$prior/PackedCoeff1.o}
if [[ $phase == carry ]]; then
  binary=/tmp/spin-inner-packet-qw5uUw/inner-fusion
fi
printf 'phase=%s\nexponent=%s\nbinary=%s\ncpu=15\n' \
  "$phase" "$exponent" "$binary" > "$logs/context.txt"
printf 'modes=%s\n' "$*" >> "$logs/context.txt"

if [[ $phase == build || $phase == link ]]; then
  flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
    -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
    -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
    -I"$root" -I"$ref/spin/src/kernels" -I"$ref/spin/include" -I"$prior" -I"$base")
  if [[ $phase == build ]]; then
    g++ "${flags[@]}" -c "$root/inner_size_probe.cpp" -o "$root/SizeDriver.o" 2>&1 | tee "$logs/compile.txt"
  fi
  command=(g++ "${flags[@]}" "$root/SizeDriver.o"
    "$coeff_object" "$mixer" "$bch" -o "$binary")
  printf '%q ' "${command[@]}" > "$logs/command.txt"
  printf '\n' >> "$logs/command.txt"
  "${command[@]}" 2>&1 | tee "$logs/build.txt"
  sha256sum "$root"/*.h "$root/inner_size_probe.cpp" "$0" "$binary" \
    "$root/SizeDriver.o" "$coeff_object" "$mixer" "$bch" > "$logs/hashes.txt"
  nm -S -C "$binary" > "$logs/symbols.txt"
  objdump -d -C "$binary" > "$logs/disassembly.txt"
  exit 0
fi

if [[ ! -x $binary ]]; then
  printf 'Missing executable: %s\n' "$binary" >&2
  exit 2
fi
sha256sum "$binary" "$0" > "$logs/hashes.txt"
if [[ $phase == check ]]; then
  for mode in "$@"; do
    for seed in 1 17; do
      taskset -c 15 "$binary" "$mode" "$exponent" "$seed" 0 \
        > "$logs/check-$mode-$exponent-$seed.txt"
      cat "$logs/check-$mode-$exponent-$seed.txt"
    done
  done
  exit 0
fi

calls=61
repeats=1
profile=0
seeds=(1 17)
if [[ $phase == confirm ]]; then
  calls=301
  repeats=2
  seeds=(1 17 43 91)
fi
if [[ $phase == profile ]]; then
  profile=1
fi
printf 'calls=%s\nrepeats=%s\nprofile=%s\nseeds=%s\n' \
  "$calls" "$repeats" "$profile" "${seeds[*]}" >> "$logs/context.txt"
modes=("$@")
for seed in "${seeds[@]}"; do
  for ((rep=0;rep<repeats;++rep)); do
    for ((slot=0;slot<${#modes[@]};++slot)); do
      index=$slot
      if ((rep)); then
        index=$((${#modes[@]}-1-slot))
      fi
      mode=${modes[$index]}
      taskset -c 15 "$binary" "$mode" "$exponent" "$seed" "$calls" "$profile" \
        > "$logs/full-$mode-$exponent-$seed-$rep.csv"
      cat "$logs/full-$mode-$exponent-$seed-$rep.csv"
    done
  done
done
