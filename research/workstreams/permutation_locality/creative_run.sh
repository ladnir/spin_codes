#!/usr/bin/env bash
# Usage: bash creative_run.sh ROOT PHASE EXPONENT [MODE ...]
# PHASE build/check/screen/confirm/profile. Modes 0..8 must be explicit for runs.
# SPIN_CREATIVE_REPO: source repo (default: derive from this script).
# SPIN_CREATIVE_CONTROL: retained objects + PacketReference/T64Reference (ROOT).
# SPIN_CREATIVE_SOURCE: overlay source directory (default: this script's dir).
# Optional run overrides: SPIN_CREATIVE_CALLS, SPIN_CREATIVE_SEEDS (space list),
# SPIN_CREATIVE_REPEATS, SPIN_CREATIVE_PROFILE (0/1). Checks always use calls=0.
# All phases take all three shared locks. CPU 15, serial builds/runs only.
set -euo pipefail

root=$(realpath -e -- "${1:?existing build directory}")
phase=${2:?phase}
exponent=${3:?message exponent}
shift 3
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo=$(realpath -e -- "${SPIN_CREATIVE_REPO:-$script_dir/../../..}")
control=$(realpath -e -- "${SPIN_CREATIVE_CONTROL:-$root}")
source=$(realpath -e -- "${SPIN_CREATIVE_SOURCE:-$script_dir}")
workstream="$repo/research/workstreams/permutation_locality"
binary="$root/creative"
case "$phase" in build|check|screen|confirm|profile) ;; *) printf 'Invalid phase\n' >&2; exit 2;; esac
case "$exponent" in 14|15|16|17|18|19|20) ;; *) printf 'Exponent must be 14..20\n' >&2; exit 2;; esac
if [[ $phase == build ]]; then
  (( $# == 0 )) || { printf 'Build accepts no modes\n' >&2; exit 2; }
else
  (( $# > 0 )) || { printf 'At least one explicit mode is required\n' >&2; exit 2; }
  declare -A seen=()
  for mode in "$@"; do
    [[ $mode =~ ^[0-8]$ && -z ${seen[$mode]+set} ]] || { printf 'Invalid/duplicate mode: %s\n' "$mode" >&2; exit 2; }
    seen[$mode]=1
  done
fi

exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock
flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock
flock -n 7 || exit 75
logs=$(mktemp -d "$root/creative-$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
printf 'phase=%s\nroot=%s\nrepo=%s\ncontrol=%s\nsource=%s\nexponent=%s\ncpu=15\nmodes=%s\n' \
  "$phase" "$root" "$repo" "$control" "$source" "$exponent" "$*" > "$logs/context.txt"
sha256sum "$0" > "$logs/hashes.txt"

if [[ $phase == build ]]; then
  for required in "$source/creative_probe.cpp" "$source/InnerComposed.h" "$source/route_packed_codegen.py" \
    "$source/outer_gather_codegen.py" \
    "$workstream/inner_size_probe.cpp" "$control/PacketReference.h" "$control/T64Reference.h" \
    "$control/PackedMixer.o" "$control/BchAvx512.o" "$control/tile5layout2.o"; do
    [[ -f $required ]] || { printf 'Missing build input: %s\n' "$required" >&2; exit 2; }
  done
  export PYTHONPATH="$source:$workstream${PYTHONPATH:+:$PYTHONPATH}"
  export PYTHONDONTWRITEBYTECODE=1
  flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw
    -mavx512vbmi -mgfni -mprfchw -mtune=znver4 -fno-strict-aliasing -fstack-usage
    -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
    -I"$source" -I"$root" -I"$control" -I"$workstream" -I"$repo/spin/src/kernels" -I"$repo/spin/include")
  g++ --version > "$logs/compiler.txt"
  generate() {
    local label=$1 destination=$2
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
  run_logged() {
    local label=$1
    shift
    printf '%q ' "$@" >> "$logs/commands.txt"
    printf '\n' >> "$logs/commands.txt"
    "$@" 2>&1 | tee "$logs/$label.txt"
  }
  generate route-header "$root/RoutePacked.h" python3 -B "$source/route_packed_codegen.py" --emit header
  generate route-outer "$root/RoutePacked.cpp" python3 -B "$source/route_packed_codegen.py" \
    "$repo/spin/src/kernels/generated/BchCircuit.h" --emit outer
  generate gather-outer "$root/Gather.cpp" python3 -B "$source/outer_gather_codegen.py" \
    "$repo/spin/src/kernels/generated/BchCircuit.h"
  run_logged compile-route g++ "${flags[@]}" -MMD -MF "$root/RoutePacked.d" \
    -c "$root/RoutePacked.cpp" -o "$root/RoutePacked.o"
  run_logged compile-gather g++ "${flags[@]}" -MMD -MF "$root/Gather.d" \
    -c "$root/Gather.cpp" -o "$root/Gather.o"
  run_logged compile-driver g++ "${flags[@]}" -MMD -MF "$root/CreativeDriver.d" \
    -c "$source/creative_probe.cpp" -o "$root/CreativeDriver.o"
  run_logged link g++ "${flags[@]}" "$root/CreativeDriver.o" "$root/RoutePacked.o" "$root/Gather.o" \
    "$control/tile5layout2.o" "$control/PackedMixer.o" "$control/BchAvx512.o" -o "$binary"
  sha256sum "$source/creative_probe.cpp" "$source/InnerComposed.h" "$source/route_packed_codegen.py" \
    "$source/outer_gather_codegen.py" "$root/Gather.cpp" "$root/Gather.o" \
    "$workstream/inner_size_probe.cpp" "$control/PacketReference.h" "$control/T64Reference.h" \
    "$root/RoutePacked.h" "$root/RoutePacked.cpp" "$root/RoutePacked.o" "$root/CreativeDriver.o" \
    "$control/tile5layout2.o" "$control/PackedMixer.o" "$control/BchAvx512.o" "$binary" >> "$logs/hashes.txt"
  nm -S -C "$binary" > "$logs/symbols.txt"
  objdump -d -C "$binary" > "$logs/disassembly.txt"
  printf 'Build complete; no checks or benchmarks run.\n'
  exit 0
fi

[[ -x $binary ]] || { printf 'Missing executable: %s\n' "$binary" >&2; exit 2; }
sha256sum "$binary" >> "$logs/hashes.txt"
calls=61
repeats=1
profile=0
seed_text='1 17'
if [[ $phase == confirm ]]; then calls=301; repeats=2; seed_text='1 17 43 91'; fi
if [[ $phase == profile ]]; then profile=1; fi
calls=${SPIN_CREATIVE_CALLS:-$calls}
repeats=${SPIN_CREATIVE_REPEATS:-$repeats}
profile=${SPIN_CREATIVE_PROFILE:-$profile}
seed_text=${SPIN_CREATIVE_SEEDS:-$seed_text}
if [[ $phase == check ]]; then calls=0; repeats=1; profile=0; fi
[[ $calls =~ ^(0|[1-9][0-9]*)$ && $repeats =~ ^[1-9][0-9]*$ && $profile =~ ^[01]$ ]] \
  || { printf 'Invalid calls/repeats/profile override\n' >&2; exit 2; }
if [[ $phase != check && $calls == 0 ]]; then printf 'Use phase check for calls=0\n' >&2; exit 2; fi
seed_text=${seed_text//$'\n'/ }
read -r -a seeds <<< "$seed_text"
(( ${#seeds[@]} )) || { printf 'At least one seed required\n' >&2; exit 2; }
declare -A seed_seen=()
for seed in "${seeds[@]}"; do
  [[ $seed =~ ^(0|[1-9][0-9]*)$ && -z ${seed_seen[$seed]+set} ]] \
    || { printf 'Invalid/duplicate seed: %s\n' "$seed" >&2; exit 2; }
  seed_seen[$seed]=1
done
printf 'calls=%s\nrepeats=%s\nprofile=%s\nseeds=%s\n' \
  "$calls" "$repeats" "$profile" "${seeds[*]}" >> "$logs/context.txt"
modes=("$@")
for seed in "${seeds[@]}"; do
  for ((rep=0;rep<repeats;++rep)); do
    for ((slot=0;slot<${#modes[@]};++slot)); do
      index=$slot
      if ((rep%2)); then index=$((${#modes[@]}-1-slot)); fi
      mode=${modes[$index]}
      command=(taskset -c 15 "$binary" "$mode" "$exponent" "$seed" "$calls" "$profile")
      printf '%q ' "${command[@]}" >> "$logs/commands.txt"
      printf '\n' >> "$logs/commands.txt"
      "${command[@]}" > "$logs/mode-$mode-k$exponent-seed-$seed-rep-$rep.csv" 2>&1
      cat "$logs/mode-$mode-k$exponent-seed-$seed-rep-$rep.csv"
    done
  done
done
