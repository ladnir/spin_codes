#!/usr/bin/env bash
# All compilation, correctness checks, and benchmarks are strictly serial.
set -euo pipefail
root=$(realpath "${1:?fresh directory}")
phase=${2:?phase}
shift 2
prior=/tmp/spin-bch-tune-fN0NiP
ref=/tmp/spin-joint-9n57TT
base=/tmp/spin-bch-tune-vKm7EZ
objects=(/tmp/spin-packed-mixer-NnskDD/PackedMixer.o /tmp/spin-bch-compare-2e2UPS/BchAvx512.o)
exec 9>/tmp/prindal-addition-encoder-benchmark.lock; flock -n 9 || exit 75
exec 8>/tmp/bare-spin-benchmark.lock; flock -n 8 || exit 75
exec 7>/tmp/hypercat-benchmark.lock; flock -n 7 || exit 75
logs=$(mktemp -d "$root/$phase-XXXXXX")
printf 'Logs: %s\n' "$logs"
flags=(-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni -mprfchw -mtune=znver4
  -fno-strict-aliasing -DSPIN_BCH_AVX512=1 -DSPIN_WORKSPACE_ROUTING_OPT=1
  -I"$ref/spin/src/kernels" -I"$ref/spin/include" -I"$prior" -I"$base")
if [[ $phase == bootstrap ]]; then
  g++ "${flags[@]}" -fstack-usage -c "$base/gfni_t64_probe.cpp" -o "$root/Driver.o" 2>&1 | tee "$logs/driver.txt"
  g++ "${flags[@]}" -c /tmp/spin-cost-profile-LaoICZ/encoder_cost_probe.cpp -o "$root/Cost.o" 2>&1 | tee "$logs/cost.txt"
  g++ "${flags[@]}" -c "$prior/packed_coeff_r4.cpp" -o "$root/Basis.o" 2>&1 | tee "$logs/basis.txt"
  g++ "$root/Driver.o" "$prior/PackedCoeff1.o" "${objects[@]}" -o "$root/full-base"
  g++ "$root/Cost.o" "$prior/PackedCoeff1.o" "${objects[@]}" -o "$root/cost-base"
  sha256sum "$root/Driver.o" "$root/Cost.o" "$root/Basis.o" "$root/full-base" "$prior/PackedCoeff1.o" > "$logs/hashes.txt"
  exit 0
fi
if [[ $phase == build ]]; then
  family=${1:?generator family}; mode=${2:?mode}; tag=${3:?tag}
  # Optional compiler controls affect only this candidate object. The shared
  # driver and correctness harnesses remain byte-identical across variants.
  extra_flags=("${@:4}")
  PYTHONPATH="$prior" python3 -B "$root/outer_${family}_codegen.py" /tmp/spin-bch-compare-2e2UPS/BchCircuit.h --mode "$mode" > "$root/$tag.cpp" 2> "$logs/generator.txt"
  printf '%q ' g++ "${flags[@]}" "${extra_flags[@]}" -fstack-usage -c "$root/$tag.cpp" -o "$root/$tag.o" > "$logs/build-command.txt"
  printf '\n' >> "$logs/build-command.txt"
  g++ "${flags[@]}" "${extra_flags[@]}" -fstack-usage -c "$root/$tag.cpp" -o "$root/$tag.o" 2>&1 | tee "$logs/compile.txt"
  g++ "$root/Driver.o" "$root/$tag.o" "${objects[@]}" -o "$root/full-$tag"
  g++ "$root/Cost.o" "$root/$tag.o" "${objects[@]}" -o "$root/cost-$tag"
  sha256sum "$root/$tag.cpp" "$root/$tag.o" "$root/full-$tag" "$root/cost-$tag" > "$logs/hashes.txt"
  nm -S -C "$root/full-$tag" > "$root/symbols-$tag.txt"
  exit 0
fi
if [[ $phase == basis ]]; then
  tag=${1:?tag}
  g++ "$root/Basis.o" "$root/$tag.o" "${objects[@]}" -o "$root/basis-$tag"
  for seed in 1 17; do
    taskset -c 15 "$root/basis-$tag" basis "$seed" > "$logs/basis-$seed.txt"
    cat "$logs/basis-$seed.txt"
  done
  exit 0
fi
if [[ $phase == check ]]; then
  for tag in "$@"; do
    for seed in 1 17; do
      for exponent in 14 20; do
        taskset -c 15 "$root/full-$tag" "$exponent" "$seed" 0 > "$logs/check-$tag-$exponent-$seed.txt"
        cat "$logs/check-$tag-$exponent-$seed.txt"
      done
    done
  done
  exit 0
fi
if [[ $phase == cost ]]; then
  for tag in "$@"; do
    taskset -c 15 "$root/cost-$tag" 1 31 2048 > "$logs/cost-$tag.csv"
    cat "$logs/cost-$tag.csv"
  done
  exit 0
fi
[[ $phase == screen || $phase == confirm || $phase == profile ]] || exit 2
calls=31; repeats=1; profile=0
if [[ $phase == confirm ]]; then calls=101; repeats=2; fi
if [[ $phase == profile ]]; then profile=1; fi
tags=("$@")
for seed in 1 17; do
  for ((rep=0;rep<repeats;++rep)); do
    for ((slot=0;slot<${#tags[@]};++slot)); do
      index=$slot; if ((rep)); then index=$((${#tags[@]}-1-slot)); fi
      tag=${tags[$index]}
      taskset -c 15 "$root/full-$tag" 20 "$seed" "$calls" "$profile" > "$logs/full-$tag-$seed-$rep.csv"
      cat "$logs/full-$tag-$seed-$rep.csv"
    done
  done
done
