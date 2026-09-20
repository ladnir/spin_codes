#!/usr/bin/env bash
set -euo pipefail
root=${1:?repository root}
cd "$root"
result=workstreams/spin_optimized/measurements/single_forward
mkdir -p "$result"
for mode in timing profile; do
  prof=OFF; [[ $mode != profile ]] || prof=ON
  build="build-single-$mode"
  cmake -S workstreams/spin_optimized -B "$build" -DCMAKE_BUILD_TYPE=Release \
    -DSPIN_BCH_AVX512=ON -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARK=ON \
    -DSPIN_FORWARD_FOUR=ON -DSPIN_FORWARD_PROFILE="$prof" \
    -DSPIN_BIDIRECTIONAL_SOURCE="$root/hypercat/native/spin" > "$result/$mode-configure.log" 2>&1
  targets=(spin_single_benchmark)
  [[ $mode != timing ]] || targets+=(spin_single_transpose_benchmark)
  cmake --build "$build" --target "${targets[@]}" -j2 > "$result/$mode-build.log" 2>&1
done
for repeat in 1 2 3; do
  directions=(forward transpose inplace)
  [[ $repeat != 2 ]] || directions=(inplace transpose forward)
  for m in 16 18 20; do
    for direction in "${directions[@]}"; do
      exe=spin_single_benchmark
      [[ $direction == forward ]] || exe=spin_single_transpose_benchmark
      "build-single-timing/$exe" "$m" "$direction" 101 > "$result/m$m-$direction-default-r$repeat.json"
    done
    for tile in 64 256 512; do
      build-single-timing/spin_single_benchmark "$m" forward 101 "$tile" > "$result/m$m-forward-t$tile-r$repeat.json"
    done
  done
done
for m in 16 18 20; do
  for tile in 0 64 256 512; do
    build-single-profile/spin_single_benchmark "$m" forward 31 "$tile" \
      > "$result/profile-m$m-t$tile.json" 2> "$result/profile-m$m-t$tile.log"
  done
done
