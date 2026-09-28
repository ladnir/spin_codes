#!/usr/bin/env bash
set -euo pipefail
root=$(realpath "${1:?root containing spin and research}")
cmake -S "$root/research/workstreams/permutation_locality" -B "$root/build-sanitize" \
  -DCMAKE_BUILD_TYPE=Debug -DSPIN_BUILD_TESTS=ON \
  '-DCMAKE_CXX_FLAGS=-O1 -fsanitize=address,undefined -fno-omit-frame-pointer'
cmake --build "$root/build-sanitize" -j2
# These are correctness checks, not performance measurements. Locks also keep
# their load off another project's benchmark intervals.
exec 9>/tmp/prindal-addition-encoder-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 9 || exit 1
exec 8>/tmp/bare-spin-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 8 || exit 1
exec 7>/tmp/hypercat-benchmark.lock
flock -w "${SPIN_LOCK_WAIT:-0}" 7 || exit 1
for seed in 1 17; do
  "$root/build-sanitize/locality" 14 16 1 verify "$seed"
  for g in 8 16 32; do
    "$root/build-sanitize/locality" 14 "$g" 1 wide-stream "$seed" 1
    "$root/build-sanitize/locality" 14 "$g" 1 wide "$seed" 1
  done
done
for mode in fused-1k fused-4k fused-16k; do
  "$root/build-sanitize/locality" 14 1 0 "$mode" 17 1
done
for c in 1 4 8; do
  "$root/build-sanitize/locality" 14 4 1 state-gather 17 1 "$c"
done
"$root/build-sanitize/locality" 14 16 1 xor-verify 17
for mode in xor-wide-stream wide-stride16 wide-stride20; do
  "$root/build-sanitize/locality" 14 16 1 "$mode" 17 1
done
"$root/build-sanitize/locality" 15 32 1 split-verify 17
"$root/build-sanitize/locality" 15 32 1 split-stream 17 1
"$root/build-sanitize/locality" 16 64 1 split-stream 17 1
"$root/build-sanitize/locality" 17 128 1 split-stream 17 1
for seed in 1 17; do
  for columns in 2 4; do
    "$root/build-sanitize/locality" 14 4 1 block-verify "$seed" 1 "$columns"
  done
done
for mode in block-greedy block-release block-reverse block-restrict block-share6 block-share8 block-gfni1 block-gfni2 block-gfni-static block-gfni-ternary block-gfni-unroll block-gfni-blend; do
  "$root/build-sanitize/locality" 14 4 1 "$mode" 17 1 4
done
ctest --test-dir "$root/build-sanitize/spin" --output-on-failure
for seed in 1 17; do
  for mode in random-gfni random-gfni-repack random-gfni-prepared random-gfni-prefetch random-gfni-copy random-gfni-prepared-prefetch random-gfni-copy-prefetch random-gfni-copy-mapped random-gfni-manual-copy random-gfni-copy-mapped-prefetch random-gfni-copy-packed random-gfni-manual-packed random-gfni-vector-packed random-gfni-vector-mapped; do
    "$root/build-sanitize/locality" 14 4 1 "$mode" "$seed" 1 4
  done
  "$root/build-sanitize/locality" 14 4 1 random-gfni-copy "$seed" 1 2
  "$root/build-sanitize/locality" 14 4 1 random-gfni-vector-mapped "$seed" 1 2
  "$root/build-sanitize/locality" 14 4 1 random-gfni-vector-packed "$seed" 1 2
  for mode in random-gfni-vector-prefetch1 random-gfni-vector-prefetch2 random-gfni-vector-prefetch4 random-gfni-vector-tile1 random-gfni-vector-tile1-prefetch2 random-gfni-physical random-gfni-physical-copy random-gfni-physical-prefetch random-gfni-vector-current-prefetch; do
    "$root/build-sanitize/locality" 14 4 1 "$mode" "$seed" 1 4
  done
  "$root/build-sanitize/locality" 14 4 1 random-gfni-physical "$seed" 1 2
  "$root/build-sanitize/locality" 14 4 1 verify "$seed" 1 4
done
