# Shallow Feistel routing experiment

Six-round table Feistel is a useful candidate, but computing every routing
address inside the encoder is not yet the fastest implementation. At K=2^18,
building a six-round route and then encoding takes about 4.05 ms with reused
scratch, versus 4.80 ms for a compact exact-shuffle control. The measured gain
over that control is approximately 16%.

This experiment is disabled by default and has no public `CodeSpec` identifier.
The ordinary library still uses its existing shuffles and descriptor version 1.
No Feistel result below carries a distance or setup-failure certificate.

## Construction

For a d-bit domain, split a coordinate into a low part L of floor(d/2) bits and
a high part R of ceil(d/2) bits. One pair of rounds updates L by XOR with F(R),
then R by XOR with G(L). The round tables have the required output widths.
The inverse undoes these operations in reverse order. This also handles odd d
without enlarging the domain or rejecting outputs.

Each row and each region receives its own round tables. Tagged random-word
streams separate row-table generation from region-table generation. The
prototype uses the baseline's batched word generator; round evaluation consists
only of lookup and XOR operations. There is no shared permutation bank, runtime
AES call, or claim that this is an encryption scheme.

Row domains have 256 coordinates. Region domains have K/128 positions. The
complete routing still sends each outer row into distinct regions, with exactly
one coordinate per row in each region. Its distribution is not the independent
uniform-permutation distribution used by the certificates.

The outer remains BCH [256,128]. IMT maps, transvection generation, and encoding
circuits are unchanged. The full-schedule experiment covers K=2^16, 2^18, 2^20.
The compact transpose experiments currently cover K=2^18 with T128S19 only.
The public library's support for other natural lengths is unchanged.

## Measured implementation choices

Ryzen 7950X, CPU 15, GCC 15.2, Release, `znver4`, AVX-512 BCH. Each call encodes
one stream of 128-bit records. The following matched comparison is at K=2^18:

| Routing and implementation | Plan setup | First encode, reused scratch | Total, reused scratch | Total, fresh scratch |
|---|---:|---:|---:|---:|
| Exact shuffle, compact materialized route | 2.854 ms | 1.899 ms | 4.797 ms | 5.053 ms |
| 4 rounds, compact materialized route | 1.401 ms | 1.905 ms | 3.305 ms | 3.483 ms |
| 6 rounds, compact materialized route | 2.073 ms | 1.975 ms | 4.051 ms | 4.153 ms |
| 8 rounds, compact materialized route | 2.907 ms | 1.921 ms | 4.831 ms | 4.871 ms |
| 6 rounds, scalar on demand | 0.044 ms | 5.888 ms | 5.933 ms | 6.069 ms |
| 6 rounds, batched on demand | 0.044 ms | 4.752 ms | 4.796 ms | 4.815 ms |
| 8 rounds, batched on demand | 0.087 ms | 5.311 ms | 5.398 ms | 5.514 ms |

Plan and first-encode columns use the reused-scratch runs. Component medians
need not sum to the median total. Reusing scratch still samples a new route and
new transvections for each call; it does not reuse the previous code.

The compact materialized plan stores one 32-bit route per coordinate and one
copy of the IMT masks: 2.03125 MiB. The six-round on-demand plan stores 296 KiB,
including its masks. Both compact implementations need one 8 MiB scratch
buffer, instead of the ordinary plan's 8 MiB buffer plus a 1 MiB tile.

The compact experiment also omits synchronous huge-page collapse during scratch
creation. Therefore its fresh-call gain over the ordinary library is not solely
a permutation-generation gain. It avoids unused forward/tiled schedules and
uses an ordinary scratch vector. The exact-shuffle control makes these same
changes and matches the ordinary encoder's output.

With the ordinary library's full execution schedules retained, the first
experiment measured 14.426 ms fresh / 7.292 ms reused for exact shuffles, and
12.914 ms fresh / 6.428 ms reused for six-round Feistel. The two experiments
show why the execution-table and workspace lifecycle must be controlled before
attributing a speedup to a new permutation family.

Scalar on-demand evaluation has dependent table loads for every address.
The batched version evaluates 16 independent addresses at a time. It reduces
six-round warm encoding from 5.849 ms to 4.645 ms, but remains slower than the
materialized route's 1.758 ms. The batch evaluator preserves the reverse inner's
existing unrolled computation and performs no encode-time allocation.

All timings are medians of three process medians. Each process constructs nine
new seeded plans, excludes two warmups, and times a first encode followed by
three warm encodes. Method order reverses in the middle repetition. Benchmarks
run serially under the shared encoder locks. Setup plus first encode excludes
input filling, PPRF expansion, communication, base OT, hashing, and the receiver's
separate choice-byte encoding. These are not end-to-end Silent OT timings.

## Structural screen

For each sampled permutation P, the screen evaluates five fixed nonzero XOR
differences a. For each a, it counts how often P(x) XOR P(x XOR a) takes each
output value, over every x. It records the largest count over these five
differences and all output values. The probes are 1, the low-half mask, the
first high-half bit, the top bit, and the all-ones difference.

Mean largest counts were:

| Domain | Uniform shuffle | 2 rounds | 4 rounds | 6 rounds | 8 rounds |
|---|---:|---:|---:|---:|---:|
| 256 positions; 1,024 samples | 8.19 | 73.88 | 14.96 | 8.43 | 8.17 |
| 2,048 positions; 256 samples | 9.91 | 285.25 | 21.63 | 9.82 | 9.87 |
| 8,192 positions; 128 samples | 11.00 | 625.00 | 25.55 | 11.08 | 11.19 |

Two rounds is a deliberately weak diagnostic control. Four rounds retains
substantially stronger XOR-difference patterns than uniform shuffles. Six and
eight rounds are close to the uniform control on this particular statistic.
This motivates starting further analysis with six or eight rounds, not adopting
four rounds solely because it has the best timing.

The screen also measures quarter-domain concentrations of two fixed 32-coordinate
supports, and their overlaps under separately generated permutations. These
supports are contiguous and evenly spaced; they are not asserted to be BCH
codewords. The screen does not search the actual SPIN code for minimum-weight
words, inspect all input differences, or estimate a 40-bit failure probability.
It does not establish that the Feistel family has adequate distance.

## Correctness and scope

Tests exhaust permutation domains for multiple seeds and check both inverse
identities. They compare materialized, scalar, and batched routing addresses,
check global bijectivity and row/region uniqueness, and compare the full-schedule
encoder with its dense reference. Every compact benchmark first checks its
encoder against the full-schedule version of the same map.

Checksums agree across materialized/scalar/batched implementations and across
fresh/reused scratch, separately for each round count. The compact exact-shuffle
control agrees with the public map. The ordinary API and descriptor-v1 known
answers also pass. Linux AVX2 and Ryzen AVX-512 checks passed; focused ASan/UBSan
checks cover the new permutation, routing, and structural-screen code. Those
focused checks link the existing release encoder; they are not a full sanitizer
rebuild of the experimental fused encoder.

## Reproduce

```sh
cmake -S spin -B out/feistel -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BUILD_EXPERIMENTS=ON -DSPIN_TUNE=znver4
cmake --build out/feistel -j4
ctest --test-dir out/feistel --output-on-failure
out/feistel/spin_feistel_screen
# Run benchmarks serially, with no competing benchmark.
out/feistel/spin_feistel_bench 6 --reuse
out/feistel/spin_feistel_direct 6 materialized --reuse
out/feistel/spin_feistel_direct 6 online --reuse
out/feistel/spin_feistel_scalar 6 online --reuse
```

Round count zero selects the exact-shuffle control. The serial comparison
scripts expect executables under COMPARISON_ROOT/build and write raw results
under that comparison directory:

```sh
bash spin/experiments/feistel/compare.sh COMPARISON_ROOT 15
bash spin/experiments/feistel/compare_batched.sh COMPARISON_ROOT 15
python3 spin/experiments/feistel/summarize.py COMPARISON_ROOT/results COMPARISON_ROOT/results-batched
```

The initial scalar measurements in this report precede the batching change.
`spin_feistel_scalar` reproduces that implementation; `spin_feistel_direct`
now selects batching for on-demand routing. The second comparison supplies
the matched table above. Raw measurement files are not committed.

## Next step

Keep six-round materialized Feistel as the current experimental candidate.
Investigate cheaper address evaluation or block-local route caching before
selecting an on-demand implementation. Separately, test actual low-weight BCH
supports and multi-row cancellations; the permutation-only screen is not enough.
The compact exact-shuffle path is also worth integrating independently because
it preserves the original map and removes much of the setup overhead.
