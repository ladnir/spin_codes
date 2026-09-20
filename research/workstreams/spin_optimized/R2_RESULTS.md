# Two-round K16 integration

`Configuration::T64S12R2` is available in both reusable libraries.
It selects the certified `(64,12,2)` inner at `K=65536`, rate 1/2,
and relative distance 10%. The full setup-failure margin is
**49.3275773868 bits**. Other message lengths are rejected.
Defaults, K20 configurations, and the submitted artifact are unchanged.

The bidirectional library supports block forward encoding, bit-packed
forward encoding, and transposed encoding, including after compaction.
The constructor and workspace APIs are unchanged. `wideView()` remains
restricted to S19. Fixed maps and routing match the earlier prototype;
two independent transvections precede each feedback addition.

## Validation

All six tests passed in each of five builds: direct AVX-512, tiled AVX-512,
AVX2-only, AVX-512 detection disabled, and ASan/UBSan. These tests cover:

- Both K16 options, two routing seeds, and two coefficient seeds.
- Dense forward and transpose references, adjoint equality, and bit-packed
  forward agreement with block forward encoding.
- Both routing layouts, multiple tile sizes, compaction, in-place transpose,
  and preservation of the unused buffer suffix.
- Frozen output hashes from the prior one-round and two-round implementations.
- Existing configurations and rejection of unsupported K16 message lengths.

The ISA audit found no EVEX instructions outside the separately dispatched
objects. The exact-map checker matched A/B to the assembled 49.33-bit
certificate; source authentication of that certificate still passes.
These checks validate the integration without changing the mathematical bound.

An initial generator error applied the reverse-update edit inside the forward
routine. The bidirectional fallback test caught it. The corrected generator
uses separate exact anchors for the forward and reverse updates; every build
reported here was run after that correction.

## Matched transpose timings

Peach Ryzen 7950X, CPU 15, GCC 15.2, Release, `znver4` tuning, AVX-512 BCH,
and direct K16 routing. Each cell is the median of three process medians,
each with 101 timed calls and three warmups. Runs were serial, with reversed
configuration order in the middle repetition. Encoding is in place;
setup, allocation, and input copying/resetting are excluded.

| Configuration | Full margin, bits | Transpose-only, seed 1 | Seed 17 | Bidirectional transpose, seed 1 | Seed 17 |
|---|---:|---:|---:|---:|---:|
| `(128,19,1)` | 41.818 | 0.343581 ms | 0.348220 ms | 0.343842 ms | 0.347429 ms |
| `(64,12,1)` | 45.382 | 0.335867 ms | 0.331409 ms | 0.333592 ms | 0.331239 ms |
| `(64,12,2)` | **49.328** | **0.343712 ms** | **0.342369 ms** | **0.344092 ms** | **0.339704 ms** |

The two-round option retains approximately the old configuration's runtime
while gaining about 7.51 bits of proved margin. The small timing differences
do not establish a speedup. Forward paths were validated, not benchmarked.
All benchmark outputs agree across the two libraries and repetitions.

After Packed24 compaction, the two-round option retains 557,056 setup bytes
in the transpose-only library and 1,900,552 bytes in the bidirectional library.
These are increases of 16 KiB and 32 KiB over the selected one-round option.
Workspace remains 3 MiB. The larger bidirectional setup preserves schedules
needed by its forward APIs.

## Reproduction and next step

Build `workstreams/spin_optimized` with `SPIN_BIDIRECTIONAL_SOURCE` pointing
to the native SPIN source snapshot. That source is copied, never edited.
The `r2_check.sh ROOT` harness assumes that snapshot is at
`ROOT/hypercat/native/spin`; it runs the build matrix before serial timings.

```sh
build-k16-r2-integrated/spin_benchmark 16 auto 101 1 0 0 6412r2
build-k16-r2-integrated/spin_bidirectional_benchmark 16 auto 101 1 0 0 6412r2
python3 workstreams/spin_optimized/r2_summary.py workstreams/spin_optimized/measurements/r2-integration
```

Raw samples, test logs, and source/binary hashes remain ignored under
`measurements/r2-integration/`. No data files were committed.
The next useful measurement is the downstream benefit of the stronger margin,
rather than another small encoder-only speed improvement.
