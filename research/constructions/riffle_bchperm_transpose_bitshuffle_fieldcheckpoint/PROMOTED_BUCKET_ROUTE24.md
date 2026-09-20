# Promoted ExactPerm-BucketRoute24 evaluator

`ExactPerm-BucketRoute24` is the promoted evaluator for
`Riffle ExactPerm FieldCheckpoint v1`.  It computes the same linear map and
samples the same exact permutation family as the frozen evaluator.

The reusable implementation is
`libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h`.  Initialization
samples the structured exact permutation, constructs its inverse route, and
precomputes two immutable schedules:

- a 24-bit source-to-bucket-slot schedule occupying 6 MiB;
- a 24-bit bucket-slot-to-tile-offset schedule occupying 6 MiB.

The production initializer consumes a caller-owned cryptographic PRNG.  The
benchmark also retains an explicit deterministic initializer so receipts can
reproduce the historical permutation and coefficients exactly.

The selected hot path walks both schedules with three-byte pointer increments.
The reverse checkpoint pass writes each finalized 128-bit value to one of four
sequential bucket streams.  It then scatters one 8 MiB bucket into outer order
and applies the outer transpose two blocks at a time.  Software prefetch is
disabled.  The caller owns a reusable 40 MiB workspace: 32 MiB for routed
values and 8 MiB for one destination tile.  No allocation occurs in the hot
call.

The outer circuit is the seed-114 exact XOR DAG found by a 128-seed search.  It
uses 2934 XORs, down from 2961.  The generator rows and the `[256,128,>=38]`
code are unchanged.

On an uncontended Ryzen 9 7950X core, 21 alternating trials measured:

- frozen tiled evaluator: 21.495972 ms median;
- promoted production evaluator: 11.333555 ms median;
- production-to-baseline ratio: 0.527241.

The complete production output was compared with the frozen evaluator before
timing.  The check passed.  Repeated production-only tuning runs reach about
10.9 ms; 11.33 ms is the conservative alternating headline.

The following exact-map alternatives were tested and rejected:

- direct random scatter: about 17.2 ms;
- one-byte bucket IDs with runtime counters: slower than precomputed slots;
- paired 21-bit and 19-bit schedule decoding: slower than scalar 24-bit loads;
- explicit tile write prefetch: slower than no prefetch;
- four-way AVX-512 outer evaluation: slightly slower than two-way AVX2.

Transparent huge-page collapse for the workspace was mildly positive but
noisy.  It remains an optional deployment policy rather than part of the
portable evaluator.

The full measurement is in
`receipts/promoted_bucket_route24_peach_7950x.json`.

## 2026-08-27 phase and counter audit

An instrumented copy of the promoted path measured the following 21-trial
medians on Peach CPU 15.  The timers were placed inside one encoding call, so
the phases retain the production cache handoff:

- checkpoint state update plus four sequential value streams: 3.999 ms;
- permutation scatter within each 8 MiB tile: 3.208 ms;
- two-way AVX2 outer transpose: 3.981 ms;
- complete instrumented call: 11.255 ms.

The first value write is therefore not a random 32 MiB scatter.  It is four
interleaved monotone streams.  The tile write is the only permutation-ordered
store in the promoted evaluator.

On this AMD processor, the generic `cache-misses` event counts L2 misses.
Production-only runs measured about 2.27 million such misses per encoding.
Instruction-address sampling assigned about 81 percent of the L2-miss samples
to the tile-scatter loop.  These misses need not leave the CCD's 32 MiB L3.
The tile scatter accounts for only about 29 percent of elapsed time despite
causing most L2 misses.

Long production-only runs retired about 2.23 instructions per cycle and had a
0.33 percent branch-miss rate.  Cycle sampling assigned approximately 33
percent to the outer XOR circuit and 20 percent to the field checkpoint
transform.  The remaining production work is primarily the two routing loops.
The generated outer and field-transform functions occupy about 53 KiB and 12
KiB of machine code.  L1 instruction-cache misses were 3.16 percent; dispatch
accounting attributed about 25 percent of unused slots to the front end and 38
percent to backend stalls.

Two code-generation probes did not justify replacing the compiler output with
hand-written assembly.  Size-oriented compilation reduced the field transform
from about 12 KiB to 2.4 KiB but slowed the profiled call from 11.26 ms to
11.37 ms.  A correct seed-114 four-way AVX-512 circuit measured 11.34 ms in a
clean production-only run, versus 11.05 ms and 11.15 ms in the surrounding
AVX2 runs.  It remains available only behind
`LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER`.

The audit points to algorithmic code generation rather than handwritten
instruction spelling: reduce the outer circuit's register pressure and live
ranges, or reduce/fuse the 8 MiB tile scatter.  Straight translation of the
current loops into assembly is unlikely to move the end-to-end result much.

### Liveness-aware circuit search

The follow-up search tested both reverse-topological scheduling of a fixed XOR
DAG and different Paar circuit seeds.  For seed 114, the historical descending
schedule has 681 peak live adjoints and live area 435945.  A greedy release
schedule reduced these to 447 and 308323 without changing the 2934 XOR count.
GCC's stack frame fell from 23584 bytes to 19776 bytes, but production-only
timing worsened to 11.274 ms versus bracketing seed-114 runs of 11.123 ms and
11.193 ms.  Reducing abstract liveness and spills was therefore insufficient.

A complete 128-seed sweep compared transpose XOR count, natural-order peak
liveness, live area, and greedy-schedule liveness.  Three natural-order
candidates were checked end to end:

- seed 104: 2936 XORs, 671 peak live, 11.465 ms;
- seed 57: 2949 XORs, 672 peak live, 10.994 ms;
- seed 112: 2969 XORs, 656 peak live, 11.167 ms;
- seed 114 controls around the run: 11.052 ms and 11.043 ms.

The small apparent seed-57 win did not survive an alternating in-process
outer-only test.  Seed 114 measured 4.840 ms and seed 57 measured 4.902 ms, so
seed 57 was 1.28 percent slower.  All candidates were compared against the
frozen evaluator and produced identical complete output.

The natural-order candidates also produced larger compiler stack frames than
seed 114 despite their lower abstract peak-live counts.  Actual GCC allocation
is affected by circuit topology and instruction scheduling in ways the simple
interval metric misses.  Seed 114 remains production.  The reusable search is
in `scripts/search_bch256_transpose_pressure.py`; the paired outer benchmark is
available behind `LIBOTE_RIFFLE_BCH_CIRCUIT_PAIR_BENCH`.
