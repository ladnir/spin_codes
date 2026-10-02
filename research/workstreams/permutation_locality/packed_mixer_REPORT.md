# Packed-Domain Mixing: First Cost Screen

Later measurement: `packed_driver_REPORT.md` identifies outlined route
callbacks in this large harness. A minimal force-inlined driver measures GL32
at 6.203 ms versus 5.805 ms legacy; the retained legacy binary is 5.829 ms.
Use that comparison for the current optimized-code tradeoff. The results
below are preserved as the original exploration and code-equivalent controls.

2026-09-30. Keeping BCH data packed makes stronger local mixing practical.
An initial build measured full GL(32,2) block mixing at **6.094 ms**, versus
**6.573 ms** unchanged. A follow-up build that also relocates the exact legacy
GF16 maps measured GL32 at **6.418 ms**, versus **6.582 ms** unchanged.
The 7.29% initial improvement is therefore build-context-sensitive, not a
robust final speedup claim. The exact-code relocation gives a smaller
**2.19% improvement** in the final matched comparison, without changing the
encoder or its setup distribution. The new GL32 construction still needs
a complete distance certificate.

The earlier retained encoder's approximately 5.83 ms timing remains a separate
measurement. This new executable has a slower inner/routing context. We have
not transplanted the new mixer into that executable and do not subtract the
observed difference to manufacture a faster result.

## Follow-Up: Relocate the Exact Existing GF16 Maps

The shared route preserves the four row-lane identities. Each packet's GF16
transpose can therefore move from the streaming route callback into packed
BCH preparation, after reindexing its control by the route destination.
This is the same encoder for the same setup seed, not a new distribution.

For packet `p`, let `base=packets.bases[p]`, `tile=base/1028`, and
`column=(base%1028)/4`. Input row lane `q` contributes to output lane `r`
when bit `8*q+2*r` of `packets.controls[p]` is set. Within each packed
eight-column byte, these coefficients form diagonal 8-by-8 binary matrices.
Their GFNI action is consequently just a bytewise AND mask.

Two controls use that observation:

- Mode 6 expands the exact old controls into diagonal GFNI matrices and
  uses the existing full-GL32 packed kernel.
- Mode 7 expands the same controls into byte masks. It replaces each of
  the four GFNI mixer terms with AND, retaining the three lane shuffles
  and XORs. Both modes omit GF16 in the route callback.

The initial implementation retains the 16 MiB expanded coefficient layout
for both controls. It does not exploit the much smaller underlying GF16
description. Exact relocation works here because the shared route leaves
each packet's four rows together; it must not be blindly applied to the
independent-row route.

Four process medians per variant, seeds 1 and 17 with forward/reverse order,
101 calls per process, all in the final executable:

| Variant | Median (ms) | Process-median range (ms) | Change from unchanged |
|---|---:|---:|---:|
| Unchanged legacy GF16 route | 6.581844 | 6.554453--6.640774 | -- |
| Exact legacy maps, packed GFNI | 6.441292 | 6.361663--6.495563 | -2.14% |
| Exact legacy maps, packed AND | 6.437951 | 6.431123--6.450519 | -2.19% |
| New full GL32, no GF16 | 6.417893 | 6.379246--6.445430 | -2.49% |

The GFNI and AND relocation timings are effectively tied at this precision.
Mode 6 and mode 7 outputs match the original `joint.cpp` shared-GF16 R2
encoder bit-for-bit, including the suffix, for both seeds and three patterns
at `K=2^14` and `K=2^20`. Their timed final checksums also match mode 0.
Both pass all 131,072 physical basis vectors and the scalar/adjoint checks.

Separate phase runs put the new GL32 inner/route at 2.985--2.990 ms and its
packed BCH at 3.417--3.424 ms. The initial build had 2.669 ms and 3.418 ms.
Thus the approximately 0.32 ms change is in the inner/routing context, not
the measured packed BCH work. Its cause remains unresolved; adding modes
changed the executable, so we do not attribute it to GFNI arithmetic.

For the exact relocation, phase means were 2.960--2.970 ms inner/route and
3.417--3.476 ms BCH with GFNI; 2.946--2.954 ms and 3.543--3.582 ms with AND.
These separate means are not a decomposition of the table's medians.
One unchanged baseline phase run had a large inner/route outlier, 6.391 ms
instead of 3.876 ms for the other seed; its raw result is retained and it
was not part of the 101-call comparison table.

No additional variants or tuning sweep followed this comparison.

## Construction and Layout

The workload is the transposed rate-one-half encoder at `K=2^20`, with
128-bit XOR elements and IMT `(t,s,R)=(128,19,2)`. Four BCH rows form one
16 KiB input tile. The routed tile stride remains 1,028 elements.

In the forward direction, each tested mixer follows the BCH outer and
precedes the existing shared-coordinate permutation. Blocks contain eight
consecutive canonical BCH coordinates, not a fresh random partition.
The transpose implementation applies the corresponding transpose map after
reverse routing and before BCH transpose.

Each prepared byte holds eight coordinates of one BCH row and one element
bitplane. The four BCH rows occupy four 128-bit SIMD lanes. The variants are:

| Mode | Packed mixer | Existing GF16 packet maps | Coefficient storage at K=2^20 |
|---|---|---|---:|
| 0 | None; unchanged GFNI BCH | Retained | None used |
| 1 | Shared GL8 on dense half only | Retained | 0.5 MiB allocated |
| 2 | Shared GL8 on all 256 coordinates | Retained | 0.5 MiB |
| 3 | Independent GL8 per row and eight-coordinate block | Retained | 4 MiB |
| 4 | Independent GL32 per four-row/eight-coordinate block | Omitted | 16 MiB |
| 5 | Same GL32 implementation as mode 4 | Retained control | 16 MiB |

Mode 1 uses only half its allocated coefficient table. Mode 0 allocates an
unused 0.5 MiB table because the harness shares setup machinery; no timed
kernel reads it. Modes 3--5 repeat each 64-bit GFNI coefficient for both
halves of a 128-bit lane. More compact coefficient representations were not
tested.

The sampler draws square binary matrices and rejects singular matrices.
With independent uniform source bits, this samples uniform GL8 or GL32.
The benchmark supplies reproducible SplitMix64-based words; it does not
implement an ideal independent-random-bit source or claim cryptographic
setup randomness. Setup and coefficient expansion are outside all timings.

The full-word variants pack all 256 coordinates once. They apply the mixer,
accumulate the sparse/systematic and parity contributions into packed BCH
outputs, then unpack those outputs once. Shared and independent-row GL8 each
use one GFNI operation per prepared vector. GL32 uses four GFNI operations,
three cyclic 128-bit-lane shuffles, and three XORs per prepared vector.
All affine constants are zero.

Why mode 4 can omit GF16: an ideal uniform GL32 map sends every fixed nonzero
32-bit input to a uniform nonzero 32-bit output. Conditional on its eight-packet
support, the nonzero four-bit labels are independent and uniform. This local
property replaces the old label randomizer; it does not bound how many
canonical blocks an outer word activates. That remaining question belongs
to the distance proof. Shared GL8 and row-independent GL8 do not have the
same packet-label property, so their GF16 maps remain.

## Initial-Build Full-Encoder Timings

Each entry is the median of four process medians: seeds 1 and 17, each with
forward and reverse variant order, 101 timed calls per process. Three warmups
precede each process's samples. The range is the minimum and maximum process
median, not a confidence interval. No sample or process was discarded.

| Variant | Median (ms) | Process-median range (ms) | Change from unchanged |
|---|---:|---:|---:|
| Unchanged GFNI + GF16 | 6.573273 | 6.544114--6.600359 | -- |
| Full-word shared GL8 + GF16 | 6.856395 | 6.847869--6.962764 | +4.31% |
| Full-word row-independent GL8 + GF16 | 6.978232 | 6.966480--7.013637 | +6.16% |
| Full-word GL32, no GF16 | 6.094287 | 6.054673--6.150832 | -7.29% |

Input initialization, permutation generation, matrix sampling, buffer
allocation, scalar checks, and final checksums are excluded. Each timed call
performs the complete in-place transpose, including routing, packed mixing,
all BCH input/output layout conversions, and output stores. As in the prior
comparison harness, repeated calls consume the preceding call's buffer;
there is no timed input reset or copy.

## Initial Cost Screen

These are separate, single-seed, 31-call process medians. They are useful for
cost diagnosis, not additional confirmation of the table above. Hot tests
repeat one tile 2,048 times per timed call. They are cache-hot, not necessarily
fully L1-resident. Bulk tests process 2,048 distinct tiles with normal routing
padding. Full timings include inner and routing work.

| Variant | Hot ns / four-row tile | Bulk mixer+BCH (ms) | Full (ms) |
|---|---:|---:|---:|
| Unchanged | 801.39 | 2.686997 | 6.549431 |
| Dense-half shared GL8 | 817.57 | 2.728275 | 6.779549 |
| Full-word shared GL8 | 953.10 | 2.886419 | 6.994600 |
| Full-word independent-row GL8 | 955.61 | 3.045765 | 6.978670 |
| Full-word GL32, no GF16 | 1238.10 | 3.599535 | 6.039844 |
| Full-word GL32, GF16 retained | 1240.40 | 3.603873 | 7.686747 |

The GL32 arithmetic is not free: its bulk mixer+BCH phase costs about 0.91 ms
more than unchanged BCH in this screen. Removing the GF16 route operation
more than offsets that cost in the measured full encoder.

Separate 31-call instrumented runs recorded mean phase times, including
warmups:

| Variant | Inner + route (ms) | Packed mixer + BCH (ms) |
|---|---:|---:|
| Unchanged | 3.877692 | 2.658068 |
| Dense-half shared GL8 | 4.086421 | 2.712609 |
| Full-word shared GL8 | 4.103256 | 2.781749 |
| Full-word independent-row GL8 | 4.093371 | 2.945597 |
| Full-word GL32, no GF16 | 2.669335 | 3.418027 |
| Full-word GL32, GF16 retained | 4.148606 | 3.476943 |

These means are not an exact decomposition of the independently measured
101-call medians. Template/code-generation and cache context also affect
the inner phase: even variants retaining the identical route operation have
different phase measurements. The matched full-encoder result includes
those effects rather than attributing every difference to an instruction.

## Checks and Measurement Conditions

All six modes passed the following before or during measurement:

- All 131,072 physical input basis vectors for a four-row tile of 128-bit
  elements, for both setup seeds. Expected images come from an ordinary-layout
  scalar mixer followed by the existing production XOR BCH circuit.
- Three dense scalar-reference inputs per mode and seed.
- Local matrix/transpose inner-product checks, 128 trials per canonical block.
- Full encoder/scalar-mixer-BCH equality, including untouched suffix, on three
  patterns at `K=2^14` for both seeds.
- Full scalar-mixer-BCH equality at `K=2^20` before every timed full run.

The full reference shares the inner/routing callback with the optimized
variant. These checks independently validate the new mixer, packed BCH
composition, and suffix behavior, not a second implementation of the inner
or permutation. No new sanitizer run was performed for this cost screen.

Machine: Peach Ryzen 7950X, GCC 15.2, pinned to logical CPU 15. Compilation:
`-O3 -DNDEBUG -std=gnu++20 -mavx2 -mavx512f -mavx512vl -mavx512bw -mgfni
-mtune=znver4 -fno-strict-aliasing`, with the existing BCH and workspace-routing
configuration macros. All variants reside in one executable. The unchanged
GFNI kernel is regenerated together with the new kernels from the same
`gfni_bch.py` and BCH matrix header.

Every check/measurement batch holds all three existing benchmark locks.
No benchmark or compiler job was run concurrently with the timings.
GCC reports 8,320-byte frames for dense-half preparation and 16,512-byte frames
for full-word preparation, including the prepared arrays. The full packed
output helpers have 8-byte frames. These figures are not spill-byte counts.

## Reproduction and Artifacts

New sources, all in `research/workstreams/permutation_locality/`:

- `packed_mixer_codegen.py`: invokes the unchanged GFNI generator and adds
  the packed mixer kernels.
- `packed_mixer_perf.cpp`: setup sampler, scalar reference, basis/adjoint
  checks, and isolated/full benchmarks.
- `packed_mixer_run.sh`: reproducible build, checks, screen, and confirmation.

No production source, existing kernel, proof input, or default was changed.
The remote build is `/tmp/spin-packed-mixer-NnskDD`; its reference-header and
previous comparison directories are `/tmp/spin-joint-9n57TT` and
`/tmp/spin-bch-compare-2e2UPS`. Those latter objects supply only the production
XOR reference and the unused comparison-harness Restrict symbol. The timed
GFNI kernels are all compiled in the new `PackedMixer.o`.

Run on Peach:

```sh
bash /tmp/spin-packed-mixer-NnskDD/packed_mixer_run.sh /tmp/spin-packed-mixer-NnskDD build
bash /tmp/spin-packed-mixer-NnskDD/packed_mixer_run.sh /tmp/spin-packed-mixer-NnskDD check
bash /tmp/spin-packed-mixer-NnskDD/packed_mixer_run.sh /tmp/spin-packed-mixer-NnskDD screen
bash /tmp/spin-packed-mixer-NnskDD/packed_mixer_run.sh /tmp/spin-packed-mixer-NnskDD confirm
```

The script also accepts explicit reference-header and comparison-object
directories as its third and fourth arguments. Their preparation is recorded
in `bch_compare_run.sh` and `bch_compare_REPORT.md`.

Raw logs are locally retained under ignored `tmp/packed-mixer-results/` and
remotely under `measurements/packed-{check,screen,confirm}-*` in this build.
The initial measured executable SHA-256 is
`cd44971552bc03c94f141f88091e0724f5bee32ce5d83c8b73285c152e5145d1`.
The generated kernel object SHA-256 is
`d73788854612f24b268beac3b56f80e0e35b2221d1454ea8145d02eaa53e079a`.

For the follow-up build, the executable SHA-256 is
`3b30f5f3fd1334cdbf734ff74d3f68ed960442041754e725635f29dad614ceac`;
the generated kernel object SHA-256 is
`ebe4ef88a815da7ceb8abfb486fc199dbd6d7a6b27d35d0d227a98a694b3a3b5`.
Reproduce that comparison with the same script's `build`, `check-relocated`,
`relocated`, and `relocated-profile` phases. Its raw logs are named
`packed-check-relocated-*`, `packed-relocated-*`, and
`packed-relocated-profile-*` in the same measurement directory.

Next: preserve the exact relocation as a code-equivalent optimization candidate,
and analyze canonical eight-coordinate block activity for the GL32 design.
Before adopting either kernel, integrate it into the retained optimized
driver and remeasure with compact coefficients if useful. A compact GF(2^32) multiplier has the
same fixed-nonzero-vector distribution, but its implementation cost was not
tested here. A cross-block shear network also remains untested.
