# GFNI Inner and the 16-State Candidate

2026-09-30. Keeping the existing 19-state IMT update in a packed GFNI
representation reduces the complete precomputed transposed encoder from
**5.889024 ms to 5.334629 ms**, a **9.41% time reduction** in a matched
comparison. The encoder and its ideal setup distribution are unchanged.
The [greater-than-10% distance / greater-than-52.05-bit certificate](packed_mixing/R4_CLOSURE.md)
therefore still applies. No production defaults or paper claims change.

The workload is K=2^20, rate one-half, with 128-bit XOR elements on Peach's
Ryzen 7950X. Setup, allocation, correctness checks, and checksums are outside
timing. The BCH/GL32 phase retains the previous compact coefficients and
wide stores. These results do not benchmark a complete 16-state encoder.

## Exact-Map Implementation

The four sampled transvections compose to an identity plus a rank-at-most-four
correction. The transposed update first computes four parities of the
19 state coordinates, then adds selected combinations to those coordinates.
The new kernel evaluates these two binary maps with zero-offset GFNI affine
instructions. It does not replace the sampled product by another random map.

Packing collects the same payload bit from eight code coordinates into one
byte. It does not mix independent payload bits. Three groups of eight
coordinate slots hold the 19-coordinate state; the five padding bits stay
zero. Each group uses two 512-bit registers to cover all 128 payload bits.

The selected mode keeps this packed state between epochs. It still unpacks
state for the original sparse output circuit and packs the syndrome from
the original zeta/finish circuit. Thus this is not an entirely packed inner:
both conversion boundaries are included in the measured gain. A separate
round-trip mode also packs the entering state on every update.

GFNI metadata occupies 48 bytes per epoch, or 768 KiB at K=2^20, versus
36 bytes and 576 KiB for the scalar-fused metadata. Packed state occupies
384 bytes, versus 304 bytes for 19 ordinary 128-bit words. The comparison
harness allocates all metadata forms in the same order for every mode.
There are no hot-path allocations or indirect calls.

## Matched Measurements

GCC 15.2, CPU 15, the existing ISA/optimization flags plus `-mgfni` for the
driver. All processes run serially under the three shared benchmark locks.
Each process uses three warmups and 101 measured calls. Seeds 1 and 17 run
in both mode orders. No samples or processes are discarded. As in the
retained harness, inputs are not reset between measured calls.

The confirmation batch is `gfni-confirm-DaquLN`:

| Seed | Order | Scalar fusion (ms) | Persistent packed GFNI (ms) | Round-trip GFNI (ms) |
|---|---|---:|---:|---:|
| 1 | scalar, persistent, round-trip | 5.911516 | 5.324269 | 5.635601 |
| 1 | round-trip, persistent, scalar | 5.892220 | 5.333346 | 5.646801 |
| 17 | scalar, persistent, round-trip | 5.859569 | 5.335911 | 5.638326 |
| 17 | round-trip, persistent, scalar | 5.885828 | 5.346892 | 5.636523 |
| | Median of process medians | **5.889024** | **5.334629** | **5.637425** |

The saving is 0.554396 ms, or 9.4140% of this matched scalar control.
Earlier measurements near 5.85 ms used a different executable and are not
the denominator for that percentage.

A separate 31-call instrumented batch, `gfni-profile-B3zUKl`, reports:

| Implementation | Inner + routing means (ms) | Packed GL32/BCH means (ms) |
|---|---:|---:|
| Scalar fusion | 2.721409--2.726578 | 3.133043--3.180746 |
| Persistent packed GFNI | 2.184281--2.188415 | 3.137693--3.175343 |
| Round-trip GFNI | 2.492356--2.495447 | 3.156816--3.184736 |

These means include warmups and do not exactly decompose the headline
medians. They locate the gain in the inner/routing phase; they do not
measure pure state-update time or establish a microarchitectural cause.

## Checks and Provenance

The compiled checks pass for both seeds:

- All 2,048 physical basis bits for S=16 and all 2,432 for S=19, in both
  map orientations; 256 fresh mask sets; post-update feedback; padding;
  and guarded unpacking at a 16-byte-offset address.
- The optional dense 16-by-16 update against original products of four
  and eight transvections, in both orientations, including all physical
  basis inputs and 64 fresh products.
- Every epoch's actual 19-state map against the original masks.
- At K=2^14 and K=2^20, all three complete encoder modes against the
  original dense inner, forward/transpose adjoint, and independent scalar
  GL32 plus production BCH reference, with dense and sparse inputs and
  the untouched buffer suffix.

Independent source review found no blocking issue. The standalone S16
tests check a linear-update implementation, not a distance certificate or
a complete S16 encoder. No sanitizer campaign was run in this experiment.

| Artifact | SHA256 |
|---|---|
| `FusedR4Gfni.h` | `d2f340910e175ded705b74a3bb5eedbfbe78964c9be49bbf00fe3db31a51c179` |
| `gfni_r4_probe.cpp` | `33a1324e5af83f583ca29b15ed80d47df6aa5381a47f57945037397f5ddd246d` |
| `gfni-r4` executable | `3b61d526a74314e23ea14880f1dede7f4d9a212cc0ebd43b54fe7770d52230a3` |

The runner is [gfni_r4_run.sh](gfni_r4_run.sh), using the same retained
reference headers and BCH objects as the previous campaign. Its fresh
remote directory is `/tmp/spin-bch-tune-Q1wko3`. Source snapshots, executable,
hash receipts, compiler reports, and all measurement logs are retained
locally under ignored `tmp/packed-hill/gfni-r4-performance/`. Raw results
remain outside version control. Original proof and implementation sources
are unchanged.

## Why Explore 16 States?

Sixteen coordinates occupy exactly two packed groups, with no padding.
A full 16-by-16 binary matrix uses four 8-by-8 blocks: eight GFNI operations
across the 128 payload bits. This equals the GFNI count for the factored
16-to-4-to-16 update, excluding packing, broadcasts, and XORs.

Consequently, setup can precompose eight or more transvections without
increasing the full-matrix hot-path cost. A uniformly sampled invertible
16-by-16 matrix is another possible distribution: it sends each fixed
nonzero state uniformly over the 65,535 nonzero states. That eliminates
the lazy component but not the probability of cancellation against a
nonzero feedback value, which is 1/65,535. No uniform-GL16 sampler or
certificate is implemented here.

Changing the state size requires new fixed maps and a new whole-code proof.
The first expansion candidate restricts the retained expansion to its
first 16 independent generator rows. Exact enumeration gives rank 16 and
nonzero weight spectrum
`48:622, 56:13840, 64:36623, 72:13808, 80:642`.

The feedback map must not be formed by simply truncating the current
19-bit columns. That projection makes columns 20 and 110 equal, creating
a weight-two zero-feedback input. The new weight-five candidate instead
uses 128 distinct columns on 16 coordinates. It has rank 16, full rank on
every four-bit packet, and exact feedback-kernel minimum weight four.
An independently defined BCH-based feedback candidate has kernel minimum
weight six; it has not been numerically screened. Definitions and exact
spectra are in [s16_maps.py](packed_mixing/s16_maps.py). All four
[structural tests](packed_mixing/s16_test_maps.py) pass.

The initial S16 comparison-point calculations are exploratory, not a
certificate. They reuse the retained outer comparison mixture and have
neither complete occupancy coverage nor freshly authenticated outer
premises. Map checks alone do not establish the requested 10% / 40-bit
whole-code bound.

The bounded weight-five screen returned the following binary64 log2
proposals after its basic variance refinement:

| Auxiliary comparison mean | Four updates | Eight updates |
|---|---:|---:|
| 0.032 | +1180.31 | +598.58 |
| 0.104 | +4634.90 | +3868.65 |

These positive values are unsuccessful bound proposals, not negative
distance margins or low-weight counterexamples. At mean 0.104, the limited
generic regional refinement gives +4167.27 for eight updates and does not
improve the basic proposal. The screen lacks the map-specific lazy-state
and joint-return refinements used in the S19 closure, and has no outward
replay. It therefore does not support a direct comparison with that final
certificate. The next proof comparison must apply equally strong bounds
to both state sizes. Diagnostic receipts are the three local
`tmp/packed-hill/s16-weight5-*.json` files.

Next: retain the measured exact-map S19 gain, then investigate S16 with
strong precomposed mixing and its own proof bound. A separate exact-map
optimization can remove conversions at the fixed-map interfaces without
changing the certified construction.
