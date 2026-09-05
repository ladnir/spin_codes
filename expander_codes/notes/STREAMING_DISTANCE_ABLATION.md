# Reduced-distance ablation for the streaming field encoder

## Question

The proof samples regular expander edges, edge labels, and convolution taps
with substantial independence.  The streaming implementation uses a signed
striped expander and one label per code coordinate.  An earlier implementation
also repeated its taps every 256 positions.  The experiment separates these
three changes.

The script `scripts/streaming_ec_distance_ablation.py` materializes reduced
generator matrices over a prime field.  It tests full rank, generator rows,
projective row pairs, contiguous shortening, and supports grouped by tap
phase.  The script compares these variants:

1. independent regular topology, edge labels, and taps;
2. signed striped topology with independent edge labels and taps;
3. signed stripes with shared column labels and independent taps;
4. the third variant with periodic taps.

These tests describe the sampled small codes.  They do not bound the distance
of the Goldilocks instance.

## Rank-based shortening at 26/13

The first experiment used `F_127`, region size 21, memory four, and code
dimensions `[546,273]`.  Three seeds gave full rank in every variant.  The
lightest projective row-pair weights were:

| seed | independent regular | striped, edge labels | striped, column labels | periodic taps |
|---:|---:|---:|---:|---:|
| 1 | 498 | 497 | 482 | 485 |
| 2 | 497 | 492 | 485 | 484 |
| 3 | 496 | 496 | 482 | 480 |

Shared labels lower these pair weights in this experiment.  If two message
rows share several expander coordinates, one projective coefficient can cancel
all shared coordinates with the same sign ratio.  Independent edge labels
usually give a different ratio at each shared coordinate.

Region size 21 exaggerates this mechanism because it is smaller than the 26
regions.  It also divides `127-1`, so every cyclic frequency is present in the
base field.  The rate-one-half Gilbert--Varshamov relative distance over
`F_127` is approximately 0.3651, or 199 coordinates at length 546.  The
shortening search found weights between 267 and 274.  It therefore found no
word near the random-code distance.

## Pair collisions at the deployed region size

The script `scripts/streaming_ec_pair_collision_audit.py` evaluates the same
shared-label mechanism without materializing the full generator.  It groups
collisions by their base displacement and sign ratio.

For the deployed region size `q=80659`, a sweep of 1,000 topology seeds gave:

| maximum edges canceled by one row pair | seeds |
|---:|---:|
| 1 | 828 |
| 2 | 172 |
| 3 or more | 0 |

Thus the large pair-weight reduction at region size 21 does not transfer
directly to the deployed size.  The implementation now resamples its compact
offset table until no slot difference repeats.  The deployed size enables this
check, which excludes the two-edge cases without changing the encoding loop.

## Exact tiny-code controls

The script `scripts/streaming_ec_exact_ablation.py` enumerates every nonzero
message for small ternary codes.  For 32 sampled `[20,10]` codes with degree
10/5, memory four, and tap period four, the mean distances were:

| variant | minimum over seeds | mean distance |
|---|---:|---:|
| independent regular | 4 | 4.46875 |
| signed topology | 2 | 4.12500 |
| shared labels | 2 | 3.81250 |
| periodic taps | 1 | 3.56250 |

When the tap period covered the complete 20-coordinate block, the last two
rows had the same distance distribution.  The periodic penalty therefore came
from repetition rather than another construction difference.

A second experiment used 128 sampled ternary `[18,9]` codes with degree 6/3,
memory two, and tap period three.  Shared labels had mean distance 3.8125.
Periodic taps lowered the mean to 3.640625 and increased the number of
distance-one samples from three to seven.

Tiny ternary codes magnify coincidences that are negligible over Goldilocks.
The repeated direction across two geometries nevertheless makes periodic taps
an unnecessary risk.

## Implementation decision

The current implementation no longer repeats convolution taps.  It regenerates
four coefficients per position from a public counter.  The paired encoder
shares each generated coefficient across both encoded values.

For the Goldilocks 26/13/m4 benchmark, this change reduced the persistent
schedule from 11.7 KiB to 3.7 KiB.  The receiver-pair time increased from
0.136 seconds to 0.154 seconds.  This 13 percent increase remains within the
20 percent target.

The remaining heuristic risk is the signed striped topology with shared column
labels.  Collision-free slot differences remove the direct projective-pair
mechanism.  A larger information-set or meet-in-the-middle search is necessary
to probe weights near the random-code distance.

## Meet-in-the-middle search

The script `scripts/streaming_ec_mitm_search.py` splits the message coordinates
in half and enumerates each half.  In each trial, it selects output coordinates
that must cancel, joins the two tables on those coordinates, and measures every
joined codeword.  The method is a probabilistic finder.  A reported word is an
upper bound on the sampled code's distance; failure to find a word is not a
lower bound.

An initial `[40,20]` experiment used degree 10/5 and region size four.  This
region is too small to give every pair of slots a distinct difference across
ten regions.  Its mean found weights over eight seeds were 8.625 for the
independent regular construction and 7.125 for signed stripes with shared
labels.  This experiment is a stress test of a deliberately compressed cyclic
geometry, not a model of the deployed topology.

The corrected `[40,20]` experiment used degree 4/2 and region size ten, where
the collision-free condition is feasible.  Over eight seeds, the mean found
weights were:

| variant | mean found weight | rank-deficient seeds |
|---|---:|---:|
| independent regular, edge labels | 6.625 | 0 |
| collision-free stripes, edge labels | 7.375 | 0 |
| collision-free stripes, shared labels | 5.750 | 2 |
| collision-free stripes, shared labels and edge-varying signs | 7.500 | 0 |

The two zero weights in the shared-label row were genuine rank defects.  Both
the independent-tap and periodic-tap variants had rank 19 rather than 20 for
seeds 2 and 8.  Convolution is invertible, so changing its taps cannot repair
this defect.  The edge-varying-sign experiment had full rank in all 128 seeds
tested at this geometry.  Over `F_3`, a random sign supplies every possible
nonzero edge coefficient, up to row and column scaling.

The low-degree defect did not persist when the number of regions increased.
For `[42,21]` codes with degree 6/3 and region size seven, four seeds gave mean
found weights 9.0 for the independent regular construction, 9.0 for
collision-free stripes with edge labels, 9.25 for shared labels with the
current base-constant signs, and 8.25 with edge-varying signs.  Every matrix had
full row rank.  Rank-only sweeps also found full row rank in all 128 degree-6/3
seeds at region size 13 and all 128 degree-26/13 seeds at region size 29.

These observations do not support adding a generated sign bit to every edge of
the deployed 26/13 encoder.  The extra signs repair a real low-degree failure,
but they did not improve the degree-6/3 distance sample.  The implementation
therefore retains its faster region-slot signs.  Degree-4/2 instances remain
outside the evidence supporting that heuristic.
