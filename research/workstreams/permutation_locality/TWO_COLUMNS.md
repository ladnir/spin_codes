# Two-column bundles: same speed, all single-group ranks covered

2026-09-24. Use two-column bundles as the current candidate. They encode
in about 5.07--5.10 ms at K=2^20, or 1.89x faster than matched production
controls. Outward replays at 192 and 384 bits bound the combined failure
contribution from every single-active-group message by 5.442104e-15.
This gives 47.3847 bits for that class, including all 2048 group locations.

The full certificate remains open because messages spanning multiple groups
are not covered. The 2x performance target also remains open. Production
defaults and the paper are unchanged.

## Why change the bundle size now?

The earlier four-column layout gained speed by writing four contiguous
cache lines together. The implementation in PHYSICAL_COLUMNS.md instead
writes individual complete cache lines directly into BCH order.
It no longer needs four source columns to share a contiguous destination.
This makes smaller bundles worth revisiting.

Two-column bundles give each row group twice as many independently placed
windows along the recursion. The additional spreading substantially improves
the current distance bound. Matched measurements show essentially unchanged
encoding time. This is a change to the permutation distribution, not merely
an alternative implementation of the four-column candidate.

## Distribution and instance

The instance has K=2^20 message bits, N=2^21 output bits, BCH[256,128],
and IMT(128,19). Fix the 8192 outer rows into 2048 groups of four.
Each group independently samples a uniform permutation of all 256 coordinates;
its four rows share this permutation.

Partition the shuffled coordinate positions into 128 consecutive pairs.
Each pair defines a macroregion. Independently in each macroregion, uniformly
permute the 2048 row groups. Each group occupies one aligned eight-input window:
two columns containing four row bits each. Independently shuffle the four row
labels inside each column. All choices are independent except for the stated
sharing within a row group.

Each macroregion contains 128 IMT epochs. A fixed group therefore occupies
a uniform epoch and a uniform one of its sixteen eight-input windows.
These two choices are independent. No window crosses an epoch boundary.

The transpose implementation is the existing `random-gfni-physical` mode,
with its column parameter set to two. It assembles each column in registers,
restores its row order, and writes one aligned cache line. BCH reads the
canonical layout directly. No new encoder arithmetic or allocation was needed.

## Outward bound for one active group

Fix a nonzero four-row tuple before setup. Its row span has rank h in
{1,2,3,4}. Let u be the number of coordinates at which at least one row is
nonzero. The shared uniform coordinate permutation sends this support to
a uniform u-subset of the 256 shuffled positions.

There are fourteen nonempty local shapes: sorted pairs (a,b), where a and b
are column weights between zero and four. The census enumerates all 255
nonempty eight-bit masks in all sixteen epoch windows. Every such mask has
nonzero feedback under the actual selected B map. Thus every nonempty window
activates a zero state. The census also reconstructs weighted cancellation
counts and the action of the actual expansion map E.

The verifier uses the existing seven-coordinate positive envelope: zero,
an arbitrary nonzero state, and density bounds on the five E-weight classes.
Those classes have weights 48, 56, 64, 72, and 80. They are domination bounds,
not assumptions that conditioned states are uniform. Setup probabilities use
the ideal uniform permutations and transvections from this model.

Let W_0 be the empty-epoch envelope and W_s the envelope for local shape s.
For one macroregion, compute

    R_0 = W_0^128,
    R_s = (1/128) sum_{j=0}^{127} W_0^j W_s W_0^(127-j).

For each occupied-column count b in {1,2}, take an entrywise upper envelope
over the shapes with b nonzero columns. This deliberately forgets their
remaining composition. The support recurrence uses R_0 + 2x R_1 + x^2 R_2.

Condition on the first occupied macroregion having b occupied columns and
l macroregions after it. Its probability is

    binom(2,b) binom(2l,u-b) / binom(256,u).

The coefficient recurrence retains the exact support size in the tail.
For each conditioning class, take the best Chernoff bound from seven fixed
tilts, cap it at one, and then average over the first-region classes.
The bad event is output weight at most 209715, the integer threshold for
10% of N.

For each rank, authenticated BCH spectrum caps and exact shortening witnesses
bound the cumulative number of tuples with support at most u. A decreasing
majorant of the conditional failure probabilities permits summation by parts
using these CDF caps. The verifier never treats differences of CDF upper bounds
as actual shell counts. Finally it sums the four ranks and all group locations.

All transfer and probability arithmetic uses Arb, with upward rounding of reused
bounds. The verifier resets precision after importing historical authentication
modules. Independent 192-bit and 384-bit runs agree:

| Rank of the four row words | Failure upper, including 2048 locations | Margin (bits) |
|---|---:|---:|
| 1 | <5.050576e-15 | 47.4925 |
| 2 | <4.758844e-24 | 77.4757 |
| 3 | <4.118429e-22 | 71.0403 |
| 4 | <3.915270e-16 | 51.1817 |
| **Combined** | **<5.442104e-15** | **47.3847** |

The combined row is rounded directly from the verified total. It leaves more
than 9.040525e-13 of the overall 2^-40 failure budget for multiple active groups.
That remaining budget is not a bound on their contribution.

This replay does not need the new OA support bounds, total-weight buckets,
or orbit-memory refinement explored in RANK_FOUR_SCREEN.md. Its outer inputs
authenticate 163 dependency files and check 67 exact rational shortening
witnesses. Authentication does not replay every historical BCH proof.

## Matched performance

Peach Ryzen 7950X, GCC 15.2, CPU 15, Release/znver4, 128-bit elements,
mask seed 2. Each process uses three warmups and 101 timed in-place encodes
without resetting the input. Setup, allocation, initialization, full-output
validation, and checksums are excluded. All runs are serial under the three
shared benchmark locks.

Entries below are medians of three process medians. Reverse mode order in
the middle repetition.

| Encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production control | 9.607856 | 9.662007 |
| Four-column candidate | 5.080054 | 5.072349 |
| **Two-column candidate** | **5.074072** | **5.098929** |

The two-column process medians range from 5.063102 to 5.078361 ms for seed 1,
and from 5.071988 to 5.125628 ms for seed 17. The matched speedups are
1.8935x and 1.8949x. Approximately another 5.3% latency reduction is needed
for 2x. The comparison does not support claiming a speed advantage over
four columns; its advantage is the stronger proof at essentially the same speed.

Every timed encoder checks the whole in-place result against its materialized
reference, including the suffix. The two-column mode already passed ASan/UBSan
at K=2^14, seeds 1 and 17, in the campaign recorded in PHYSICAL_COLUMNS.md.
The executable is unchanged from that campaign; only runtime geometry differs
from the previous K=2^20 performance comparison.

## Reproduction and next step

From the repository root:

    python -B research/workstreams/permutation_locality/two_column_verify.py --precision 192
    python -B research/workstreams/permutation_locality/two_column_verify.py --precision 384

The local census exhausts all 4080 nonzero window/mask choices. The zero mask
maps to zero by linearity. Its overlap checks cover every eight-bit state-window pattern.
The verifier also tests exact small coefficient recurrences, first-region
counting, zero-tilt mass, and agreement with independently implemented
binary64 region transfers. The latter is a consistency check, not the basis
of the outward result.

Use `confirm_two_columns.sh ROOT` for the matched benchmark. Raw measurements
remain under `/tmp/spin-locality-JjY7yR/measurements/two-columns` on Peach.
Release executable SHA-256:
`78f60158b994f8601451fb0b10e6e0e59c1f33ba5e1263e92b0da656c932e6fd`.
Unchanged production library SHA-256:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

Next, analyze two active groups, including their possible placement in the
same IMT epoch, then extend the occupancy range. Retain two columns while
seeking the final kernel improvement. A full certificate must include those
multiple-group messages; the single-group result alone does not prove
minimum distance for the sampled code.

Follow-up: [TWO_GROUPS.md](TWO_GROUPS.md) now covers the two-group cases
where both row spans have rank at most two. Higher ranks and larger
occupancies still require bounds.
