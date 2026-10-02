# First complete certificate for four-bit independent-row packets

At K=2^20 and N=2^21, the four-bit independent-row construction has
minimum distance at least 10486, hence relative distance greater than
0.5%, except with setup probability below 2^-43.744.

This is the first full-occupancy certificate for this four-bit route.
It is a starting point for increasing the proved distance, not a proposed
replacement for the stronger two-bit operating point. The encoder has
not changed: BCH[256,128], four adjacent independently shuffled rows per
group, independent regional packet permutations and packet-lane shuffles,
and two-update IMT(128,19) with the authenticated expansion and feedback
maps. The state starts at zero and is not flushed. Setup is sampled once;
the distance guarantee holds simultaneously for all nonzero messages.

## Complete coverage

Let q count groups containing at least one nonzero message row. There are
2048 groups. Every nonzero message belongs to exactly one q in 1,...,2048.

For q=1,...,58, reuse the [existing complete sparse certificates](../LOW_OCCUPANCIES.md).
They give an expected bad-message count below 6.786362e-14 at output-weight
cutoff 209715. This also bounds the smaller bad event at cutoff 10485.

For q=59,...,2048, the new positive row-mixture comparison and full
composition cover give an upper below 2^-103.426. Its 812 terminal cells
include 808 outward-bounded cells and four exactly excluded cells. There
are no unresolved cells. The cover includes all component assignments,
including mixed group types and all-one components. It does not infer
unexamined occupancies from selected examples.

Adding the two ranges gives the conservative upper

    E[number of nonzero messages with output weight <= 10485]
        < 6.786362000000000007339e-14
        < 2^-43.744.

Markov's inequality bounds the probability that this count is nonzero.
Outside that event, minimum distance is at least 10486; and
10486 > 0.005 * 2097152. The total bound is limited by the reused sparse
upper, not the new dense contribution.

The dense cover was generated with 192-bit Arb arithmetic and then
independently replayed and assembled at 256-bit precision. The latter
reported a total-margin enclosure centered at
43.74435493968343722904693849160856205243 bits. This is a bound on the
proved margin, not an estimate of the actual setup failure rate.

## What made the dense proof work

The [proof interface and implementation](README.md) retain the packet
weight mixture rather than taking the worst packet shape everywhere.
There are 70 four-row component types but only four independent mean
packet-weight coordinates to cover.

A single comparison tilt was inadequate: one that made light packets
tractable could badly overcharge all-one packets. Cell-specific positive
tilts, justified by rational affine weight majorants, resolved this
problem. The final cover has 507 such witnesses, 24 density-tangent
witnesses, two coupled affine-moment witnesses, and 275 basic witnesses.
These counts describe proof witnesses, not alternative code designs.

The run also verifies its supporting halfspaces and can enumerate exact
integer compositions near a boundary. No exact-composition witness was
needed in the final cover. The complete split tree and fresh outward
replay, rather than floating search scores, establish coverage.

## Reproduction

The retained local witness file is
`tmp/four-bit-dense-adaptive-half-final.json` (325380 bytes), SHA-256

    0a99fc6157c793e0689ef439d2560484d334cf4b8d79e91adcd08dec486944e6

It remains experiment data under `tmp/`, not a production dependency.
From the repository root:

```sh
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/assemble.py tmp/four-bit-dense-adaptive-half-final.json --precision 256
```

The assembly driver recomputes the dense bound; it relies on the earlier
sparse theorem rather than rerunning those 58 support covers. To generate
a new complete cover from scratch, with no saved atlas input:

```sh
python -B research/workstreams/permutation_locality/independent_rows/dense_closure/atlas.py --threshold 10485 --input-tilt 1/32 1/64 1/128 1/256 --max-cells 10000 --max-depth 64 --output tmp/four-bit-dense-half-regenerated.json
```

Search refinements may produce a different valid partition. Generation
is not a certificate unless its final unresolved list is empty and all
retained cells were verified outward. The assembly command fails closed
on incomplete or screening-only data.

## Hill-climbing

Keep this certificate as a fallback. Increase the cutoff with
`atlas.py --retarget`, which rechecks all old cells at the new distance.
It retains passing cells, retunes existing comparison witnesses where
possible, and subdivides the remaining cells. A new cutoff is not closed
until its whole partition passes and the total is independently replayed.
The first higher target is 1%; the [initial attempt](HILL_CLIMB.md) made
progress but did not complete. The stronger two-bit certificate remains
unchanged throughout this work.
