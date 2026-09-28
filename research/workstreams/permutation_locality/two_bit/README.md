# Two-bit replay on the current BCH/IMT construction

This is a separate proof experiment, not an encoder default. It retains
BCH[256,128], K=2^20, N=2^21, IMT(128,19), two independent transvection
updates per step, and 256-region routing. It replaces each group of four
adjacent outer rows by a pair. The original target bad event is a nonzero message
whose output has weight at most 209715 (10% of N, rounded down).

Complete mixed-support covers now verify every q=1,...,400 active pairs.
Their
combined bound remains below 2^-49.11, dominated by the one-pair case.
At 10%, higher occupancies still require coverage. Lowering the distance
target has now closed the full range through 9.25%, with margin above
49.11 bits; see [FIRST_CLOSURE.md](FIRST_CLOSURE.md).
[COVER_STATUS.md](COVER_STATUS.md) records the reproduction commands and
the remaining gap at the original distance.

The first complete points use the same encoder parameters. The dense
verifiers record the exact output-weight threshold in every new witness.
The 9.25% cover passed independent 256-bit replay. As of 2026-09-28,
the 9.5% effort is paused in favor of joint route--inner tuning.
Sparse coverage reaches 424 at that cutoff; the q>=425 dense checkpoint
has 2733 terminal leaves and 199 unresolved cells. No complete 9.5%
certificate is claimed. The checkpoint and earlier complete witnesses
remain local under ignored `tmp/`.

## What the old two-bit proof supplies

The archived Riffle g=2 calculation really did close a complete finite
proof: K=2^20, N=2^21, 9% relative distance, and 62.7194713852 bits of
failure margin. Its local layers used BCH[128,64,22], not our current
BCH-256/IMT construction. Its packet permutation was uniform globally,
not independently within our 256 regions.

That distinction matters. Under the old global shuffle, conditioning on
the packet-weight profile a=(a0,a1,a2) leaves a uniform class of size

    Q(a) = M! / (a0! a1! a2!) * 2^a1,     M=N/2.

For a fixed valid witness, the log bound is affine in a minus log2 Q(a).
It is convex, so the same witness at the three vertices bounds an entire
triangle of profiles. An exact triangle/singleton cover then gives a
finite certificate. Separately optimized witnesses at the vertices do
not justify that convexity argument.

Our regional route is not uniform conditional on one global profile.
We therefore reuse the counting and fixed-witness discipline, not the
old numerical certificate or its global normalization. The first replay
retains group-support conditioning and exact regional placement from
the [independent-row bridge](../independent_rows/BRIDGE.md).

Historical sources are in the separate `permute_conv` checkout:
`explorations/g2_triangular_cover_math.md` and `PROOF_STATUS.md`.
The recorded complete artifacts were located and their SHA-256 hashes
checked; the historical full calculation was not rerun:

| Artifact under its `out/` directory | SHA-256 |
|---|---|
| g2_triangle_hybrid_complete.json | d45cf8ffacdabf7d74669db60eea958b01a5ee30fc72b290b4cab8ca986e2333 |
| g2_triangle_hybrid_outward_parallel.json | c03d188708ee5dbc65c586cd7d35eb5cb554ed12af47b89f77d4d032f1580875 |
| g2_triangle_hybrid_outward_cell_sum.json | e7245efa0c7ff6bf6cf25189f69501bbc26805ee8e7a27e5e7178d18ebebfd1d |

## The new two-bit ensemble

There are 8192 BCH rows. Group adjacent rows into 4096 pairs. Independently
shuffle the 256 coordinates of each row. Region r receives coordinate r
from both rows in each pair, as one two-bit packet. Independently permute
the 4096 packets in each region and independently swap the lanes of each
packet with probability 1/2. Each region still contains 8192 bits, processed
as 64 inner steps with 64 two-bit slots per step.

All these choices and both transvections per step are sampled independently
at setup and then fixed for all messages. The bounds average over this
ideal setup distribution, not a particular seed or heuristic permutation.
Independent lane swaps are part of this ensemble, not an assumption about
the existing four-bit encoder.

For a fixed message, condition on each pair's packet-weight histogram.
Its active-region set is then uniform among subsets of its union size u;
the histograms of distinct pairs depend on independent row shuffles.
Within an epoch, conditional on the assignment of labeled packets to
epochs, their slots are a uniform injection independent of the entering
state. This permits the local operator to average slots and lanes before
taking a maximum over packet weights.

For local shape (a,b), with a single-bit packets and b double-bit packets,
the exact number of possible inputs is

    D(a,b) = binomial(64,a+b) binomial(a+b,a) 2^a.

There are only 44 nonempty shapes with a+b<=8, versus 494 for four-bit
packets. We rebuilt the character and cancellation counts for these new
slots; we did not reuse four-bit counts. In particular, all 192 nonzero
one-packet feedback values are distinct, and two occupied slots cannot
produce zero feedback. Three slots can: shape (2,1) has zero-feedback
probability 1/124992. Thus the analysis must still account for cancellation.

## Local enumerator and transfer bound

For a state character, let n0,n1,n2 count two-bit slots with zero, one,
or two negative coordinate signs. Its feedback character polynomial is

    product_{h=0}^2 (1 + (2-2h)y1 + (-1)^h y2)^n_h.

Coefficient (a,b), followed by exact inverse Walsh transformation, counts
each feedback state for that shape. Signed coefficients and inverse counts
are checked against D(a,b). The production calculation has 677 distinct
character histograms. Exact forced-coordinate character calculations also
give the second overlap moments used by the existing quadratic bounds.

For a fixed expanded state word of weight v, write h for its weight inside
one two-bit slot and z=exp(-lambda). Its output-moment polynomial is

    z^v product_slots (1 + ((2-h)z + h/z)y1 + z^(2-2h)y2).

Dividing coefficient (a,b) by D(a,b) gives the exact-law tilted output
moment. The implementation evaluates this positive polynomial outward.

The state envelope has nine coordinates: zero-state mass, fresh mass,
mature mass and its separate pointwise density bound, and five uniform
expansion-weight-class bounds. The density coordinate is auxiliary and
is excluded from terminal mass. Two updates give lazy probability 1/4;
the remaining 3/4 refreshes uniformly among nonzero states before feedback.
The fresh lazy branch retains a pointwise tilt bound when its feedback-law
form must be preserved. It is not replaced by an averaged tilt there.

Through the census cutoff we use the exact input moments, feedback counts,
and quadratic overlap bounds. Above it, a conditional atom bound follows
by exposing all but one packet: the last feedback value has at most one
representation among the remaining slots and its lane masks. Odd total
input weight cannot have zero feedback because every feedback column has
odd weight. The new `conditional_atom` refinement instead leaves any
censused subset of packets unexposed, dividing its exact atom cap by the
probability that its slots avoid the exposed positions. It takes the best
valid subset bound, including the original one-packet bound. These bounds
cover every local occupancy through 64;
the census cutoff does not truncate the construction.

If T_j is the complete shape-maximized epoch operator, the regional
operator for r active pairs is

    R_r = [x^r] (sum_{j=0}^64 binomial(64,j) x^j T_j)^64
                  / binomial(4096,r).

This is exact sampling without replacement. For q active pairs and fixed
support witnesses p_g, use the Poisson-binomial mixture of R_r and raise
it to the 256th power. The support denominator for pair g is
binomial(256,u_g) p_g^u_g (1-p_g)^(256-u_g). The location factor is
binomial(4096,q), not binomial(2048,q).

Outer shell and CDF caps come from the independently permuted pair of BCH
words, using the authenticated BCH spectrum caps and its exact total
2^128 words. The shared support routines are called with `rows=2`.
If a factor rho is applied to every double-bit packet in the inner bound,
the outer measure uses its reciprocal for the pair's intersection size.
CDF differences are used only in a valid decreasing-majorant fold; they
are not treated as individual shell caps.

The separate `profiles.py` also retains both profile coordinates instead
of collapsing to union support. Let A_w be the outer spectrum and put
c_w=A_w/binomial(256,w). For u active regions, b double-bit packets, and
a=u-b single-bit packets, the averaged pair count is exactly

    binomial(256,u) binomial(u,b)
        * sum_{v=0}^a binomial(a,v) c_{b+v} c_{b+a-v}.

Here v of the single-bit packets belong to the first row. Every term is
nonnegative, so replacing A by its componentwise caps gives a valid cap
for each profile. This is not a claim to know the true BCH spectrum.
The full-size calculation produces 15109 nonzero capped entries, including
the zero pair. Summing over b agrees exactly with the independent support
enumerator at all 257 union sizes. This finer table is available for the
next proof refinement. The complete occupancy covers through q=400 still
maximize over packet shapes and use the coarser support counts. The
[profile-measure branch](PROFILE_MEASURE.md) now uses the finer table and
an independently derived iid-reference envelope for selected dense boxes.
The [row-mixture cover](MIXTURE_BOUND.md) additionally bounds mixtures of
light, central, and heavy rows through a two-dimensional projection.

## Completed results

All results below use rho=1, the authenticated current maps, and 192-bit
Arb arithmetic. Reported bounds are upper bounds on expected bad-message
counts and therefore on the corresponding existence probabilities.

For q=1, the verifier extracts every support coefficient directly from
`e_zero (R_0+x R_1)^256`, without a Cauchy support-denominator loss.
It takes the best valid output tilt separately at each support, then
folds the BCH pair counts and includes all 4096 locations. The log2 upper is

    -49.1128773825041051875572910542268819581503867673541909159.

The stated conservative margin is 49.11287738250 bits. This includes
one or both active rows in the pair, all row weights, all shuffled support
sizes, and every nonzero message on those rows. It does not include messages
spanning multiple pairs.
An independent 256-bit-precision rerun confirms the same stated margin.
The driver restores its requested Arb precision after calling the older
BCH authentication code, which changes the global arithmetic context.

The higher-occupancy replay examined q in {32,64,80,128}, with every active
pair having the same support u in {128,160,192,200,224,256}. All 24 events
passed. Examples of outward log2 uppers, rounded upward:

| Active pairs q | Common support u | Log2 upper |
|---:|---:|---:|
| 32 | 128 | -2749.871958 |
| 64 | 192 | -11067.431509 |
| 128 | 128 | -9538.947016 |
| 128 | 192 | -15750.129758 |

A second completed replay examined q in {160,256} and common support
u in {38,64,96,128,160,192,224,256}, with tilts .032, .064, and .096.
All 16 additional events passed outward replay. Examples, rounded upward:

| Active pairs q | Common support u | Log2 upper |
|---:|---:|---:|
| 160 | 38 | -5424.213648 |
| 256 | 38 | -7612.523917 |
| 256 | 128 | -10021.128519 |
| 256 | 192 | -14162.450997 |

These 40 very small selected-event bounds are not whole-code margins.
They omit unequal supports, intervening support values, and other
occupancies. The pair occupancy q is not directly comparable to the same
q in a four-row experiment: the number of possible active rows differs.

## Reproduction and next step

From the repository root:

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/two_bit -p 'test_*.py'
python -B research/workstreams/permutation_locality/two_bit/profiles.py
python -B research/workstreams/permutation_locality/two_bit/probe.py --one-group
python -B research/workstreams/permutation_locality/two_bit/probe.py --cutoff 8 --tilts .032 .048 .064 --groups 32 64 80 128 --supports 128 160 192 200 224 256 --outward
python -B research/workstreams/permutation_locality/two_bit/probe.py --cutoff 8 --tilts .032 .064 .096 --groups 160 256 --supports 38 64 96 128 160 192 224 256 --outward
```

The initial ten tests check character coefficients and inverse counts against full
small-input enumeration, all translated second moments on a small map,
the output polynomial for every small state word, regional placement,
all one-group support coefficients, CDF folding, and the exact pair-profile
enumerator. Profile tests also check monotonicity under componentwise
spectrum caps and agreement with the weighted support enumerator. The large-map run
also checks map authentication, exact shape totals, and single-packet
feedback injectivity. These checks support, but do not replace, the stated
distribution and transfer argument.

Next, improve the distance of the complete relaxed-target proof.
The existing 10% sparse bounds also apply to every smaller distance target
by inclusion of bad events. Keep the local (a,b) profile available for
finer bounds if shape maximization becomes the obstruction.
A pair's complete profile has just two degrees of freedom,
so this is a substantially smaller test case for carrying profile-resolved
bounds back to four-bit packets. The old global triangular cover cannot
be substituted for these regional conditioning obligations.

No production encoder, four-bit verifier, or proof cache is modified.
Two-bit encoding performance has not been measured in this experiment.
