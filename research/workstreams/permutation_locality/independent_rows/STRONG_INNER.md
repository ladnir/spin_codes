# Independent rows: stronger inner bounds and outward replay

The independent-row route now has complete outward coverage for every
occupancy from 1 through 58. Their combined margin exceeds 43.74 bits;
see [LOW_OCCUPANCIES.md](LOW_OCCUPANCIES.md) for the later batch and
reproduction commands. The earlier selected runs and refinements below
record how the proof developed. These use the new averaged outer counts,
not an old-route certificate. The full-code target is still open.

## Construction and scope

The construction uses BCH[256,128], K = 2^20, N = 2^21, four adjacent rows
per packet group, independent coordinate permutations per row, and the
existing uniform region-group and packet-lane permutations. The inner is
IMT(128,19) with two independently sampled transvections per step, initial
state zero, and no flush. The bad event is output weight at most 209715
for some nonzero message in the stated occupancy class.

An occupancy q means exactly q of the 2048 four-row outer groups have a
nonzero message. A result for q alone does not cover smaller occupancies,
larger occupancies, or their union. All probabilities average over ideal
setup randomness fixed for all messages afterward.

## Reused inner bounds; new outer measure

`verify.py` builds the eleven-coordinate envelope from the existing
fresh-state, window-average, multi-window, fresh-collision, zero-return,
mature-tail, and pair-tail routines. It transforms the update coefficients
to two transvections and checks 330 direct empty/single-window component
inequalities per parameter choice. Region placement uses outward matrix
coefficient extraction without replacement. Optional switches expose the
later full-feedback, window-histogram, and joint-cancellation refinements.

The exact rational independent-row CDF supplies the outer measure. Each
entry is rounded upward to an integer for the existing CDF-fold interface;
monotonicity and domination are checked. No rank-specific shared-shuffle
count is used. Weighted all-one penalties remain coupled: the inner
receives rho^J, and the outer measure receives rho^(-J).

The universal cover helper searches using binary64 proposals. A successful
cover checks every support-vector box and its label multiplicity, then
recomputes the complete sum using outward Arb arithmetic. Proposal scores
alone are never a certificate. The original shared-route entry point and
its default construction are unchanged.

## Verified occupancy fourteen

With 192-bit precision, the full support cover uses 85 retained leaves
after 25 splits. It bounds the entire occupancy-fourteen bad-event
contribution by

    U_14 < 1.147336234632e-180 < 2^(-597.7487).

The witness grid is tilts .005, .008, .016 and penalties .75, 1. The
full-domain identity covers 219^14 ordered support vectors, corresponding
to every union size from 38 through 256 in each active group. The location
factor binomial(2048,14) and all nonzero message choices are included.
The margin is for this contribution only, not the whole code. The later
batch independently covers the other occupancies through 32.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 14 --tilts .005 .008 .016 --penalties .75 1 --max-splits 100 --target-bits 40
```

## Verified occupancy thirty-two

At 192-bit precision, a separate full-domain cover uses 529 retained
leaves after 25 splits. It includes every ordered support vector in
{38,...,256}^32 and all binomial(2048,32) active-group locations. Its
outward bad-event contribution is

    U_32 < 2.552383037069e-88 < 2^(-290.9778).

The witness grid is tilts .008, .016, .024 and penalties .75, .95. The
largest retained box allows eight groups to have support 38–147 and
twenty-four to have support 148–256, including all assignments through
the cover's multiplicity. Thus this is not merely a homogeneous-support
check. This particular run covers only occupancy 32; the later batch
in LOW_OCCUPANCIES.md fills the intervening occupancies.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 32 --tilts .008 .016 .024 --penalties .75 .95 --max-splits 100 --target-bits 40
```

## Initial 64-group screen with eleven coordinates

The base eleven-coordinate pipeline, without the later optional
refinements, was screened at tilts .032 and .048 and penalty .95.
Each row below concerns the homogeneous vector with the listed union size
in all 64 active groups. These are binary64 proposals, not outward claims:

| Union size | Best proposed log2 contribution upper |
|---|---:|
| 128 | -2507.55 |
| 144 | -1081.82 |
| 160 | +422.90 |
| 176 | +1719.58 |
| 192 | +2751.31 |
| 224 | +3649.97 |
| 256 | +1046.25 |

The positive entries identify an unresolved bound, not a bad word. The
largest remaining tested gap lies at larger supports. The outer counting
improvement therefore does not by itself complete the inner proof.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .032 .048 --penalties .95 --probe-supports 128 144 160 176 192 224 256
```

## Refined 64-group exact-support results

The later run adds full-feedback distributions through six occupied
windows, exact window-histogram moments through eight windows, and joint
cancellation/output coefficients through three windows. It uses tilts
.032 and .048 and penalties .75 and .95. The input enumerations pass their
positivity, total-count, direct-enumeration, and marginal cross-checks.

This run also uses exact-support shell caps where they improve on the
CDF. The shell caps come from the positive subset measure in `support.py`,
not differences of a CDF upper. After selecting rational coefficient
witnesses in binary64, `--shell-points` replays the complete point bound
outward at 192-bit precision. At support 256, the coefficient witness
p = 1 is evaluated directly without dividing by a vanishing probability.

For clarity, these bounds concern the expected number of bad messages
with q = 64 whose **shuffled** union has size u in every active group.
The support condition depends on setup. Markov's inequality bounds the
probability that any such message exists, but these seven support events
do not cover all support vectors at occupancy 64.

| Common union size u | Outward log2 contribution upper, rounded | Below 2^-40? |
|---|---:|---|
| 128 | -3475.292 | Yes |
| 144 | -1907.071 | Yes |
| 160 | -910.105 | Yes |
| 176 | -169.702 | Yes |
| 192 | +284.705 | No |
| 224 | -25.214 | No |
| 256 | -4400.041 | Yes |

The decimal entries are rounded summaries, not standalone directed-rounding
certificates. The verifier prints their Arb enclosures and tests the budget
against the outward value. Every selected witness uses penalty .95;
supports 128 and 144 use tilt .032, and the others use .048.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .032 .048 --penalties .75 .95 --full-feedback 6 --window-histogram 8 --joint-cancellation --shell-points --probe-supports 128 144 160 176 192 224 256
```

The dense-endpoint problem partly came from using a CDF as a shell bound.
At penalty .95, the log2 count at exactly 256 is bounded by 489.458 rather
than the CDF's 513.213. Across 64 groups, this alone removes about 1520.34
bits of counting slack. At supports 160,176,192,224 the corresponding
improvements are only 10.70,14.35,20.06,53.52 bits across all groups. Thus
the middle-support gap cannot be explained by this CDF issue alone.

## Later tuning and the remaining obstruction

Additional outward probes tested larger output tilts, different all-one
penalties, and total-input-weight tilts. None closed the middle support
range. For the homogeneous support-192 event at q = 64:

| Input-weight tilt a | Best tested output tilt / penalty in that run | Outward log2 upper, rounded |
|---|---|---:|
| 1 | .048 / .9 | +277.892 |
| 1.02 | .048 / .9 | +528.610 |
| .96 | .048 / .95 | +642.710 |
| .94 | .048 / 1 | +980.955 |

These are upper bounds for restricted events, not estimates of actual
failure probabilities or optimized limits of the method. At a = 1,
support 224 improved slightly to -28.428 bits but still missed even the
40-bit point budget. The grid above does not rule out every other tilt;
it does show that these inexpensive refinements have not removed the gap.
[WEIGHT_TILT.md](WEIGHT_TILT.md) proves the reciprocal weighting used here.

Reproduce the input-weight experiments with the commands below. The optional
cache is omitted so each run reconstructs the operators independently.

```sh
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .04 .044 .048 --penalties .85 .9 --full-feedback 6 --window-histogram 8 --joint-cancellation --shell-points --probe-supports 176 184 192 200 208 216 224 --target-bits 52
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .044 .048 --penalties .9 .95 --weight-tilt 1.02 --full-feedback 6 --window-histogram 8 --joint-cancellation --shell-points --probe-supports 176 184 192 200 208 216 224 --target-bits 52
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .044 .048 --penalties .95 1 --weight-tilt .96 --full-feedback 6 --window-histogram 8 --joint-cancellation --shell-points --probe-supports 160 176 184 192 200 208 216 224 240 256 --target-bits 52
python -B research/workstreams/permutation_locality/independent_rows/verify.py --groups 64 --tilts .048 --penalties 1 --weight-tilt .94 --full-feedback 6 --window-histogram 8 --joint-cancellation --shell-points --probe-supports 160 176 184 192 200 208 216 224 240 256 --target-bits 52
```

## Next proof milestone

The independent-row distribution is now used inside the inner analysis,
rather than only in outer counts with worst-case packet shapes.
[ROW_CONDITIONING.md](ROW_CONDITIONING.md) gives the coefficient argument
and outward eleven-coordinate replay. At q = 64 it closes all messages
whose active groups each contain four rows with weights in [64,192],
at 79.92 bits. Two selected mixtures with three-active-row groups also
close. This is progress on different message classes, not a full cover
of the union-support events tabulated above.

The next milestone is to cover the remaining row types and occupancies.
The exact typed-placement implementation preserves each group's type
across regions. Hölder envelopes can combine related types, but the
coarse all-exceptions envelope loses too much. Larger occupancies also
need stronger bounds or further tilt tuning; mixed types are not the
only remaining issue.

Exact shell bounds are now integrated into interval folding through
`--shell-cover`; `--prefix-rank` jointly exploits shell and CDF constraints.
Those options reduce counting slack but cannot close a support point that
already fails with its individual shell cap. The old shared-shuffle proof
remains a separate active track, with its source and certificates unchanged.

## Regression checks

The independent-row exact support tests pass, including upward integer
CDF conversion. The driver also replays 2430 exact interval-fold checks,
adaptive-cover geometry and multiplicities, retained-parent checks, and
28 noncommuting placement-prefix identities. These support safe reuse of
the general cover machinery; they do not replace the probability argument
for the new outer measure in [README.md](README.md).
