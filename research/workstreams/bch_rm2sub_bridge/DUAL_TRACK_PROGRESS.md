# BCH-256: certified improvement and growth-curve evidence

This work separates two questions. Can we improve the existing distance
certificate? Can BCH-256 inform extrapolation from smaller, known spectra?
Both calculations use the fixed BCH [256,128] outer and the selected t64/s20
map in `T64_S20_FULL_CLOSURE.md`. This map differs from the nested map used
in the finite-theory agent's parameter study.

## Track 1: an outward-certified Q1 improvement

At K=2^20, the full certified margin improves from 50.4390685433 to
50.4872982775 bits. Only the occupation-one contribution is replaced.
Every higher-occupation bound, the selected map, and cutoff 209716 are unchanged.
The new exact rational bound is `combined_full_upper` in
`generated/refresh_q1_outward_v1.json`. Its 512-bit replay checks all 92
shell coefficients and the rational aggregation. No BCH spectrum assumption
enters this result. The randomness model and limitations remain those of
`T64_S20_FULL_CLOSURE.md`.

The new checker implements the uniform-refresh argument from
`../finite_asymptotic_theory/landscape_db/BCH_GROWTH_ANALYSIS.md` in positive
Arb arithmetic. It also uses finer tilt witnesses. Thus the total 0.048-bit
gain is not attributed solely to the new state classes.

For completeness, fix one epoch input u and nonzero prestate x. The encoder
outputs u+A(x), then updates to alpha*x+B(u). The fresh multiplier alpha is
uniform over the nonzero elements of GF(2^20). The current output weight does
not depend on alpha. Consequently, weighting that output does not bias alpha.

Let M=2^20-1 and kappa=M/(M-1). We retain four classes of normalized state
laws, each with an unnormalized weighted mass: zero (Z), arbitrary nonzero
(D), uniform nonzero (U), and nonzero density at most 1/(M-1) (L).
For z in (0,1), let a_v count nonzero inner words of weight v and let d_A
be their minimum weight. Define

    m0 = sum_v a_v z^v / M,
    m1 = sum_v a_v [v z^(v-1) + (t-v) z^(v+1)] / (t M).

With zero input, the upper weighted transitions are

    Z -> Z: 1; D -> U: z^d_A; U -> U: m0; L -> U: kappa*m0.

With one uniformly positioned input bit, Z goes to D with weight z.
From D, U, or L, respectively, put r=z^(d_A-1), m1, or kappa*m1.
The outgoing weights are r/M to Z and r(M-1)/M to L.
The one-bit syndrome is nonzero for the fixed map. Exactly one multiplier
cancels it. Conditional on survival, the updated state is uniform outside
zero and that syndrome. Mixtures preserve the stated density bound.
For zero input, multiplication restores exact uniformity. These facts justify
positive composition across epochs, including after weighting past outputs.

The checker deliberately omits the diagnostic evaluator's optional minimum
between the density and pointwise bounds. Both versions are upper bounds;
omitting that minimum simplifies outward arithmetic. Region matrices average
the active epoch over its positions. A positive polynomial recurrence averages
outer supports by binomial(256,w), and the existing BCH dual bounds the shell sum.

The producer uses 256-bit Arb and binary epoch powering. The replay uses
512-bit Arb and linear epoch iteration. They share the epoch formulas and
outer-support recurrence, so the replay is not an independent proof of those
formulas. Exact rational tests additionally compare the new envelope against
an actual small GF(4) encoder and enumerate support placements.

## Track 2: separate transfer loss from unknown spectrum

The following Q1 diagnostics use cutoff floor(2K/10), complete epochs, and
the same selected map at every size. The first two columns use the same
certified BCH constraints, evaluated in nearest binary64. They are not new
outward certificates at those sizes. The third column substitutes the
conditioned-even binomial spectrum, not the true BCH spectrum.

| log2 K | Old transfer, BCH constraints | Refresh transfer, BCH constraints | Refresh, modeled spectrum |
| ---: | ---: | ---: | ---: |
| 16 | 54.156 | 54.156 | 74.217 |
| 20 | 50.479 | 50.488 | 70.554 |
| 24 | 46.364 | 46.509 | 66.576 |
| 28 | 40.202 | 42.510 | 62.577 |
| 32 | 1.577 | 38.510 | 58.577 |

The comparison uses identical tilt witnesses for the old and new transfers.
At large K, repeated density penalties in the old transfer create a substantial
artificial decline. The refreshed spectrum-constrained curve instead loses
approximately one margin bit per doubling over the larger tested sizes.
The modeled-spectrum curve has the same local behavior. Higher occupancies
have not been recomputed for these larger message lengths.

The spectrum model reserves one word at weights 0 and 256 and assigns
(2^128-2) times normalized binomial mass to even weights 38 through 218.
It predicts weight 38 as dominant at every tested K. Removing the row-count
factor gives C=M1+log2(K/128), which approaches 83.577 bits in these samples.
The persistent-state onset model gives a Chernoff intercept of 83.611 bits.
Their agreement supports the transfer explanation conditional on the same
spectrum model; it is not independent evidence that BCH-256 has that spectrum.

The roughly 20.07-bit separation between the constraint and modeled curves
therefore concerns uncertainty in the outer spectrum in this regime, not
an observed deterioration in the actual BCH code. Neither endpoint measures
the true failure probability. In particular, fitting the constraint curve as
though it came from an exact spectrum would confound code behavior with proof loss.

As a sensitivity check, multiplying modeled counts at weights 38 through 50
and their complements by 2^b reduces the modeled margin by almost b bits,
for b in {1,4,8,12,16}. These inflated arrays are stress envelopes, not
necessarily realizable code spectra; their total mass is not renormalized.
This translates the estimator's desired accuracy into a concrete requirement
on the dominant low-shell counts. For example, an eight-bit tail inflation
reduces the K=2^32 modeled Q1 margin from 58.577 to 50.577 bits.

The onset model's exact probability has intercept 87.423 bits, about 3.81
bits above its Chernoff bound. This is a possible direction for tightening
the SPIN analysis, not 3.81 additional proved bits for the real encoder.

Existing search floors and algebraic evidence remain relevant to testing the
spectrum assumption. The discovery-only partial sum in the diagnostic output
is not an upper or lower bound on actual setup failure. Biased orbit searches
do not certify completeness. This pass quantifies the sensitivity; it does
not upgrade those searches to confidence intervals or fit a new size law.

## Reproduction and next work

Restore the local BCH proof inputs using `MIGRATION.md`, then run sequentially:

```text
python -B workstreams/bch_rm2sub_bridge/test_refresh_q1.py
python -B workstreams/bch_rm2sub_bridge/dual_track_q1.py --output workstreams/bch_rm2sub_bridge/generated/dual_track_q1_new.json
python -B workstreams/bch_rm2sub_bridge/certify_refresh_q1.py --discovery workstreams/bch_rm2sub_bridge/generated/dual_track_q1_v1.json --output workstreams/bch_rm2sub_bridge/generated/refresh_q1_outward_new.json
python -B workstreams/bch_rm2sub_bridge/certify_refresh_q1.py --discovery workstreams/bch_rm2sub_bridge/generated/dual_track_q1_v1.json --output workstreams/bch_rm2sub_bridge/generated/refresh_q1_outward_new.json --verify
```

Existing result files are write-once. The diagnostic uses cutoff 209715 at
K=2^20; the certificate retains the original, stronger cutoff 209716.
The three small v1 diagnostic/producer/replay files are retained for audit.
Their source dependencies include local inputs; a fresh clone still requires
the migration restore step for a complete replay.

Next on track 1: certify selected larger-K Q1 points and optimize the BCH
dual for the refreshed shell coefficients. Coordinate higher-occupation
refresh work with the finite-theory agent; multiple-bit zero syndromes require
different transitions from Q1.

Next on track 2: compare the predicted low-shell counts against the BCH-256
search and incidence evidence, and assess which multiplicative error ranges
are defensible. Use held-out known spectra to test that assessment procedure.
Do not infer a full-distance usable-K limit from the Q1 curves alone.
