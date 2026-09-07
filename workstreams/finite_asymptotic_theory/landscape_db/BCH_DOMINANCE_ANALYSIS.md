# Checking whether Q1 represents the BCH failure bound

The BCH-64/128 calculations support a useful Q1 approximation in part of
the engineering surface, and expose substantial higher-occupation effects
elsewhere. At T=64,S=20, BCH-128 has negligible aggregation loss throughout
the tested message-size curve. BCH-64 develops a visible loss as K grows.
Reducing S makes that loss appear earlier. Increasing T can introduce a
much larger first-moment obstruction that the Q1 curves do not reveal.
BCH-256 is excluded; both constituents here use exact spectra.

Here K is the message size in bits, B is the BCH block length, T is the
number of input bits per inner epoch, and S is the state size in bits.
The relative-distance target is 1/10. The nested inner maps are fixed to
the recorded constructions. Expectations are over encoder permutations
and fresh nonzero multipliers. Q counts nonzero outer rows, with
L=K/(B/2) rows in total. A full margin is minus the base-two logarithm of
the selected upper bound on the expected number of bad words. A positive
margin gives a candidate Markov failure bound; the numerical evaluations
below are audited binary64 diagnostics, not outward arithmetic certificates.

## Coverage and what it establishes

The original grid has 130 geometries, comprising the following slices for
each BCH block. Their intersections are counted only once.

| Slice | Parameters |
| --- | --- |
| Message size | T=64, S=20, every integer log2 K from 12 through 26 |
| Message/state interaction | T=64, S in {10,12,16}, log2 K in {16,18,20,22,24} |
| Epoch/state interaction | K=2^20; T=64 with S=7..20, T=128 with S=8..20, T=256 with S=9..20 |

The complete-grid audit passes for all 130 geometries: 84 have replayed
complete upper bounds, of which 77 are useful and seven remain weak; the
other 46 have positive first-moment lower obstructions. No sparse-only
points remain. Among the 77 useful bounds, 68 have higher/Q1 upper-bound
ratios at most 0.1, and nine have larger corrections. Q1 is a majority of
the selected bound at 76 of these geometries, but not at BCH-64 T64/S10,
K=2^22. Complete evidence coverage does not assert a useful upper bound
throughout a continuous surface.

The seven unresolved upper bounds are:

| BCH block | T | S | log2 K | Full margin (bits) | Weak part of selected bound |
| --- | ---: | ---: | ---: | ---: | --- |
| 64 | 64 | 9 | 20 | -17824.2664 | Dense Q257..L |
| 128 | 64 | 9 | 20 | -5127.0658 | Dense Q257..L |
| 64 | 128 | 17 | 20 | -25469.7139 | Dense Q257..L |
| 64 | 128 | 18 | 20 | -10331.1951 | Dense Q257..L |
| 128 | 128 | 17 | 20 | -12756.4797 | Dense Q257..L |
| 64 | 64 | 10 | 24 | -139.6375 | Q1..4 union, Q5..256, and dense tail |
| 64 | 64 | 12 | 24 | -1.6440 | Q5..256 union |

These are inconclusive bounds, not proved code failures. At the final
S=12 point, Q1..4 has a useful 3.6064-bit margin and the dense tail has
55.2544 bits; the sparse interval contributes -1.6056 bits and determines
the full deficit. Its weakest individual witness is Q62, at 0.13484 bits.
At S=10,K=2^24, Q1..4 already has a -0.71676-bit union margin, so fixing
only the dense tail cannot make the selected full bound useful. Fresh
searches improved both points, but did not resolve them.

## How the engineering curves change

At T=64,S=20, selected full margins and the loss from adding Q>=2 are:

| BCH block | log2 K | Full margin | Full-bound loss relative to Q1 |
| --- | ---: | ---: | ---: |
| 64 | 12 | 16.639117327 | 0.000320032 |
| 64 | 20 | 10.371924411 | 0.015023447 |
| 64 | 24 | 6.171264034 | 0.222935387 |
| 64 | 26 | 3.574424844 | 0.820138074 |
| 128 | 12 | 37.661913344 | 6.12408e-8 |
| 128 | 20 | 32.769598132 | 5.77436e-9 |
| 128 | 24 | 28.782539084 | 9.10589e-8 |
| 128 | 26 | 26.783189456 | 3.63971e-7 |

All values are in bits. Every integer message exponent from 12 through 26
now has a useful full bound at this T,S setting for both BCH blocks.
The small-K curves include a transfer transition; they should not be fit
to the same large-K slope. Once the per-row terms stabilize, doubling K
roughly costs one Q1 margin bit and doubles the relative pair contribution.
The latter makes the full BCH-64 curve bend away from Q1. For BCH-128 the
normalized pair coefficient is about twenty bits smaller, so its analogous
correction remains tiny over the tested range.

The unchanged two-term engineering model is calibrated only at K=2^20.
It is checked at all twelve larger-K comparison points, with no refitting.
At K=2^26 its margin estimate is optimistic by 0.07283 bits for BCH-64
and conservative by 0.01360 bits for BCH-128. These are finite checks of
a model of selected bounds, not a theorem about extrapolation or actual
failure probabilities. The normalization and error decomposition follow below.

The state comparison has a second effect. At K=2^20,T=64, reducing S
from 20 to 10 changes the full BCH-64 margin from 10.3719 to 6.6743 bits,
and increases the aggregation loss from 0.0150 to 0.5844 bits. At the same
parameters BCH-128 changes from 32.7696 to 27.0410 margin bits, but its
aggregation loss is only 2.14e-6 bits at S=10. Thus both blocks lose margin
when S shrinks, while BCH-64 also incurs a substantial higher-occupation
correction. At BCH-64 S=10,K=2^22 the full margin is 3.5754 bits and the
loss is 1.6949 bits: the selected higher-occupation upper bounds together
exceed the selected Q1 bound. This statement compares upper bounds; it
does not establish the ratio of the true contributions.

The middle state S=16 retains a modest margin cost. At K=2^24, reducing
S=20 to S=16 costs 0.27004 full-margin bits for BCH-64 and 0.27442 bits
for BCH-128. The corresponding BCH-64 aggregation losses are 0.22294
and 0.34559 bits. This supports an engineering tradeoff between state
size and margin, without using a fixed 40-bit pass/fail threshold.

## Why T needs its own evidence

The exact integer kernel scan finds 46 positive first-moment lower bounds,
all replayed at 90 digits. For both BCH blocks at K=2^20:

| T | Tested states with a positive first-moment lower bound |
| ---: | --- |
| 64 | 7 through 8 |
| 128 | 8 through 16 |
| 256 | 9 through 20 |

Multi-bit inputs can lie in the kernel of the syndrome map and keep the
state zero. Their output then equals their input. The number of such
low-output trajectories can overwhelm the Q1 contribution. For BCH-128
T=256,S=20, the lower exponent on the expected bad-word count is about
102,883.94 bits, despite a Q1 margin near 32.75 bits. A small first-moment
bound is impossible at that geometry under the analyzed construction.
This does not imply a large code-failure probability: many bad words may
coincide in rare setups.

The CSV quantifies this discrepancy. If the total first moment is at least
2^a and its Q1 contribution is at most 2^(-m), then the true higher/Q1
first-moment ratio is at least 2^(a+m)-1, whenever the denominator is
positive. The zero-denominator case has only higher-occupation mass.
The report stores the logarithm of this lower ratio and the resulting
upper limit -a on any possible first-moment margin. These lower witnesses
are different evidence from ratios between selected upper bounds.

## Plots, validation, and search quality

`bch_full_bound_k_scaling_v8.png` shows every S=20 message-size anchor,
the Q1 curve, and the fixed-calibration model. Its loss panels label the
absolute size of the corrections so a log axis does not make tiny BCH-128
corrections look practically large. `bch_full_bound_k_state_scaling_v8.png`
compares S=10,12,16,20 across K; weak upper bounds interrupt those curves.
`bch_full_bound_state_scaling_v8.png` shows the integer state sweep.
`bch_full_bound_coverage_v8.png` identifies useful bounds, weak bounds,
and first-moment obstructions at K=2^20. White cells are outside the grid.
`bch_q1_vs_dense_tradeoff_v8.png` explains why the nearly flat Q1-versus-T
curves do not describe the full first moment.

Each full reference replays Q2..4 compositions, every selected sparse
witness, and a complete dense integer partition when needed. The intervals
cover every Q from 1 to L without gaps. The verifier uses an independent
log-domain composition evaluator and selected 90-digit checks. The original
121-test core pass remains applicable because its numerical sources were
not changed. Model arithmetic and the full grid receive additional checks
through `verify_bch_counting_model_v2.py` and `audit_bch_evidence_grid_v1.py`.
Both pass. The model check has maximum absolute discrepancy 5.35e-15 bits;
the grid audit authenticates 602 files and reconciles aggregate margins and
ratios within 2.54e-15 bits. All five regenerated figures were visually
inspected. These arithmetic checks do not measure the slack of the bounds.

Reusing witnesses can save search effort, but the target bounds must be
recomputed. Fresh target-specific dense searches made 22 of 24 selected
weak full bounds useful. Fresh short-K searches also repaired weak
transports near Q=L. At low S and larger K, sparse witnesses need their
own optimization. These improvements measure search slack; a poor
transported bound is not a code counterexample or a reliable scaling trend.

Numerical producers and verifiers are versioned and frozen once receipts
refer to their hashes. Parent references used by transports remain immutable;
replaced targets are archived locally. Generated receipts, plots, CSVs,
and caches stay ignored. Git contains source and explainers only.

## What would justify using the Q1 surface?

Fix one geometry (K,t,s,B), relative distance 1/10, and the recorded inner
maps. Let Q count the nonzero outer rows and let L=K/(B/2). For each Q,
let U_Q be an upper bound on that occupation's first-moment contribution.
The expectation is over the encoder permutations and fresh multipliers,
for the fixed BCH constituent and inner maps.

Write U_rest=sum_(Q=2)^L U_Q. The complete first-moment bound is
U_full=U_1+U_rest. The loss relative to the Q1 margin is

    Delta = log2(1 + U_rest/U_1).

If U_rest/U_1 <= 0.1, the loss is at most 0.138 bits. This establishes
that adding the remaining upper bounds barely changes the Q1 certificate
candidate. It does not establish the relative sizes of the true failure
events: upper bounds alone need not preserve their ratios.

Every reported ratio must identify its occupation coverage. A Q2..4 ratio
does not control Q5 and higher. A loose full upper bound neither establishes
Q1 dominance nor refutes the distance property. Evidence at finitely many
geometries does not establish dominance between them without an additional
uniform argument.

## Separating message-size counting from the transfer

For one fixed constituent and inner map, write each chosen occupation
bound as U_Q=choose(L,Q) H_Q(L), where L=K/(B/2). This defines H_Q(L)
by removing the choice of active outer rows. The dependence of H_Q on L
still includes the transfer, the distance cutoff, and the optimized
counting witnesses. In particular,

    U_2/U_1 = (L-1)/2 * H_2(L)/H_1(L),
    U_3/U_1 = (L-1)(L-2)/6 * H_3(L)/H_1(L).

These identities explain what a message-size extrapolation must check.
If H_2/H_1 has stabilized, doubling K approximately doubles the Q2/Q1
ratio, even when Q1's own margin follows a nearly straight line in log2 K.
The corresponding full-margin loss grows as log2(1+U_rest/U_1), rather
than as a fixed offset. The small-K points can have substantial transfer
effects, so the counting factors alone do not establish their slope.

This is a conditional engineering model. Similar values at selected K
do not prove a uniform bound between them. A fixed-Q model also cannot
exclude a contribution with Q proportional to L; that requires the
complete-tail evidence retained separately in this audit. Comparisons
should use similarly refined witnesses, since a changing search budget
can otherwise look like a change in the code's scaling.

The matched t64/s20 references give the following measured coefficients:

| BCH block | log2(H2/H1), K=2^20 | log2(H2/H1), K=2^26 |
| --- | ---: | ---: |
| 64 | -20.578141 | -20.589313 |
| 128 | -40.896370 | -40.918440 |

Their changes are only 0.0112 and 0.0221 bits. Calibrating H1 and H2 at
K=2^20 and retaining the exact row-count factors predicts Q1+Q2 margins
of 3.647249954 and 26.769597768 bits at K=2^26. The complete bounds there
have 3.574424844 and 26.783189456 bits. The first estimate is optimistic
by 0.07283 bits; the second is conservative by 0.01360 bits. The comparison
is recorded in `bch_k_counting_model_v4.json`.

This supports a local engineering explanation: the per-row transfer terms
are nearly stable over these large-K anchors, while the count of pairs
grows quadratically. BCH-64 therefore develops a visible aggregation loss
at the large endpoint. The normalized pair coefficient for BCH-128 is
about twenty bits smaller, so its loss is still tiny in the studied range.
The short-K coefficients differ substantially and should not be fit to
the same plateau. No untested large-Q contribution is excluded by this
model; the complete endpoint replays supply that evidence separately.

## Intermediate message sizes and the state interaction

The model retains its K=2^20 coefficients. Complete bounds at the two
intermediate message sizes give the following comparisons at t64/s20;
negative error means the estimate is conservative relative to the
selected full upper bound.

| BCH block | log2 K | Two-term estimate | Full margin | Estimate minus full |
| --- | ---: | ---: | ---: | ---: |
| 64 | 22 | 8.327781119 | 8.333885368 | -0.006104249 |
| 64 | 24 | 6.163579314 | 6.171264034 | -0.007684720 |
| 128 | 22 | 30.769598114 | 30.780002354 | -0.010404240 |
| 128 | 24 | 28.769598045 | 28.782539084 | -0.012941039 |

All four discrepancies are below 0.014 bits without refitting. The report
separates the error in predicting the Q1+Q2 sum from the penalty for adding
Q3 and higher. At BCH-64 K=2^26, those terms are -0.01209 and 0.08491 bits,
respectively, for a net optimistic error of 0.07283 bits. The latter
penalty includes the slack in the selected higher-occupation bounds; it
is not a measurement of the corresponding true event probabilities.
Tiny BCH-128 corrections are computed with log1p rather than subtracting
nearly equal margins. `verify_bch_counting_model_v2.py` independently
checks the predictions by scaling U1 by L/L0 and U2 by
choose(L,2)/choose(L0,2), using 90-digit arithmetic.

The full K/state comparison uses t64 and relative distance 1/10:

| BCH block | log2 K | Full margin at s16 | Full margin at s20 | Cost of using s16 |
| --- | ---: | ---: | ---: | ---: |
| 64 | 20 | 10.215545930 | 10.371924411 | 0.156378481 |
| 64 | 22 | 8.151754083 | 8.333885368 | 0.182131286 |
| 64 | 24 | 5.901222091 | 6.171264034 | 0.270041943 |
| 128 | 20 | 32.494912065 | 32.769598132 | 0.274686066 |
| 128 | 22 | 30.505581685 | 30.780002354 | 0.274420669 |
| 128 | 24 | 28.508121661 | 28.782539084 | 0.274417424 |

The BCH-128 state cost is nearly constant across these three K values.
For BCH-64 it increases with K, as the higher-occupation correction grows.
At K=2^24, the correction is 0.3456 bits at s16 versus 0.2229 bits at s20.
This gives a quantitative state-size tradeoff supported by complete
bounds at every displayed point. It does not extend the s16 conclusion
to K=2^26 or to smaller unverified states. The curves are in
`bch_full_bound_k_state_scaling_v8.png`.

## Multi-bit transfer with uniform refresh

An epoch takes a field state x and binary t-vector u. It outputs u+A(x),
then updates x to alpha*x+J(u), where alpha is fresh and nonzero. The state
starts at zero. All studied geometries use complete epochs and no flush.
Let z=exp(-lambda), with lambda>0. The current input is uniform among the
weight-j t-vectors, independently of the incoming weighted state within
the corresponding coefficient calculation.

Put M=2^s-1 and kappa=M/(M-1). The transfer retains four classes for the
normalized state distribution after weighting earlier output moments:

- Z: zero state;
- D: arbitrary nonzero state;
- U: uniform nonzero state;
- L: nonzero state with point probabilities at most 1/(M-1).

Each class also carries its unnormalized weighted mass. Let a_w count the
nonzero inner words of weight w, and let d_A be their minimum weight. Define

    h_j(w) = sum_v choose(w,v) choose(t-w,j-v)
                   z^(w+j-2v) / choose(t,j),
    r_D,j = max_(w:a_w>0) h_j(w),
    r_U,j = sum_w a_w h_j(w)/M,
    r_L,j = min(r_D,j, kappa*r_U,j).

Let k_j count the weight-j vectors in ker(J), and put
beta_j=k_j/choose(t,j). For j>0 the zero-state transitions are

    Z -> Z : beta_j z^j,
    Z -> D : (1-beta_j) z^j.

These transitions retain nonzero inputs with zero syndrome explicitly.
For a live class with moment upper bound r, define

    k = min(r, beta_j z^max(0,d_A-j)),
    n = min(r, (1-beta_j) z^max(0,d_A-j)).

The following live-state transitions are upper bounds:

    live -> Z : n/M,
    live -> L : (1-1/M) r + k/M.

To see this, let R be the actual output moment and K_0 its zero-syndrome
part. Fresh multiplication sends every nonzero prestate uniformly over
the nonzero field elements, independently of the already emitted output.
A nonzero syndrome cancels with probability 1/M. A zero syndrome never
cancels. Thus the actual zero and surviving masses are

    (R-K_0)/M,    (1-1/M)R + K_0/M.

The displayed bounds follow from R<=r, K_0<=k, and R-K_0<=n. A surviving
nonzero-syndrome update has density at most 1/(M-1); a zero-syndrome update
is uniform. Their mixture therefore belongs to L. This argument requires
no independence between the output weight and the syndrome.

For j=0, the transfer uses the stronger exact uniform refresh:

    Z -> Z : 1,
    D -> U : z^d_A,
    U -> U : sum_w a_w z^w/M,
    L -> U : min(z^d_A, kappa sum_w a_w z^w/M).

Positive matrix-polynomial composition averages the j active positions
within a region and preserves the state between regions. The existing
fixed outer-weight composition and typed-count inequalities then apply
to these four-state matrices. Their counting measures remain fixed for
each witness; they are not resampled constituents.

## Search and evidence discipline

### Activation density for the larger-t follow-up

The t128 follow-up uses additional information about the state created by
J(u). Here J is the transpose of the recorded generator of A. Let K_j(w)
be the coefficient of x^j in (1-x)^w(1+x)^(t-w). Character inversion gives
the following upper bound on the probability of any particular syndrome
when u is uniform among the weight-j inputs:

    pmax_j = [choose(t,j) + sum_(w>0) a_w |K_j(w)|]
             / [2^s choose(t,j)].

The signed sum without absolute values recovers the exact zero-syndrome
probability beta_j. The implementation checks this identity against
every recorded kernel coefficient using integer arithmetic.

For nonzero syndromes, cap pmax_j by 1-beta_j. With M=2^s-1, define

    W_j = max(1-beta_j, (M-1) min(pmax_j,1-beta_j)).

The nonzero syndrome measure has total mass at most W_j and every point
has mass at most W_j/(M-1). It can therefore be dominated by W_j times
a distribution in class L: add the missing mass while respecting that
pointwise cap. This gives the alternative zero-state transitions

    Z -> Z : beta_j z^j,
    Z -> L : W_j z^j.

All live-state transitions remain as above. The current policy uses this
alternative when W_j <= 4(1-beta_j), and retains class D otherwise.
The factor four is a choice in the bound, not an encoder parameter.
`syndrome_density_v1.py` implements the alternative. Exhaustive small
examples check every syndrome probability, future-state potential, and
three-epoch path. This extension changes no previously recorded bound.

The t128 search also separates the all-one outer word into its own band
and optimizes the other band probabilities jointly with the tilt and
type proposal. Smoothed costs guide that search; every retained bound
uses the original unsmoothed costs and direct transfer evaluation.

Three individual types selected from weak BCH-128 t128/s20 boxes improve
from approximately -65,000 margin bits to 20,822.75, 47,773.90, and
234,369.67 bits under direct witness refinement. These are selected
binary64 point diagnostics. They demonstrate search slack at those
types, not closure of their surrounding boxes or the full interval.
The refinement isolates zero-count faces and uses relative count
uncertainty to choose subdivisions. Every split is checked as an exact
disjoint integer partition. The five-category t128/s20 search remained
weak after 2,500 refinements; the smaller cover below closes that geometry.

The bulk dense-tail alternative replaces the nontrivial weight bands by
one pointwise binomial majorant. For a fixed p in (0,1), set

    Gamma(p) = max_(0<w<B, A_w>0)
                 A_w / [choose(B,w) p^w (1-p)^(B-w)].

The counting measure of every nonzero, non-all-one BCH row is dominated
by Gamma(p) times the Bernoulli-p row law. The zero row is kept separate,
as is the all-one row when the recorded spectrum contains it. Thus the
dense type vector has two coordinates for BCH-64 and three for BCH-128.
The recorded BCH-64 spectrum has maximum weight 56 and no all-one word;
BCH-128 has one all-one word. This is an upper bound from the exact fixed
spectra, not a random-code ensemble assumption. Its weaker weight
information is offset by the smaller type cover at t128/s20 for both
blocks, as the complete bounds above demonstrate. The producers are
`seed_bch_dense_v3.py` and `close_bch_dense_v11.py`.

The same bulk bound for BCH-128 at t128/s17 remains inconclusive:
1,000 refinements leave a complete dense upper-bound margin near
-12,756.48 bits. Subdivision has stalled at that value. Its exact-kernel
zero-state lower bound is also inconclusive, near -9,512.48 log2 units.
Neither calculation establishes an obstruction or a useful full bound
at s17. The lower convex witness already uses adjacent regional weights
1638 and 1640, so replacing its broad-range convex hull by a variance
constraint alone has no evident benefit at that witness.

For BCH-64 at t128/s18, the bulk cover similarly stalls near
-10,331.20 bits. Both this geometry and BCH-128 t128/s17 have now had
their complete unions replayed; they are weak full upper bounds, rather
than missing tail calculations. No failure-probability conclusion follows.
In the BCH-64 cover, the limiting types have about 11,432 active rows,
and the bulk probability is near 0.5055. A future refinement can separate
the exceptional weight-56 shell before increasing the type subdivision
budget: the bulk majorant at p=1/2 pays 33.542 bits per active row,
whereas the largest remaining shell cost is 32.211 bits. Whether that
change closes s18 requires a complete new cover and replay.

`occupation_refresh_v1.py` is a separate evaluator. Historical producer
sources and their receipts remain unchanged. Exhaustive GF(4) checks cover
every epoch occupation, extreme distributions in all four classes, and
three successive epochs. Positive-arithmetic checks also compare region
coefficients and serialized moments against direct products.

The sparse contraction also has a positive-arithmetic implementation in
`occupation_composition_positive_v1.py`. Its envelope is a maximum of
positive linear forms, so scaling one matrix entry commutes with every
envelope step. Direct positive correlations replace the temporary array
of logarithmic summands. Very small correlation outputs are recomputed
in log space; coefficient-scale guards, support tracking, and periodic
normalization similarly send unsafe envelope entries to the original
log evaluator. The retained bounds are still replayed by the original
log implementation in the full verifier.

Three tests cover ordinary and rare coefficients, deterministic support,
folds through 512 steps, and complete composition matrices. Twelve saved
BCH boxes agree within 1.1e-11 log units. A selected-case profile supports
using the new kernel in subsequent searches; it is not an encoder benchmark.
`close_bch_sparse_tail_v5.py` uses this kernel and atomically saves each
completed occupation. `verify_bch_sparse_resume_v1.py` confirms a completed
252-occupation run resumes without changing any saved evidence.

The sparse search retains separate outer-band compositions. It explores
a wider range of auxiliary Bernoulli probabilities than the historical
search. These probabilities belong to the counting inequality, not to
the encoder. Optimizing them does not change the code or its randomness.

The dense search partitions all integer row-type counts, including zero
rows. Each selected box uses one fixed positive coefficient proposal and
one tilt. A lookup table helps choose proposals, but every reported box
bound uses a direct transfer evaluation. Interpolated values never enter
the selected cover's sum.

`study_bch_dominance_v1.py` evaluates Q2..4 at the 130 existing BCH-64/128
engineering geometries. It records composition-level witnesses, source
hashes, and the added margin penalty. Its full-dominance status remains
unresolved until Q5 and higher are controlled. Generated CSV, JSON, logs,
and figures stay local and ignored by Git.

The next step is to aggregate authenticated complete tail covers at
representative extremes and around the observed Q1 state-size knees.
Until those checks succeed, the earlier curves describe Q1 contributions.

## Distinguishing search slack from cancellation effects

The 130-point sparse pass gives Q2..4 at most one tenth of the Q1 bound
at 98 points. Twelve selected composition witnesses replay at 90 decimal
digits, with log discrepancies below 1e-10. These remain sparse bounds;
they do not establish full dominance at 98 points.

Coarse witness searches can produce misleading spikes. At BCH-64,
K=2^26, t64/s20, continuous refinement reduces the Q1..4 aggregation
penalty from 5.54 to 0.743 bits. For one nearly fixed BCH-128 dense type
at Q150..151, joint tilt and proposal optimization improves its margin
from -1,249.49 to 1,784.88 bits. Neither initial deficit was evidence of
an intrinsic cancellation problem. Every improved bound is evaluated
directly; interpolation only guides the search.

The initial complete BCH-128 reference combines explicit region coefficients
through Q256 with a disjoint dense cover above that cutoff. The initial
BCH-64 reference uses explicit region coefficients through Q1024 because its weak type boxes
were concentrated at occupations in the hundreds. Clipping and refining
the saved cover above that cutoff closes the remaining interval. Both
complete aggregations have passed replay.

A positive lower bound on the first moment is different evidence. It
cannot be removed by optimizing an upper-bound witness. The zero-state
calculation below supplies that distinction for the larger t settings.

## An independent test for a first-moment obstruction

A separate lower bound can distinguish upper-bound slack from a limitation
of the first-moment method. It restricts attention to trajectories whose
state stays zero at every epoch. Their output equals their input, so any
such message of weight at most floor(BL/10) is bad for the distance target.
A large lower bound on the expected number of these messages prevents a
small first-moment certificate. It does not imply a large probability of
code failure; many bad messages may occur in the same rare setup.

Fix an exact outer weight w with multiplicity N_w, and an even occupation
Q satisfying Qw <= floor(BL/10). There are choose(L,Q) N_w^Q such messages.
After the independent outer-row permutations, each row support is a
uniform weight-w subset of the B regions. Let q_i count active inputs in
region i. Each marginal q_i is Binomial(Q,w/B), and sum_i q_i=Qw.

The probability that every q_i is even is at least 2^(1-B). To establish
this, take the Fourier transform of the row-support distribution on the
binary B-cube. Its character averages are real numbers. The return
probability after Q rows is 2^-B times the sum of their Qth powers.
Every term is nonnegative for even Q. The zero character contributes one,
and the all-one character contributes ((-1)^w)^Q=1.

Choose a regional count interval [a,b] whose total outside probability,
summed over all B regions, is at most 2^(-B-4). Binomial Chernoff bounds
supply this interval without assuming independence between regions.
The probability of even counts throughout that interval is then at least

    g = (1-1/32) 2^(1-B).

Let k_j be the exact kernel weight counts of J. Conditional on q active
positions in one region, the probability of zero syndrome at every epoch is

    R_q = [x^q] (sum_j k_j x^j)^(L/t) / choose(L,q).

Different regional shuffles are independent. Let h be the lower convex
hull of the points (q, log R_q) for the even q in [a,b]. For every allowed
profile, convexity and the fixed total give

    sum_i log R_(q_i) >= B h(Qw/B).

Consequently the bad-word first moment for this restricted class is at least

    choose(L,Q) N_w^Q g exp(B h(Qw/B)).

`bch_zero_state_lower_v1.py` evaluates this expression using the exact
outer spectrum and inner kernel counts. Its positive coefficient powers,
parity bound, concentration charge, and convex-hull step have separate
combinatorial tests. Numerical lower bounds remain binary64 diagnostics
until independently replayed or rounded outward.

The exact-coefficient refinement uses integer arithmetic before taking
logarithms. With F(y)=sum_j k_(2j)y^j, E=L/t, and b_n=[y^n]F(y)^E,
the identity F(F^E)'=E F' F^E gives

    n b_n = sum_(j=1)^min(t/2,n) ((E+1)j-n) k_(2j) b_(n-j).

Here b_0=1. The implementation checks exact divisibility by n and
nonnegative coefficients. This short recurrence preserves cancellation
exactly and avoids quadratic polynomial convolution. Separate tests
compare its coefficients with direct integer polynomial products.
The logarithms, supporting lines, and concentration charge still need
high-precision replay; exact coefficients alone do not make the complete
lower bound an outward certificate.
