# Checking whether Q1 represents the BCH failure bound

The updated analysis has thirty-two useful full bounds for BCH lengths 64
and 128, using their exact spectra and the recorded nested inner maps.
Twenty-nine have small aggregation losses. At the largest BCH-64 message size,
higher occupations cost 0.8201 bits: Q1 still contributes more than half
of the complete bound, but the rest are no longer negligible. BCH-256
remains outside this primary analysis.

All useful references use the same Q2..4 refinement procedure.
`refine_bch_full_anchors_v1.py` archived the eleven preceding references
and replayed their complete unions after refinement. The two new state-size
checks use the same refinement through the v3 batch. Eleven neighboring-state
references reuse those covers and independently replay every target bound.
The intermediate message-size checks and their state comparisons extend
these references. The ratios below compare upper-bound contributions,
not true event probabilities.

Selected anchors at K=2^20 and relative distance 1/10:

| BCH block | t | s | Q1 margin | Full margin | Higher occupations / Q1 | Margin loss |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 64 | 13 | 9.539167386 | 9.464552721 | 0.0530797268 | 0.0746146644 |
| 64 | 64 | 20 | 10.386947858 | 10.371924411 | 0.0104678685 | 0.0150234467 |
| 64 | 128 | 19 | 10.369172623 | 10.353299778 | 0.0110629650 | 0.0158728454 |
| 64 | 128 | 20 | 10.380476474 | 10.365380374 | 0.0105187561 | 0.0150960997 |
| 128 | 64 | 14 | 31.824890700 | 31.824890671 | 2.0155080e-8 | 2.9077636e-8 |
| 128 | 64 | 20 | 32.769598137 | 32.769598132 | 4.0024810e-9 | 5.7743582e-9 |
| 128 | 128 | 18 | 32.699950881 | 32.699950874 | 4.7960779e-9 | 6.9192782e-9 |
| 128 | 128 | 19 | 32.739116968 | 32.739116962 | 4.2924275e-9 | 6.1926642e-9 |
| 128 | 128 | 20 | 32.759208744 | 32.759208738 | 4.0336251e-9 | 5.8192933e-9 |

The message-size endpoints use t64/s20:

| BCH block | log2 K | Q1 margin | Full margin | Higher occupations / Q1 | Margin loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64 | 12 | 16.639437359 | 16.639117327 | 0.0002218537 | 0.0003200317 |
| 64 | 26 | 4.394562919 | 3.574424844 | 0.765574961 | 0.820138074 |
| 128 | 12 | 37.661913405 | 37.661913344 | 4.2448887e-8 | 6.1240796e-8 |
| 128 | 26 | 26.783189820 | 26.783189456 | 2.5228518e-7 | 3.6397053e-7 |

Margins and losses are in bits. At the small endpoint every occupation
fits in the explicit coefficient calculation. The largest BCH-64 union
includes Q5..512 with 7.834040415 margin bits and Q513..2,097,152 with
65.517619946 bits. Its refined Q2..4 interval has 4.965027813 bits.
Thus the observed loss is backed by complete occupation coverage.

The smaller t64 state settings now have useful full evidence: reducing s
from 20 to 13 costs about 0.9074 full-margin bits for BCH-64; reducing it
to 14 costs about 0.9447 bits for BCH-128. At s20, doubling t to 128 costs
only 0.00654 and 0.01039 bits, respectively. These are comparisons of
replayed bounds at the displayed geometries, not an implementation-cost
claim or a uniform interpolation over the missing points.

At K=2^20 and t64, every integer state from 13 through 20 for BCH-64
and 14 through 20 for BCH-128 now has a useful full bound. The eleven
new references retain the source integer partitions, re-evaluate all
sparse and dense witnesses at the target map, and pass the full v5
verifier. No monotonicity in s is assumed. The full curves and their
aggregation losses are in `bch_full_bound_state_scaling_v7.png`. BCH-128
s20 also replays the same source cover: its older loose dense tail had
created a 1.6e-9-bit upward artifact in the loss curve. The matched cover
reduces that penalty from 7.40e-9 to 5.77e-9 bits. This is a bound
refinement, not a smoothing operation or a change in the encoder.

Twenty complete upper bounds remain weak. At K=2^20 these are BCH-64
t128/s17..18, BCH-128 t128/s17, and the t64 states s9..12 for BCH-64
and s9..13 for BCH-128. At K=2^22 and 2^24, the transported t64/s10
and t64/s12 bounds remain weak for both constituents. Their zero-state
lower bounds are inconclusive. These are unresolved bounds, not code
counterexamples or first-moment obstructions.

Fixed-witness transport alone also gave weak bounds at t64/s16 for the
four intermediate-K cases. Fresh dense witness searches close all four,
while retaining their already useful transported sparse intervals.
This demonstrates why a poor transported upper bound should trigger a
target-map search before being interpreted as a parameter limitation.

The exact integer kernel scan separately finds 46 positive zero-state
first-moment lower bounds, all replayed at 90 digits. For both BCH blocks
at K=2^20, the detected settings are:

| t | Tested states with a positive first-moment lower bound |
| ---: | --- |
| 64 | 7 through 8 |
| 128 | 8 through 16 |
| 256 | 9 through 20 |

For example, the BCH-128 t256/s20 lower exponent is about 102,883.94 bits,
while Q1 alone has a positive margin near 32.75 bits. These trajectories
show why the Q1 surface cannot describe a small complete first-moment
bound throughout the t/s grid. No failure-probability lower bound is inferred.

`bch_engineering_evidence_v7.csv` labels all 130 geometries: twenty-nine full
bounds with small loss, three useful full bounds with larger losses, twenty
weak full upper bounds, 46 first-moment obstructions, and 32 sparse-only points.
`bch_full_bound_coverage_v7.png` shows the t/s coverage and continuous losses.
`bch_q1_vs_dense_tradeoff_v7.png` contrasts Q1 with the zero-state lower
exponents. White cells are outside the recorded grid; gray cells have
only sparse evidence.

`bch_full_bound_k_scaling_v7.png` compares the message-size anchors with
an explicit two-term engineering model calibrated at K=2^20. The model
predicts the K=2^26 full margins within 0.073 bits for BCH-64 and 0.014
bits for BCH-128. Its dashed curves are estimates; full-bound markers
appear only at replayed geometries. The normalization and its limitations
are explained below. The report uses each reference's selected sparse
components and retains the older coarse penalty separately.

All 121 tests pass. The full-reference replays, selected 90-digit checks,
and 46 positive lower-bound replays are separate from the test suite.
The new positive composition evaluator agrees with twelve retained BCH
witnesses within 1.1e-11 absolute log units. The v5 sparse producer saves
every completed occupation; a replay of 252 saved occupations leaves
both their files and the full interval byte-identical. These are audited
binary64 diagnostics, not outward arithmetic certificates.

The same-map transport checks reproduce both source sparse intervals
within 4.6e-13 absolute log units and both dense covers within 4.7e-10.
Every transported target then passes the independent full log-domain
replay and selected 90-digit checks. The target receipts authenticate
the unchanged source references and the transport implementation.

The four intermediate-K references and thirteen additional transports
pass full replay, including selected 90-digit checks. Four fresh dense
searches then strengthen the s16 references. The 90-digit count-scaling
check validates all six model comparisons with maximum absolute error
5.4e-15 bits, including tiny positive corrections from Q3 and higher.

Next, fill the log2 K=16 and 18 slices at s20, then s10/s12/s16. The
remaining s20-only gaps are exponents 13,14,15,17,19,21,23,25 for each
block. Weak-state refinements remain a separate task. The goal stays
active; no complete surface is asserted over unresolved points.

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
is recorded in `bch_k_counting_model_v3.json`.

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
nearly equal margins. `verify_bch_counting_model_v1.py` independently
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
`bch_full_bound_k_state_scaling_v7.png`.

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

The complete BCH-128 reference combines explicit region coefficients
through Q256 with a disjoint dense cover above that cutoff. BCH-64 uses
explicit region coefficients through Q1024 because its weak type boxes
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
