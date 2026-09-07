# Checking whether Q1 represents the BCH failure bound

The Q1 engineering curves do not yet establish the full failure margin.
The present task checks BCH lengths 64 and 128, using the exact spectra
and nested inner maps from the engineering study. BCH-256 remains outside
this primary analysis.

Seven useful complete reference bounds are verified at K=2^20: both BCH
sizes at t64/s20 and t128/s19..20, plus BCH-128 at t128/s18.
Q1 dominates the bound at these seven geometries. The verified
ratios compare upper-bound contributions, not actual event probabilities.

| BCH block | t | s | Q1 margin | Full margin | Higher occupations / Q1 | Margin loss |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 64 | 64 | 20 | 10.386947858 | 10.305616126 | 0.0579942093 | 0.0813317312 |
| 64 | 128 | 19 | 10.369172623 | 10.281995558 | 0.0622895556 | 0.0871770651 |
| 64 | 128 | 20 | 10.380476474 | 10.298674360 | 0.0583392186 | 0.0818021138 |
| 128 | 64 | 20 | 32.769598137 | 32.769596142 | 1.3829333e-6 | 1.9951496e-6 |
| 128 | 128 | 18 | 32.699950881 | 32.699942698 | 5.6715649e-6 | 8.1823154e-6 |
| 128 | 128 | 19 | 32.739116968 | 32.739113643 | 2.3051477e-6 | 3.3256213e-6 |
| 128 | 128 | 20 | 32.759208744 | 32.759206705 | 1.4132063e-6 | 2.0388243e-6 |

Three additional full references check message-size endpoints at t64/s20:

| BCH block | log2 K | Q1 margin | Full margin | Higher occupations / Q1 | Margin loss |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64 | 12 | 16.639437359 | 16.638888171 | 0.0003807405 | 0.0005491879 |
| 128 | 12 | 37.661913405 | 37.661913342 | 4.3832597e-8 | 6.3237069e-8 |
| 128 | 26 | 26.783189820 | 26.783066006 | 8.5825291e-5 | 0.0001238144 |

At the small endpoint, all occupations fit in the explicit calculation:
Q5..128 for BCH-64 and Q5..64 for BCH-128. No dense interval is omitted.
For BCH-128 at the large endpoint, the full union includes Q5..256 and
a complete Q257..1,048,576 dense cover. Its selected 90-digit dense
replays have maximum absolute log error 1.61e-7. The v4 verifier also
replays the dominant Q2, Q3, and Q4 components at 90 digits for every
new reference, and authenticates transitive source dependencies.

The BCH-64 log2 K=26 result is still pending. Its Q513..2,097,152
dense interval has 65.517619946 margin bits. The Q5..512 calculation
uses lazy composition evaluation with a 25-bit search target; every Q
must still be included and the complete union replayed. The earlier
Q5..256 run was explicitly stopped after spending substantial effort
on a 55-bit target. Its retained tilt cache is not a completed interval.
The saved refined Q2..4 checkpoint is required for the eventual full
reference. After this interval closes, use the same sparse refinement
method at the other K anchors before interpreting their relative slopes.

Margins and losses are in bits. BCH-128 has the following complete
occupation decomposition at t64/s20:

| Covered occupations | Contribution margin, bits |
| --- | ---: |
| Q1 | 32.769598 |
| Q2..4 | 52.234612 |
| Q5..64 | 140.960729 |
| Q65..256 | 1601.386989 |
| Q257..16384 | 62.493745 |
| Complete union | 32.769596 |

For BCH-64, Q5..256 has 37.935720241 margin bits, Q257..1024 has
1,528.683863309 bits, and Q1025..32768 has 44.108822637 bits. Q2..4
accounts for essentially all the measured aggregation loss.

`verify_bch_full_reference_v2.py` replays every selected Q2..4 composition,
every bound in the higher intervals, and their aggregation. It checks
complete disjoint dense covers of 733,273,355,136 integer types for
BCH-128 and 5,864,955,809,280 for BCH-64. The new positive coefficient
products agree with the log-domain implementation. Three selected dense
witnesses per block replay at 90 digits, with maximum log errors below
7.11e-9. These are verified binary64 diagnostics, not outward arithmetic
certificates. Other geometries require their own evidence.

At t128/s20, the bulk-spectrum cover described below gives dense margins
of 93.647787109 bits for BCH-64 and 94.484899898 bits for BCH-128.
It covers Q257 through L using 37 boxes over 32,512 integer types for
BCH-64, and 94 boxes over 134,209,152 types for BCH-128. The Q5..256
margins are 37.920582616 and 140.929206629 bits, respectively.
`verify_bch_full_reference_v3.py` replays these complete unions, including
the character-based activation transfer and selected 90-digit witnesses.
Doubling t from 64 to 128 at s20 therefore costs only 0.00694 full-margin
bits for BCH-64 and 0.01039 bits for BCH-128 at this K. This comparison
has complete-tail evidence; the same conclusion does not follow at
smaller s or other K values.

The t128 comparison now extends below s20. Reducing s from 20 to 19
costs 0.01668 full-margin bits for BCH-64 and 0.02009 bits for BCH-128.
For BCH-128, reducing s once more to 18 costs another 0.03917 bits.
At all these tested settings the higher-occupation penalty remains small.
The current full bound for BCH-64 at s18 remains weak, so the apparent
one-bit difference in the available state-size boundary may be slack in
the bound. It is not evidence of an intrinsic difference between the codes.

The exact integer kernel refinement finds 46 positive zero-state
lower-bound witnesses in the 130-point engineering grid. Every positive
witness replays at 90 digits, with maximum log error 1.39e-10. At K=2^20,
the detected first-moment obstructions are:

| BCH block | t | Tested states with a positive first-moment lower bound |
| --- | ---: | --- |
| 64 and 128 | 64 | 7 through 8 |
| 64 and 128 | 128 | 8 through 16 |
| 64 and 128 | 256 | 9 through 20 |

For example, BCH-128 at t256/s20 has a bad-word first moment at least
approximately 2^102883.94 under the exact-kernel diagnostic below. Q1 alone
has a positive margin near 32.75 bits there. Thus these Q1 curves cannot
describe a small complete first-moment bound throughout the tested t/s
range. A nonpositive or unavailable lower bound does not establish that
another setting closes. No failure-probability lower bound is claimed.

The earlier two-weight shortcut found 37 of these settings. Retaining all
kernel coefficients strengthens that screen: in particular, t128 remains
obstructed through s16. Thus the next useful t128 checks start at s17;
this is a necessary restriction from this lower bound, not a sufficient
condition for closure. At t256 every available state setting is obstructed.

`bch_engineering_evidence_v4.csv` labels all 130 engineering geometries:
ten have complete Q1-dominant bounds, 46 have first-moment obstructions,
two have weak full upper bounds, and 72 have only sparse evidence in this
updated analysis. The companion figure `bch_q1_vs_dense_tradeoff_v4.png`
puts the Q1 curves above the
zero-state lower exponents. Its stars identify the complete bounds.
No full-margin surface is interpolated through the unresolved points.
`bch_full_bound_coverage_v4.png` maps the evidence at K=2^20 and prints
the continuous margin loss at every useful full-bound setting. White
cells are outside the recorded grid; gray cells have only sparse evidence.
`bch_full_bound_k_scaling_v4.png` shows full-margin anchors against the
Q1 grid, with their margin losses on a logarithmic scale. It draws no
full-bound interpolation through the missing geometries. The report
uses each full reference's selected Q2..4 components and keeps the old
coarse sparse penalty in a separate column.

Validation: all 116 workstream tests pass, including the finite-field
transfer, positive coefficient, exact integer kernel, lazy composition,
joint witness, activation density, and integer subdivision checks.
The complete reference replays and all
46 positive lower-bound replays are separate from that test suite.
Two additional interval-ledger tests pass for v4, including rejection
of missing, overlapping, reversed, and mismatched dense intervals.

The distinction matters when choosing the state dimension s. Q1 has one
nonzero outer row, so each region contains at most one active input bit.
It includes cancellation of the state by that bit. Higher occupations
also admit multiple input bits within an epoch and nonzero inputs with
zero syndrome. These events can change the apparent state-size knee.

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
