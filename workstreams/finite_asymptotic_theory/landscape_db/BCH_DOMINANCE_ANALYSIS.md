# Checking whether Q1 represents the BCH failure bound

The Q1 engineering curves do not yet establish the full failure margin.
The present task checks BCH lengths 64 and 128, using the exact spectra
and nested inner maps from the engineering study. BCH-256 remains outside
this primary analysis.

Complete reference bounds are verified for both BCH sizes at K=2^20,
t=64, s=20. Q1 dominates the bound at these two geometries. The verified
ratios compare upper-bound contributions, not actual event probabilities.

| BCH block | Q1 margin | Full margin | Higher occupations / Q1 | Margin loss |
| --- | ---: | ---: | ---: | ---: |
| 64 | 10.386947858 | 10.305616126 | 0.0579942093 | 0.0813317312 |
| 128 | 32.769598137 | 32.769596142 | 1.3829333e-6 | 1.9951496e-6 |

Margins and losses are in bits. BCH-128 has the following complete
occupation decomposition:

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

`bch_engineering_evidence_v2.csv` labels all 130 engineering geometries:
two have complete Q1-dominant bounds, 46 have first-moment obstructions,
and 82 have only sparse evidence in this updated analysis. The companion
figure `bch_q1_vs_dense_tradeoff_v2.png` puts the Q1 curves above the
zero-state lower exponents. Its stars identify the two complete bounds.
No full-margin surface is interpolated through the unresolved points.

Validation: all 102 workstream tests pass, including the new finite-field
transfer, positive coefficient, exact integer kernel, lazy composition,
and direct dense-witness checks. The complete reference replays and all
46 positive lower-bound replays are separate from that test suite.

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
