# Checked BCH-sandwich bounds for the dual spectrum

The distance-derived dual kernel left useful BCH structure unused.
`dual_sandwich.py` applies the existing sandwich constraints to selected
dual weights and checks every accepted upper bound over the rationals.
The result improves outer counts, not yet the full SPIN distance proof.

## Model and objective

The existing half-spectrum variables are q_w=A_w(Q) and h_w, the common
weight-w count in a nonzero coset of Q inside P. The fixed constituent
has spectrum q_w+31h_w. These identities use the checked irreducible
eight-dimensional cyclic quotient and Q<=C<=P, not a random choice of C.

For even j, the dual shell is the following linear objective:

    A_j(C^perp) = 2^-128 sum_{w=0,2,...,128} m_w K_j(w)(q_w+31h_w),
    m_w = 2 for w<128, and m_128=1.

The script replays the original constraint audit, including the published
anchor checks and quotient algebra. It removes every Wambach and
orbit-search lower constraint. The independently replayed dual-distance
certificate allows the q and h orthogonal-array equalities through degree
28. Authenticated bounds on q_w+31h_w provide individual upper bounds
q_w<=cap_w and h_w<=cap_w/31.

Zero variables and q_0=1 are eliminated exactly. Remaining variables are
scaled by integers near their binomial sizes; each row is divided by a
positive rational. The resulting model has 91 variables, 332 inequalities,
and 30 equalities. These transformations do not strengthen the model.

## Checking a proposed witness

Write its constraints as Ay<=b, Ey=d, and 0<=y<=U. Include eliminated
variables in the objective constant c_0, so the objective is c_0+c^T y.
For any proposed prices mu>=0 and arbitrary lambda, define

    r_i = max(0, c_i-(A^T mu+E^T lambda)_i).

Then an exact upper is

    c_0 + mu^T b + lambda^T d + r^T U.

The script converts proposed floating prices to exact rational numbers,
checks their signs, computes every residual exactly, and includes its
entire correction. Explicit upper-bound prices can also be included.
Twelve toy tests use deliberately inaccurate prices to check the repair.

An optional reconstruction solves selected active-column equations over
the rationals. It is only another price proposal. The same complete
sign and residual checks apply afterward; no numerical feasibility test
is trusted.

## Results

Ordinary HiGHS solves reported numerical infeasibility. A bounded search
in the LP dual produced some useful prices, but also invalid numerical
objectives. The latter are not proof results. The table lists only the
upper bounds after exact correction, intersected with the prior caps.

| Dual weight | Prior log2 cap | Checked log2 cap | Search price limit |
|---|---:|---:|---:|
| 30 | 52.449927 | 40.739557 | 100 |
| 38 | 57.501531 | 54.464856 | 100000 |
| 56 | 71.643167 | 70.066010 | 1000000 |
| 64 | 80.024612 | 79.323191 | 1000000 |

The price limits restrict witness search only. They are not assumptions
about the code or new LP constraints on its spectrum. Symmetry supplies
the same caps at complementary weights. Solver failure retains the old
cap, so it cannot manufacture an improvement.

These four checked caps feed the existing positive dual-moment identity:

| Four-row union support | Previous log2 count upper | With checked dual shells |
|---|---:|---:|
| 112 | 303.676619 | 303.676619 |
| 128 | 341.782566 | 337.195182 |
| 144 | 358.518115 | 354.738926 |
| 160 | 379.400345 | 376.661446 |
| 176 | 406.861110 | 405.741690 |

The bounds are exact integers; displayed logarithms are rounded. No
selected-point screen establishes complete support or occupancy coverage.

    python -B research/workstreams/permutation_locality/dual_sandwich.py --propagate

The optional cover flag `--dual-sandwich` recomputes and checks the four
witnesses, then enables the dual-moment refinement. It does not change
the default proof path, production encoder, or paper.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --positive-shortening --dual-shortening --dual-sandwich --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --full-feedback 6 --mature-tail 64 --pair-tail --joint-witness --penalties .75 --tilts .032 --probe-supports 128 144 160

## Solver limitation and next step

For 64 active groups, the selected-point screen above gives the following
log2 upper bounds. Each point assigns the same union support to every
active group; the tilt is .032 and the penalty is .75.

| Union support per group | Previous screen | With checked dual shells |
|---|---:|---:|
| 128 | 3046.142612 | 2839.547802 |
| 144 | 2593.176889 | 2408.779214 |
| 160 | 2396.423683 | 2259.108135 |

These are binary64 diagnostics, not outward-rounded certificates.
Although the shell caps improve the screen, all three bounds remain
above one and therefore give no useful failure-probability guarantee.
They neither close the intermediate regime nor show that the code fails.

One 60-second SoPlex attempt used the previously installed runtime under
WSL, rational input, and multiprecision settings. It stopped at its time
limit after repeated numerical failures. It returned no rational dual
solution and supplied no bound. The `--soplex` option remains a diagnostic
only; its output is not accepted as a certificate.

## Dual-coordinate experiment

The original objective extracts a small dual shell by cancellation between
large primal-spectrum terms. `dual_coordinates.py` now implements an exact
change to dual-spectrum coordinates. Define q'=A(P^perp) and
h'=(A(Q^perp)-A(P^perp))/255. Containment makes h' nonnegative: it averages
the spectra of the 255 nonzero P^perp-cosets inside Q^perp.
No assumption that these dual cosets have identical spectra is needed.
Then

    A(C^perp)=q'+7h',
    q = 2^-133 K(q'+255h'),
    h = 2^-133 K(q'-h').

Here K is the full MacWilliams transform, or its folded even half-spectrum
matrix. These formulas follow by applying the transform twice and using
dimensions 123 and 131 of Q and P. The script checks K^2=2^n I exactly
for n=4,8,12,16,256. Exhaustively enumerated repetition/RM sandwiches at
n=8,16 also check the spectrum and intermediate-code identities.
Together these tests check 4303 entries. The source audit and production
generator checks still run before constructing the BCH model.

Every retained original row is transformed over the rationals, including
primal nonnegativity and fixed low-weight coefficients. The model adds
the existing dual kernel caps and eliminates the known dual coefficients
below weight 30. It has 100 variables, 476 inequalities, and 39 equalities.
The objective q'_j+7h'_j has no signed cancellation.

This change did not produce stronger checked bounds in the bounded trial.
HiGHS failed numerically on the primal problems at weights 30,38,56,64.
Direct dual search with price limit 100 failed at the first three weights.
At weight 64 it produced a repaired bound weaker than the prior cap.
A 60-second SoPlex trial at weight 30 timed out after numerical failures,
with no rational solution. None of these solver statuses proves mathematical
infeasibility. The old weight-30 witness reproduces its prior checked bound
after the shared model refactor.

    python -B research/workstreams/permutation_locality/dual_coordinates.py --self-test-only
    python -B research/workstreams/permutation_locality/dual_coordinates.py --weights 30 38 56 64
    python -B research/workstreams/permutation_locality/dual_coordinates.py --direct-dual --weights 30 38 56 64
    python -B research/workstreams/permutation_locality/dual_coordinates.py --soplex --weights 30

A follow-up `--soplex-tight` diagnostic lowered the internal zero,
factorization, update, and pivot tolerances to 1e-80, retaining 100-digit
arithmetic. SoPlex aborted in `betterThreshold` with assertion
`th < R(1.0)`. It returned no rational solution or usable bound.
This is a solver failure, not evidence of mathematical infeasibility.
Neither SoPlex diagnostic contributes a certificate.

Before investing in another solver, test whether substantially smaller dual
shells would fix the intermediate-occupancy bound. `dual_sensitivity.py`
performs this counterfactual test separately from the certificate driver.
It substitutes rounded-up even-binomial counts at allowed dual weights;
these counts are not proved upper bounds for the BCH dual.

The test uses ceil(binomial(256,w)/2^127), intersected with the checked
caps, at even weights 30 through 226. It retains the known endpoints and
zeros. This artificial vector is not asserted to be a realizable spectrum.
All other inputs remain unchanged. It is passed through the same positive
dual-moment and all-one-column weighting formulas, with penalty .75.

| Support per group | Checked weighted log2 count | Counterfactual count | Projected 64-group log2 upper |
|---|---:|---:|---:|
| 128 | 344.094303 | 283.870567 | -1014.771284 |
| 144 | 368.421425 | 319.665579 | -711.594931 |
| 160 | 394.528954 | 367.944352 | +557.693603 |

For a singleton support interval u, the count enters the score only through
64 log2(count[u]). Therefore each projected score equals the earlier score
plus 64 log2(counterfactual_count[u]/checked_count[u]). No new inner solve
is needed for this identity. The calculation retains tilt .032 and the
earlier Bernoulli-conditioning witness; it does not optimize new parameters.

At support 176, the same substitution saves only 335.32 bits in the
64-group score; at support 192 it saves 4.99 bits. These observations
suggest splitting the next effort: stronger outer bounds for intermediate
supports, and tighter coupling between outer classes and feedback
cancellation for larger supports. They do not prove a limit on what other
tilts, penalties, spectrum bounds, or proof methods could achieve.

    python -B research/workstreams/permutation_locality/dual_sensitivity.py

The sensitivity script is not imported by the certificate driver. It emits
no certificate and changes no production configuration or paper claim.
