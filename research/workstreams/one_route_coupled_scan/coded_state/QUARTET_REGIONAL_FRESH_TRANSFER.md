# Carry fresh counts through a complete region

This continues [the fresh-count transfer](QUARTET_FRESH_COUNTS.md), with
the construction and full45%/40-bit goal unchanged. The preceding turn
made progress on clustered count paths. This turn sums every count path
inside specified bands, rather than selecting a few regional words.

The matrix formulas below are deterministic inequalities for the positive
upper transfer. Numerical results remain FLOAT diagnostics. They are not
directed certificates, whole-code margins, or failure-probability lower bounds.

## Do not maximize over the fresh count after every update

The state label is the previous input count a and a bucket for the actual
type1 count r1. The previous note bounds transitions between these labels
uniformly over all five-type histograms in the incoming bucket. Write
M_(ij) for this positive upper-transfer matrix. The destination label j
determines its input count a_j. Restrict a_j to a finite band C.

A first shortcut chooses a positive vector v_a on each count's buckets,
then replaces each count-to-count matrix by its maximum weighted row sum.
It loses too much: for C=392..424 its FLOAT expression is+148295.79.
The following calculation retains the bucket distribution for all sixteen
updates of a region before taking a maximum.

For z>0 define

    T(z)_(ij) = M_(ij) z^(a_j-408).

Choose any strictly positive vector v and set

    rho = max_i (T(z)v)_i/v_i,
    P_(ij) = T(z)_(ij) v_j/(rho v_i).

Then P has row sums at most1. No eigenvector assumption is needed. Power
iteration is only a way to propose v; failure to converge does not invalidate
the maximum-row-ratio inequality. The FLOAT scripts report proposal residuals.

Add a single initial-state row using the true all-zero initial state. If
I(z) is its tilted first-update row, define v_start=I(z)v/rho and use
P_start,j=I(z)_j v_j/(rho v_start). No transition enters this extra label.
This does not grant a reset at later regional boundaries.

Let p_s(i) be the sum of the P-weights of all sixteen-update paths starting
at i with total input-count deviation s=sum_t(a_t-408). A backward DP
computes all these coefficients. It shifts by the destination label before
multiplying by P. Consequently the corresponding coefficient of M^16,
applied to v, is exactly

    rho^16 z^(-s) v_i p_s(i).

This identity is for the upper transfer M, not the actual histogram process.
M already dominates that process coefficientwise. The DP does not assume
independence between state coordinates or between regional endpoints.

## Sum all regional payment allocations

Keep the diagnostic inactive-row count2100 and covered count c=111572.
A region contains8192 positions and6092 active outer rows. If it contains
6528+s allowed inputs, its payment count is436+s, with normalization

    Q_s = binom(6528+s,436+s)/binom(8192,2100).

Include precisely those s with0<=436+s<=2100. Other s have zero weight.
For a positive payment marker y, define

    gamma(y) = max_i sum_s p_s(i) (Q_s/Q_0) (y/z)^s.

The normalized regional operator therefore maps v below
rho^16 y^436 Q_0 gamma(y) v. This inequality composes across all256
regions. Let beta_max majorize the actual terminal factor divided by
the underlying histogram potential and v. Its cost is incurred only once.
The coefficient at total payment c is at most

    v_start beta_max rho^4096 Q_0^256 y^(256*436-c) gamma(y)^256.

Here256*436-c=44. Multiply by binom(8192,2100)b^(-165556) to obtain the
same row-union and geometric-tail expression used in the preceding note.

Every regional payment allocation is included, subject to the input-count
band. The scalar y bounds the global coefficient; it does not replace the
regional binomial factors by Bernoulli probabilities or a tangent estimate.
The initial marker proposal uses log2(y/z)=-log2(6528/436), motivated by
binomial log-concavity. The script evaluates every Q_s and then searches y.
Any positive y is valid independently of whether this search converges.

## Numerical results and what fails

All entries below use the unchanged selected field weights and base potential
from the preceding note. They include all count paths inside the band and
all regional payment allocations. The total-payment coefficient is bounded
at111572. They omit paths leaving the band and other covered/row counts.

| Input-count band | Fresh-count bucket width | FLOAT log2 expression |
| --- | ---: | ---: |
| 400..416 | 8 | -12064.53 |
| 396..420 | 8 | -3152.64 |
| 392..424 | 8 | +9357.30 |
| 392..424 | 4 | +7566.84 |

The first two values used the initial payment marker. Optimizing that marker
in the widest band saves only0.105bits with width8 buckets and0.038bits
with width4 buckets. The remaining loss is not a payment-marker search issue.

For the narrower problem with exactly44 payment435 regions and212 payment436
regions, the width8 regional bound is+8040.45. Twelve regional potential
iterations improve it to+7756.10 but do not close it. The iterations have
not converged; do not interpret their final residual as a proof obstruction.
The following lower diagnostic addresses a different and stronger question.

## Distinguish boundary slack from local-transfer slack

For the chosen positive matrix P, define

    alpha_s=min_i p_s(i),
    beta_min=min_i (terminal_i/v_i),

Here terminal_i is the same bucket terminal UPPER vector used in the
computed expression, already divided by the histogram potential. Its value
at the artificial initial label is immaterial because no path returns there.
Positivity gives a lower bound on the UPPER-TRANSFER expression for every
ordered sequence of regional totals. Summing all balanced assignments gives

    v_start beta_min rho^4096 z^44
    * (Q_-1 alpha_-1)^44 (Q_0 alpha_0)^212
    * binom(256,44) binom(8192,2100) b^(-165556).

A positive value for this lower expression means that no improved boundary
potential can make this particular local-transfer expression close. It does
not lower-bound the true failure probability. Nor does it rule out another
base potential, a sharper local coefficient bound, or another proof method.

For C=392..424 and width8 buckets, its FLOAT log2 value is+5130.99.
The analogous lower bound for the marked sum of all payment allocations
is+6329.74. The regional weighted row ratios span9.259bits; their maximum
occurs at the artificial initial label. Even replacing that maximum by
the minimum does not remove the diagnostic failure.

With width4 buckets, the balanced lower diagnostic is+5529.10 and the
marked lower diagnostic is+6844.21. The weighted row ratios now span only
0.238bits; their maximum is at input count394 and fresh-count bucket36..39.
Thus neither the artificial initial row nor the coarser width8 buckets
explain this failure. These are lower bounds on the chosen upper-transfer
expressions, not on the actual ensemble moment or bad-codeword probability.

These lower evaluations still require directed verification before becoming
certified obstructions. They are sufficient evidence to stop expanding the
same boundary-potential search and inspect the local relaxation instead.

## Next local constraint to inspect

The fresh-count transfer remembers r1 but still maximizes over all splits
of the remaining coordinates between types2,3,4. These types encode different
run ages; their numbers cannot be chosen independently after every update.
In the exact transfer, if v1 selected positions came from type1, then

    outgoing type2 count = v1.

Thus the transition that determines the new fresh count also determines
the next age cohort. Retaining this relation is the next discriminating
test. Begin with exact small transitions and the widest failing band,
rather than allocating a large unverified two-cohort state space. Keep the
outer, routing, inner distribution, and prior certificates unchanged.

## Reproduction and checks

- `quartet_fresh_band.py` preserves the rejected one-update scalar shortcut.
- `quartet_fresh_regional.py` implements the sixteen-step coefficient DP,
  regional normalization, positive-vector bounds, and payment marker.
  It checks64 exact small paths with both constant and nonconstant terminal
  vectors, two exact weighted regional sums, and a non-eigenvector example.
- `quartet_count_band.py` accepts explicit initial/terminal potential values;
  its original defaults are unchanged.
- `quartet_fresh_count_screen.py` caches active-coordinate coefficients by
  their actual degree instead of recomputing them for every(k,v0). The
  preceding rotation result replays identically at-3133.023799672046.
  Its numerical computation falls from7.57s to1.08s; this is not an encoder
  performance claim.

Runs are serial, single-threaded, remote on Peach, with a90s watchdog and
no generated files. The regional screen bounds working arrays by768MiB;
the widest width4 run takes about45s. An early width8 run stopped because
the proposed Perron vector did not converge. The corrected calculation
uses a maximum row ratio for any positive candidate, not an assumed
eigenvalue. No timed-out or failed run is counted as a passing result.
