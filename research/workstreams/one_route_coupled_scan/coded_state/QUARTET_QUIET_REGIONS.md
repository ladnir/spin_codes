# Long quiet intervals: regional counting diagnostic

2026-10-27. The active goal is whole-code closure at 45% distance and at
least 40 bits of margin. The construction and all frozen certificates remain
unchanged. This note develops bounds for longer zero-state intervals; its
numerical evaluations are FLOAT screens, not new finite certificates.

## What is being counted

Retain the construction and conditional message count in
[the band proof](QUARTET_HISTORY_BAND.md). Setup is sampled once and fixed
for every message. Here every nonzero state has weight 216--232, and each
transition between two nonzero states has 24--40 deaths. Zero-state runs
may have any length, including an initial prefix or terminal suffix.
There is no limit on the number of resets. A *nonquiet update* has at least
one nonzero endpoint state; a quiet update has both endpoints zero.

These families extend the short-reset experiment, but still omit other
nonzero weights and turnover. The three screened inactive-row counts are
1850, 1974, and 2100. No interpolation between those counts is certified.

The question is whether long quiet intervals fail because the construction
is weak or because the scalar continuation bound discards outer constraints.
A quiet update forces 512 input symbols to zero. It also contributes no
local restriction on which of those positions carry inactive-row labels.
Both facts should remain in the count.

## A regional quota bound

Fix an exact inactive-row set of size r. Each region has M=8192 positions,
exactly r of which carry its labels. Use independent Bernoulli-p labels as
a counting device, then condition their total on r. Write

    pi_r(p) = binom(M,r) p^r (1-p)^(M-r).

Fix a regional history with h nonquiet updates. Its other M-512h positions
belong to quiet updates. Their labels do not occur in the local integrand.
For any nonnegative integrand F on the nonquiet positions,

    E[F | regional total=r] <= kappa_h(p,r) E[F],

where

    kappa_h(p,r) = max_s Pr[Bin(M-512h,p)=r-s] / pi_r(p),
    max(0,r-(M-512h)) <= s <= min(r,512h).

To prove this, condition on the nonquiet labels. Their total is s; the
remaining independent labels supply the displayed binomial probability.
Maximize that probability, then average F. In particular kappa_0=1:
a wholly quiet region incurs no label-conditioning penalty. For h=16,
the bound reduces to pi_r(p)^-1.

This removes only a routing-label penalty. All forced input zeros remain
in the shared shortened outer message space. The argument does not sample
a new message or outer code after a reset.

The proof may choose p separately for each regional endpoint class and h.
These are alternate positive upper bounds on the same conditioned count,
not changes to the routing distribution. The state potential must remain
the same across regions for their operators to compose.

## Combine output cancellation with the zero-selection-free count

Let A be total state weight and D total deaths. Put W=943718, q=2^32,
theta=7/40, and retain the positive output marker ell=13592/664659.
Let B_W denote the conditional counting expression obtained by evaluating
the selected-zero polynomials at ell. Before taking its theta moment,
the zero-selection-free expression B_0 satisfies

    B_0 <= B_W (q ell)^(A+D-W).

Every feasible local selected-zero polynomial has constant coefficient one.
Removing the output-zero requests therefore gives this inequality. When
A+D<W, the output-weight restriction is automatic; use B_0 rather than
claim that a negative number of selected zeros is required.

For any proof parameter 0<=tau<=1, combining the two counts gives

    min(B_0,B_W)^theta
      <= B_W^theta sigma^[tau(A+D-W)],
    sigma=(q ell)^theta.

The inequality follows from min(1,x)<=x^tau for x>0. It applies before
averaging the same route and inner randomness. It is not an interpolation
of separately averaged probabilities for different histories.

The transfer for a transition of next weight v and d deaths receives
sigma^[tau(v+d)], and the global prefactor receives sigma^(-tau W).
Use the common potential

    f_tau(w)=X^w ell^(theta w) sigma^(-tau w),
    X=97905/432928,   f_tau(0)=1.

The next-weight factor cancels in the transformed row. The remaining
multiplier is sigma^[tau(u+d)] for incoming weight u. Thus the retained
incoming-weight/death coefficients can be reused for each tau.

## Composing regions

Aggregate states into zero and the live band. Bound the transformed
zero-to-live, live-to-zero, and live-to-live row sums by a, b, and g.
The zero-to-zero entry is one. Grade every other edge by a marker z,
giving the two-state operator

    T(z) = [[1, a z], [b z, g z]].

Extract each coefficient of T(z)^16 and multiply its degree-h entries
by kappa_h. Take the smallest available complete bound per entry and h,
then sum h. Raise the resulting regional matrix to power 256, start at
the zero state, and allow either terminal class. The terminal factor is
1 for zero and 1/min_w f_tau(w) for the band. The original outer prefactor
and binomial row-set union appear once, outside this product.

The cheap g calculation retains incoming weight and deaths jointly but
relaxes the outgoing weight. It uses a 233-by-41 logarithmic coefficient
array. A second diagnostic retains outgoing weights for p=0.6 and the
original p=112627/167924, using the preceding directed coefficient helper.
Subsequent logarithms and matrix operations are still floating point.

## Selected results and interpretation

All entries below are log2 upper-bound *screens* for one specified r,
including its row-set union. Positive values are inconclusive, not
counterexamples. Negative values require directed replay before becoming
finite certificates.

| Bound, unrestricted quiet intervals | r=1850 | r=1974 | r=2100 |
|---|---:|---:|---:|
| Regional quota cap only, tau=0 | +23367 | +69478 | +120196 |
| Also stratify nonquiet count; use A+D<=272h | -6788 | +14236 | +36347 |
| Interpolate counts, tau=0.06 | -11255 | +4880 | +43515 |
| Same, retain outgoing weights at the two largest p values | -17093 | +4857 | +43492 |
| Interpolate counts, tau=0.10 | +9238 | +11380 | +17421 |
| Same outgoing-weight refinement | +8011 | +10163 | +16204 |

The total-nonquiet-count calculation is tightest near h=3470, the crossover
at W/272. That location partly reflects the coarse 272-per-update bound,
not an observed low-weight codeword. Interpolating the two conditional
counts avoids imposing that maximum on every update.

The outgoing-weight refinement helps the first row-count slice, but does
not remove the other gaps. It is therefore not the sole missing ingredient.
The regional cap still maximizes quiet-label occupancy independently of
the live-transition count. That is the next explicit relaxation to remove:
retain the live inactive-label total, convolve it with the quiet binomial,
and impose the exact regional quota before collapsing the regional operator.

Do not redesign the inner based on these positive screens. Also do not
declare the long-quiet family closed from the negative r=1850 screen. Other
state weights, turnover, and outer-row counts remain whole-code obligations.

## Replay and checks

The source-only helper is `quartet_quiet_regions.py`. Its default mode runs
the cheap screens; `--retain-next` runs the two expensive coefficient boxes;
`--tests-only` runs the independent checks. Exact small-instance checks cover
81 regional-conditioning inequalities and 21 interpolation inequalities.
Exhaustive short paths check 56 coefficients, including the nonquiet grading
and both terminal states. All checks pass. The previous excursions replay
also reproduces its frozen band, 25-exception, and short-reset bounds exactly;
the frozen pair-defect helper's SHA256 remains unchanged.

The cheap screen took about 6.5 seconds; the two-box comparison about
31.5 seconds. An initial eight-box comparison exceeded the 60-second guard
without producing a result and was reduced to the two discriminating boxes.
All runs were serial and file-free, with at least 20 GiB free. No native
benchmark, construction/library/paper change, commit, or push was made.
