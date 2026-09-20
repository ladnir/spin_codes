# Selected (128,19) ladder through K=2^24

All three requested rungs now have full distance/setup margins above 40 bits
for the fixed BCH-256 / selected t128_s19 construction. Each rung has complete
occupancy coverage, numerical verification, and an exact full union. The setup
and distance event remain those in
[CURRENT_UNDERSTANDING.md](CURRENT_UNDERSTANDING.md).

## Certified checkpoint

| log2(K) | Full certified margin | New Q1-only margin |
|---:|---:|---:|
| 16 | 53.9443672720 | Existing certificate |
| 18 | 52.3463883689 | Existing certificate |
| 20 | 50.4482033129 | 50.4482042479 |
| 22 | 48.4706838718 | 48.4706841092 |
| 24 | 46.4762221823 | 46.4762222458 |

The three new Q1 bounds have 256-bit producers and 512-bit replays using
linear epoch iteration. They use freshly computed coefficients and the exact
BCH weighted inequality. Q1 is only one summand of the required full union.
None of these Q1 bounds alone establishes full closure or the true margin.
The full certificates additionally cover every higher occupancy; see
`T128_S19_M20_CLOSURE.md`, `T128_S19_M22_M24_CLOSURE.md`, and the ledgers
`generated/t128_s19_m{20,22,24}_ladder_full_v1.json`.

Records: `generated/t128_s19_m{20,22,24}_ladder_q1_v1.json` and matching
`_replay.json` files. `ladder_instance.py` extends the authenticated instance
identity through m24 without changing the frozen m<=20 interface.
`ladder_q1.py` implements the size-specific calculation. Four tests passed,
including agreement with the frozen Q1 backend and preservation of existing
instance identities.

## Scaling decision

At K=2^24, L=K/128=131072. Enumerating every (Q,h) all-one-count case would
require billions of cases. Partitioning arrays bounds memory but does not
remove that quadratic work. The completed ladder instead uses interval
certificates whose matrix dimension does not grow with L.

The older t64_s20 proof already uses a shuffled Poisson-binomial comparison
and continuous-density intervals. Its numerical bounds do not transfer to
t128_s19. The reusable ingredients are the comparison proof in
`K30_CLOSURE_PROGRESS.md` and the interval arithmetic in `frontier_dense.py`.
The new calculation must use the selected (128,19) four-state transfer.

## Candidate dense bound: explicitly sum the band labels

This section records the initial candidate. It is not by itself a numerical
certificate. The completed calculation uses its two-tilt extension below.

Fix Q occupied rows among L rows and a positive output tilt lambda. Use the
12 ordinary weight bands, each with reference bit probability p_g in (0,1)
and exact row-density factor Gamma_g. Add one separate all-one label with
p_*=1 and Gamma_*=1. Zero rows are outside the occupied set.

For an ordered label assignment a=(a_1,...,a_Q), define

\[
\nu_a=Q^{-1}\sum_{i=1}^Q p_{a_i},\qquad
G_a=\prod_{i=1}^Q\Gamma_{a_i}.
\]

The reference density in a region is r=Q nu_a/L. The shuffled
Poisson-binomial comparison gives a factor (L+1)^256 relative to independent
Bernoulli(r) bits. Let T(r,z) be the four-state epoch transfer under those
independent bits, with z=exp(-lambda). Define

\[
M_Q(\nu,\lambda)=e^{\lambda H}
 e_Z^{\mathsf T}T(Q\nu/L,e^{-\lambda})^{256L/128}\mathbf1.
\]

For a density interval I, choose a real auxiliary slope eta and set

\[
S(\eta)=e^{-\eta}+\sum_{g=1}^{12}\Gamma_g e^{-\eta p_g}.
\]

The first term is exactly the all-one label's contribution. Expanding the
Qth power sums all ordered label assignments, including every all-one count.
Consequently,

\[
\sum_{\boldsymbol a:\nu_{\boldsymbol a}\in I}
 G_{\boldsymbol a}M_Q(\nu_{\boldsymbol a},\lambda)
\le S(\eta)^Q\sup_{\nu\in I}
 \{e^{\eta Q\nu}M_Q(\nu,\lambda)\}.
\]

To see this, bound each moment by the displayed supremum times
exp(-eta Q nu_a), then enlarge the remaining nonnegative label sum to all
assignments. Its value is S(eta)^Q by the product expansion.

Multiply by binomial(L,Q)(L+1)^256 to obtain a bound for messages assigned
to this density interval. Summing bounds for intervals covering every feasible
nu covers every label assignment. Shared interval endpoints may overcount;
they do not omit assignments. Each interval may use its own fixed eta, lambda,
and associated verified calculation. Fix the row probabilities across an
interval cover so that every assignment has a consistently defined nu.

Unlike a pooled 13^Q maximum, this expression counts the all-one label as
one unit term. It does not require enumerating h, nor does it put the all-one
word into an ordinary band envelope. No sampled density can stand in for the
supremum: the intended checker uses outward interval bounds on T and the
scalar exponent. All accepted interval contributions must be summed exactly.

## Dense range closed at K=2^20

Every occupancy Q=512,...,8192 now has a combined certified upper bound
of approximately 2^-577.1425238434. This is a dense-range contribution, not
the full margin. The certificate is
`generated/t128_s19_m20_ladder_dense_refined_v1.json`, with its exact replay
in the adjacent `_replay.json` file.

The source cover has 1,011 rectangles. Every retained upper bound passed
a 512-bit replay, and its binary partition tree covers the full Q,nu domain.
Source: `generated/t128_s19_m20_ladder_fixed_reference_v2/cover_0000.json`;
numerical replay: `replay_cover_0000.json` in that directory.

The source calculation used the loose factor (L+1)^256 and consequently
labelled 780 leaves unresolved. Their bounds were still valid and finite.
The final certificate replaces that factor by D_L^256, where
D_L=min(L+1,4 ceil(sqrt(L)))=364. Multiplying the complete source union by
the exact rational (364/8193)^256 closes the dense range. The refinement
does not discard unresolved leaves and does not require their old search
thresholds to pass. `POISSON_DENSITY_REFINEMENT.md` proves the smaller factor;
`refine_ladder_density_cover.py` checks the exact conversion and its receipts.

### What made the scalable calculation work

The mean-only label-sum bound stayed weak at some sampled densities, even
after varying its row probabilities and using singleton weight bands.
A bounded five-category fixed-type search also remained weak. Neither result
established a code obstruction.

The successful bound adds an input-count tilt to the output-weight tilt.
It then holds the matrix's reference probability fixed within a box and
allows the input tilt to depend on the box coordinate. This requires an
explicit factor binomial(Q+12,12) for the possible band compositions.
The scalar variation is bounded by a convex secant and concave tangents.
The derivation, including the exact all-one factor and the extra composition
count, is in `DENSE_TWO_TILT.md`.

The first subdivision policy repeatedly split a tiny density interval while
leaving occupancy wide. Balanced subdivision removed that stall. A coarse
tilt grid and repeated subdivision of weak singleton occupancies caused a
second stall; the retained version uses finer witnesses and records tiny
unresolved boxes without endless subdivision. An earlier fixed-reference
run exited without a receipt after a slope proposal became nonfinite.
The retained checker uses an arbitrary valid slope in that case and verifies
its endpoint majorization outwardly. The complete retained cover, including
every parked box, passed numerical replay before the density refinement.

`ladder_dense_labels.py`, `ladder_dense_tilt.py`, and the earlier search
drivers remain intact as bounded alternatives and diagnostic history.
The certificate uses `ladder_dense_fixed_reference.py` through
`search_ladder_dense_fixed.py`. Cached four-state matrices and finite
rectangle covers replace the quadratic (Q,h) enumeration.

## Larger rungs and completion audit

K=2^22 and K=2^24 are complete, with freshly recomputed sparse and dense
bounds. Their dense covers contain 241 and 245 rectangles respectively,
each numerically replayed at 512 bits. In both cases Q1..511 are covered by
size-specific sparse workers with matching 512-bit replays. The exact full
unions have margins 48.47068387175386 and 46.476222182258425 bits.
`T128_S19_M22_M24_CLOSURE.md` records the full scope, receipt hashes, and
reproduction commands. Q1 dominates every rung.

`audit_ladder_completion.py` checks all three rungs together, including a
fresh 768-bit evaluation of every dense box, exact component reconstruction,
instance agreement, and preservation of earlier ledgers. Its retained output
is `generated/t128_s19_ladder_completion_audit_v1.json`.

Next: compare these proved choices against the completed performance
measurements, then return to the heuristic spectrum/growth-curve track.
The certificates cover the three specified sizes; they do not claim a
continuous curve of true failure probabilities. Preserve the hash-bound
numerical sources and proof note when extending the work.
