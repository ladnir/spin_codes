# Stronger shell caps and the next occupancy range

The full-range goal remains the fixed BCH [256,128] outer with RM2Sub
t128_s15 at message length 2^20. The construction and setup experiment are
those in `README.md`. In particular, each epoch uses a fresh independent
nonzero field multiplier, and output precedes the state update.

Two changes strengthen the analysis without changing that construction.
The first bounds termination mass more tightly. The second uses the already
certified dual distance to bound BCH shells that previously had loose caps.

## Verified checkpoint: every occupancy from 1 to 128

Let Z_Q count messages with exactly Q nonzero outer rows and output weight
at most H = 209716. Probability and expectation for Z_Q concern the random
setup, shared by all messages. The new receipt proves

\[
\sum_{Q=65}^{128}\mathbb E[Z_Q]
\le \sum_{Q=65}^{128}U_Q < 2^{-148}.
\]

The aggregate diagnostic margin is 148.8632719886 bits. Every integer in
this interval was computed with 256-bit Arb and replayed at 512 bits.
The replayed bounds were no larger than the stored bounds.

The previous Q1--Q64 certificates are retained. Exact rational aggregation
in `verify_coverage_128.py` checks that the complete partial sum remains
below 2^-49. This is not a statement about occupancies above 128.

The receipt is `generated/tightened_q65_q128_extension_outward.json`.
The earlier coarse screen for Q65--Q128 is also retained; it does not
establish the target throughout that interval.

## Intersect two bounds on termination mass

Use the three state-law classes, spectra, and emission moments from
`GENERAL_OCCUPANCIES.md`. Write M = 2^15-1 and kappa = M/(M-1).
For input weight j, beta_j is the probability of zero syndrome.
The full emission moments are bounded by d_j for an arbitrary nonzero
entering state and by kappa m_j for a law in class L.

Fix a nonzero entering state q and an input X. The next state is zero
with probability 1/M if B(X) is nonzero, and with probability zero otherwise.
The output is emitted before this multiplier-dependent update. Therefore,
for any distribution on the entering state independent of the fresh input,
its weighted termination mass is

\[
\frac1M\mathbb E\left[z^{\operatorname{wt}(X+Aq)}
                 \mathbf1_{B(X)\ne0}\right].
\]

The expectation is at most both the full emission moment and 1-beta_j.
Thus valid termination entries are

\[
\widetilde T_j[D,Z]
  :=\frac{\min(d_j,1-\beta_j)}M,
\qquad
\widetilde T_j[L,Z]
  :=\frac{\kappa\min(m_j,1-\beta_j)}M.
\]

The extra kappa on the second event-probability bound is conservative.
All other entries remain those of T_j. These minima intersect bounds on
the same nonnegative mass; they do not combine conditional laws.
The surviving weighted measures retain the original class-L invariant.
Consequently, positive products and the adaptive group argument remain valid.

`tightened_occupancy.py` implements this envelope separately from the frozen
implementation. Its Arb path takes minima of outward upper endpoints.
Exact tests repeat the 8,736 weighted-prefix checks, including nonzero kernel
inputs, and check every region weight on the small instance.

## A shell cap from the certified outer moments

Here the auxiliary random object is a uniform codeword of the fixed code C,
not a random setup or a random BCH code. The retained construction audit
establishes Q contained in C, where Q is the dimension-123 BCH subcode.
The independently replayed rank certificate proves d(Q dual) >= 30.
Hence C dual is contained in Q dual and has minimum distance at least 30.

For a uniform C-codeword, let W be its weight. Every set of at most 29
coordinates is uniform: failure of surjectivity onto such coordinates would
give a nonzero dual word supported there. Thus, for every polynomial f of
degree at most 29,

\[
\mathbb E[f(W)] = \mathbb E[f(B)],
\qquad B\sim\operatorname{Bin}(256,1/2).
\]

One justification is to express f in the basis of falling factorials.
Each factorial moment is a sum over products of distinct coordinate bits,
whose expectations are fixed by the stated independence.

Define the binary Krawtchouk polynomials by

\[
K_j(x):=\sum_i(-1)^i\binom xi\binom{256-x}{j-i}.
\]

Their binomial orthogonality is
E[K_i(B)K_j(B)] = 1_{i=j} binom(256,j).
For a target weight w, set

\[
P_w(x):=\sum_{j=0}^{14}
          \frac{K_j(w)K_j(x)}{\binom{256}j},
\qquad
S_w:=P_w(w)=\sum_{j=0}^{14}\frac{K_j(w)^2}{\binom{256}j}.
\]

The polynomial P_w squared has degree at most 28. Moment matching and
orthogonality give E[P_w(W)^2] = S_w. Its nonnegative contribution on W=w
is Pr[W=w] S_w^2. Since S_w >= 1,

\[
\boxed{A_w(C)\le \frac{2^{128}}{S_w}.}
\]

The implementation evaluates S_w exactly, rounds this bound down to an
integer, and takes its minimum with the retained deterministic shell cap.
No estimate of the unknown full spectrum enters this calculation.

| Weight | Previous cap, log2 | Strengthened cap, log2 |
|---:|---:|---:|
| 54 | 95.276 | 77.122 |
| 62 | 114.828 | 82.733 |
| 70 | 128.000 | 89.550 |
| 78 | 128.000 | 98.202 |
| 90 | 128.000 | 110.271 |
| 128 | 128.000 | 126.335 |

These logarithms are diagnostics; the certificate stores exact integer caps
and rational pre-rounding bounds. The low-shell caps are retained whenever
they are stronger. Complement symmetry is preserved.

`christoffel_caps.py` replays the 18-leaf rank certificate before producing
`generated/christoffel_oa29_caps.json`. Separate tests check Krawtchouk
orthogonality exactly, the reproducing-polynomial calculation, and shell
domination for a small even-parity code.

## Remaining work

The new shell caps are available for all weights but do not by themselves
close any new occupancy interval. Floating screens above 128 remain
discovery calculations. They must be recomputed outward and aggregated
before extending the coverage ledger.

The next witness search chooses a separate Bernoulli probability for each
weight group. Pure-group calculations select these probabilities only.
The final adaptive maximum still covers every mixed group sequence and
still pays the full group-assignment factor. Checking pure groups alone
would not justify a claim about all messages.

No frozen certificate was replaced. All new work is confined to this
integration directory; the other worktree remains read-only.
