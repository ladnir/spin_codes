# A kernel-aware bound for every occupancy

The low-occupancy argument now has two general components: an epoch transfer
that includes every kernel weight, and a counting bound for arbitrary numbers
of nonzero outer rows. The formulas apply to all occupancies. Numerical
coverage is a separate question.

## Current coverage

Keep the fixed BCH [256,128] code and RM2Sub t128_s15 construction from
`README.md`. There are L = 8192 rows, n = 256 regions, and N = 2^21 output
bits. The cutoff is H = 209716. The state starts at zero; output precedes
the update; each epoch uses an independent nonzero field multiplier.
All permutations and multipliers are sampled once and shared across messages.

For 1 <= Q <= L, let Z_Q count messages with exactly Q nonzero rows and
output weight at most H. The new computation covers every Q from 4 to 16:

\[
\sum_{Q=4}^{16}\mathbb E[Z_Q]
\le\sum_{Q=4}^{16}U_Q<2^{-162}.
\]

Together with the retained Q1–Q3 bounds, it gives

\[
\Pr\left[\sum_{Q=1}^{16}Z_Q>0\right]
\le\sum_{Q=1}^{16}U_Q<2^{-49}.
\]

The remaining range is Q = 17,...,8192. No numerical certificate for that
range is claimed. In particular, a valid formula for all Q is not yet a
proof of the full distance target.

## Epoch transfer for arbitrary input weight

Write t = 128, s = 15, M = 2^s - 1, and kappa = M/(M-1). Let a_w count
nonzero A-codewords of weight w. Let k_j count weight-j words in ker(B).
The exact map and MacWilliams audits supply both spectra. Define

\[
\beta_j:=\frac{k_j}{\binom tj},\qquad 0\le j\le t.
\]

For X of uniform weight j, beta_j is exactly Pr[B(X)=0]. Let

\[
f_j(w;z):=\frac1{\binom tj}
\sum_v\binom wv\binom{t-w}{j-v}z^{w+j-2v},
\]

with impossible binomial terms omitted. For z in (0,1), define

\[
m_j(z):=\frac1M\sum_w a_w f_j(w;z),\qquad
d_j(z):=\max_{w:a_w>0} f_j(w;z).
\]

Here m_j is the uniform-nonzero-state emission moment. The maximum d_j
bounds the emission moment for every nonzero entering state, including a
state created by deterministic activation. Unlike z^(d_A-j), this bound
remains useful when j exceeds the minimum A distance.

Retain three classes of normalized state laws:

- Z: point mass at zero;
- D: any distribution on nonzero states;
- L: a distribution on nonzero states with point probabilities at most kappa/M.

The classes describe decompositions of measures weighted by z to the
accumulated output weight. Define termination numerators

\[
(h_j^D,h_j^L):=
\begin{cases}
(d_j,m_j),&\beta_j=0,\\
(1-\beta_j,1-\beta_j),&\beta_j>0.
\end{cases}
\]

One valid row-to-column envelope is

\[
T_j(z):=\begin{pmatrix}
\beta_j z^j&(1-\beta_j)z^j&0\\
h_j^D/M&0&d_j\\
\kappa h_j^L/M&0&\kappa m_j
\end{pmatrix}.
\]

To justify the zero row, observe that zero remains zero exactly when B(X)=0;
its output has weight j in either branch. Otherwise it enters class D.

For a nonzero entering state q, condition on q and X. If B(X)=0, the next
state alpha q is uniform nonzero and cannot terminate. If B(X) differs from
zero, termination has probability 1/M, independently of the emitted weight.
Conditional on survival, the next law is uniform outside zero and B(X).
Both kinds of surviving law belong to L. Their output-weighted mixtures
remain in L.

When beta_j=0, the full emission moments bound the terminating mass after
division by M. Otherwise, z to the emitted weight is at most one, giving
the bound (1-beta_j)/M for an arbitrary entering law. The extra kappa in
the L row is conservative. Surviving mass is bounded by the full moment;
duplicating termination mass only enlarges the envelope.

This proves the weighted-measure invariant for every j, including kernel
weights and the endpoints j=0,t. It does not assume that kernel inputs are
rare, nor that they behave like nonkernel inputs.

## Region transfers for all occupancies

There are E = L/t epochs in a region. For 0 <= j <= L, define

\[
R_j(z):=\frac1{\binom Lj}[u^j]
\left(\sum_{a=0}^t\binom ta u^a T_a(z)\right)^E.
\]

An independent region permutation makes the positions of its j ones a
uniform j-subset of [L]. Conditional on the epoch counts, within-epoch
supports are independent and uniform. The binomial factors count this law
exactly, including every same-epoch collision pattern.

At occupancy Q, no region contains more than Q ones. Truncating the region
polynomial after degree Q therefore preserves all required coefficients.
The state continues across region boundaries.

## A maximum-density counting bound

The Q3 calculation paid sums of conditioning costs. For a general occupancy,
we can instead dominate an entire group of outer words pointwise.

Let U_w be the retained deterministic BCH shell caps. Partition the supported
nonzero weights into five groups:

```
38..54, 56..86, 88..128, 130..184, {186,188,...,218,256}.
```

Ranges advance by two. For a group G and 0 < p < 1, define

\[
\Gamma_G(p):=\max_{w\in G}
\frac{U_w}{\binom nw p^w(1-p)^{n-w}}.
\]

To see what this dominates, let nu_G be the unnormalized measure obtained
by summing the shuffled support law over all nonzero C-codewords in G.
For an n-bit support x of weight w in G,

\[
\nu_G(x)=\frac{A_w(C)}{\binom nw}.
\]

Outside G its mass is zero. If mu_p is the product Bernoulli(p) measure,
then the shell caps give the pointwise inequality

\[
\nu_G(x)\le\Gamma_G(p)\mu_p(x)
\quad\text{for every }x\in\{0,1\}^n.
\]

These are counting measures, not two normalized code ensembles. Tensoring
the pointwise inequality for Q independently permuted rows and integrating
any nonnegative function preserves it. We apply it to the conditional
output-weight moment over the remaining region permutations and multipliers.
The original code and its fixed-weight row permutations are unchanged.

## General composition bound

Let c=(c_1,...,c_5) be nonnegative counts with sum Q. Choose a probability
p_g in (0,1) for each group. Define pi_c(j) by

\[
\sum_{j=0}^Q\pi_c(j)v^j
:=\prod_{g=1}^5(1-p_g+p_gv)^{c_g}.
\]

Under the auxiliary product laws, pi_c is the number-of-ones distribution
in each region. Put

\[
M_c(z):=\sum_{j=0}^Q\pi_c(j)R_j(z).
\]

With e_Z=(1,0,0) and the three-entry column 1, a contribution bound is

\[
B_{Q,c}:=\binom LQ\frac{Q!}{\prod_g c_g!}\,
z^{-H}e_ZM_c(z)^n\mathbf1
\prod_g\Gamma_{G_g}(p_g)^{c_g}.
\]

Choose the Q row positions and order them by index. The multinomial factor
counts group assignments to those positions. Local message values may
coincide. The groups partition all nonzero outer words, so

\[
\boxed{\displaystyle
\mathbb E[Z_Q]\le\sum_{c_1+\cdots+c_5=Q}B_{Q,c}}
\qquad(1\le Q\le8192).
\]

Every composition may use its own tilt and auxiliary probabilities. The
pointwise domination argument remains valid for each contribution separately.
Neither weight-deletion monotonicity nor an independent-bit approximation
to the original setup enters this formula.

## Checked batch: Q=4 through Q=16

The batch evaluates every one of the 20,293 compositions in this range.
Binary64 selects witnesses; 256-bit Arb recomputes all moments and maximum
density factors. Density factors are rounded upward to dyadics with 192 bits
of relative precision. Final aggregation is exact rational arithmetic.

| Q | Diagnostic -log2(U_Q) |
|---:|---:|
| 4 | 162.0632 |
| 5 | 213.5637 |
| 6 | 253.1444 |
| 7 | 278.4231 |
| 8 | 333.3178 |
| 9 | 380.8339 |
| 10 | 402.8738 |
| 11 | 413.8230 |
| 12 | 432.5992 |
| 13 | 493.3020 |
| 14 | 548.8335 |
| 15 | 599.6100 |
| 16 | 615.9999 |

Every composition was replayed at 512-bit precision, producing aggregate
bounds no larger than the stored ones. This is a higher-precision replay of
the same implementation, not a second independent derivation of the full
batch. Separate exact tests provide 8,736 weighted-prefix comparisons with
the full GF(16) kernel, covering every epoch weight 0 through 8. These tests
include nonzero kernel inputs of weights four and eight. Further exact tests
cover all region weights on a toy instance, the five-group composition
counts, the auxiliary distributions, and pointwise density domination.

The receipt `generated/general_q4_q16_outward.json` is about 75 KB. It stores
per-occupancy bounds and hashes of the discovery witnesses, not large
coefficient tensors. Inputs and generated results remain ignored by Git.
Read-only checks from this worktree root are:

```powershell
python -B workstreams/bch_rm2sub_bridge/test_general_occupancy.py
python -B workstreams/bch_rm2sub_bridge/general_batch_certificate.py --verify
```

The prior Q1–Q3 source files and certificates were not changed. No file in
the other worktree was modified.

## What is still needed for the full theorem?

The transfer is now general; direct composition enumeration is not scalable.
Even with five groups, occupancy 8192 has
binom(8196,4) = 187,879,147,280,385 compositions. Checking only a few examples
cannot replace a bound on that whole set.

The next step is a certified bound over ranges of compositions, followed by
ranges of occupancies. A possible starting inequality is

\[
\sum_c B_{Q,c}\le\binom{Q+4}{4}\max_c B_{Q,c},
\]

but the maximum itself must be bounded over the complete integer simplex.
The present results do not establish convexity or justify checking vertices
only. For example, the largest recorded terms at Q=6,10,11,16 use mixed
groups, not a single group.

Coarse weight groups may also become too costly in the dense range. They can
be refined without changing the general proof, but that increases the need
for certified bounds over whole ranges. Those are the remaining analytic
and computational obligations; the observed margins through Q=16 are not
an extrapolation certificate.
