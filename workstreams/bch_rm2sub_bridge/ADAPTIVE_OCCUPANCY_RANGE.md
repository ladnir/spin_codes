# Bounding all group mixtures with an adaptive process

The general formula previously required a sum over weight-group compositions.
An entrywise adaptive bound removes that enumeration. The adaptive choices
belong only to the bounding process; the RM2Sub construction is unchanged.

## Coverage and model

Use the fixed BCH [256,128] code, RM2Sub t128_s15, and random setup from
`README.md`. There are L = 8192 rows, n = 256 regions, N = 2^21 output bits,
and cutoff H = 209716. Initial state, output chronology, independent field
multipliers, and independent permutations remain as previously defined.

For occupancy Q, let Z_Q count bad messages with exactly Q nonzero rows.
The new bounds cover every Q = 17,...,64 and satisfy

\[
\sum_{Q=17}^{64}\mathbb E[Z_Q]
\le\sum_{Q=17}^{64}U_Q<2^{-142}.
\]

With the prior certificates,

\[
\Pr\left[\sum_{Q=1}^{64}Z_Q>0\right]
\le\sum_{Q=1}^{64}U_Q<2^{-49}.
\]

The remaining range is Q = 65,...,8192. This is still not the full distance
theorem or a certificate for seeded multiplier generation.

## Distribute the counting cost across regions

Retain the five weight groups G_g, Bernoulli probabilities p_g, and maximum
density bounds Gamma_g = Gamma_{G_g}(p_g) from `GENERAL_OCCUPANCIES.md`.
For each g choose a positive number rho_g satisfying

\[
\rho_g^n\ge\Gamma_g.
\]

Fix an arbitrary ordered sequence of group labels g_1,...,g_Q. Let J be the
sum of Q independent Bernoulli variables with probabilities p_{g_i}, and
let R_j be the kernel-aware region envelopes. Define

\[
M_{\boldsymbol g}:=\mathbb E[R_J],\qquad
W_{\boldsymbol g}:=\left(\prod_{i=1}^Q\rho_{g_i}\right)M_{\boldsymbol g}.
\]

The scalar factors commute with matrix multiplication. Consequently,

\[
\left(\prod_i\Gamma_{g_i}\right)
e_ZM_{\boldsymbol g}^{\,n}\mathbf1
\le e_ZW_{\boldsymbol g}^{\,n}\mathbf1.
\]

Taking nth roots has distributed each row's counting cost across the n region
matrices. It has not changed the underlying counting measure or introduced
new randomness into the construction.

## Adaptive entrywise envelope

For a fixed target Q, set K_j^(0) = R_j for 0 <= j <= Q. Define recursively

\[
K_j^{(r+1)}[a,b]
:=\max_{1\le g\le5}
\rho_g\left((1-p_g)K_j^{(r)}[a,b]
+p_gK_{j+1}^{(r)}[a,b]\right),
\]

for 0 <= r < Q and 0 <= j <= Q-r-1. The indices a,b refer to the three
state-law classes. The maximum is taken separately for each matrix entry.

For every fixed ordered group sequence, induction on its length gives

\[
W_{\boldsymbol g}\le K_0^{(Q)}
\quad\text{entrywise}.
\]

For the induction step, condition on the next Bernoulli variable. Its two
outcomes yield a positive combination of the previous bounds at j and j+1,
weighted by the fixed group's rho. Maximizing over all groups can only
increase that entry. The base case is R_j itself.

The maximizing group may depend on r, j, and the matrix entry. Thus the
bounding process is stronger than any fixed group assignment. We do not
claim that K is an actual Markov transition matrix or that its choices can
be made by the original encoder. Only entrywise domination is needed.

All matrices are nonnegative. If A <= B entrywise, then A^n <= B^n by
repeated positive multiplication. Hence every fixed group sequence obeys

\[
\left(\prod_i\Gamma_{g_i}\right)
e_ZM_{\boldsymbol g}^{\,n}\mathbf1
\le e_Z\left(K_0^{(Q)}\right)^n\mathbf1.
\]

This step covers mixed compositions as well as pure groups. It does not
assume convexity in the composition or reduce a maximum to simplex vertices.

## Sum all group assignments

There are 5^Q ordered assignments of the five disjoint groups to Q selected
row positions. The previously established counting-measure domination and
Chernoff inequality now give

\[
\boxed{\displaystyle
\mathbb E[Z_Q]\le
\binom LQ5^Qz^{-H}
e_Z\left(K_0^{(Q)}\right)^n\mathbf1.}
\]

Any admissible p_g, rho_g, and z in (0,1) give a valid bound. Different
occupancies can use different witnesses. The price of 5^Q is explicit;
it is not omitted when replacing the composition sum by a maximum.

Once R_0,...,R_Q are available, the adaptive recurrence requires O(5Q^2)
scalar updates per matrix entry and linear working storage in Q. It does
not enumerate group compositions. For shared p_g and z, the same recurrence
can yield bounds for several occupancies by reading K_0^(q) at each depth.

## Witness search and verification

The first search used a coarse common log-odds grid. Some resulting bounds
were vacuous, including at Q=64. Those values are retained in the diagnostic
`adaptive_q1_q64_screen.json`; they are not the successful certificate.

Refinement optimizes the common log-odds shift for each occupancy and tilt.
The tilt grid is z = exp(-exp(j/10)) for j = -60,...,-30. Search arithmetic
is binary64 only; its stored probabilities are then treated as exact dyadic
rationals. The certificate does not require the search to find a global
optimum.

The outward calculation uses 256-bit Arb. It obtains upper bounds on Gamma,
then encloses exp(log(Gamma)/256) upward to form rho. Every candidate in an
adaptive maximum is replaced by its exact upper endpoint before comparison.
All remaining operations are positive matrix products and exact rational
aggregation. The region calculations retain all zero-syndrome branches.

Selected diagnostic margins are:

| Q | -log2(U_Q) |
|---:|---:|
| 17 | 142.4819 |
| 24 | 153.0809 |
| 32 | 157.0517 |
| 40 | 156.6527 |
| 48 | 148.3830 |
| 56 | 154.6816 |
| 64 | 150.5627 |

Every intervening integer occupancy is also covered. The added range's
aggregate diagnostic margin is 142.1575542586 bits.

All 48 occupancies were replayed at 512 bits, producing bounds no larger
than the stored bounds. This is a higher-precision replay of the same
implementation. Separate exact rational tests check entrywise domination
against every one of 120 fixed group sequences on toy instances. Further
tests check the root-weighted counting inequality, the assignment factor,
normalized region coefficients, and the exact range sum. The earlier
all-weight finite-field tests continue to justify the kernel-aware transfer.

The roughly 68 KB receipt is `generated/adaptive_q17_q64_outward.json`.
Read-only checks from the worktree root are:

```powershell
python -B workstreams/bch_rm2sub_bridge/test_adaptive_range.py
python -B workstreams/bch_rm2sub_bridge/certify_adaptive_range.py --verify
```

No earlier certificate or source file was modified. New code and notes
remain in our integration directory; generated artifacts stay ignored.

## Next step

Extend the adaptive method above Q=64 and reuse witnesses across longer
occupancy intervals where possible. This addresses the previous composition
enumeration bottleneck, but does not guarantee that the numerical bound
stays strong in the dense range. The factor 5^Q and the freedom granted to
the adaptive process may eventually be too costly. If so, restrict that
freedom or refine the weight groups while preserving the proved domination.
