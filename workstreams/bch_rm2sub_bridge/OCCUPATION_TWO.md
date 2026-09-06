# Occupation two for BCH [256,128] with RM2Sub

The occupation-one integration left open messages with multiple nonzero BCH
rows. This note closes occupation two for the same model with t = 128 and
s = 15. It introduces no independent-bit approximation for the permutations.

## Scope and result

Use the fixed code, setup distribution, and output-before-update recurrence
defined in `README.md`. There are L = 8192 outer rows, 256 regions, N = 2^21
output bits, and cutoff H = 209716. Multipliers are independent nonzero field
elements, sampled once per epoch and shared across messages.

Let Z_2 count messages with exactly two nonzero rows and output weight at most
H. The computed rational upper bound U_2 satisfies

\[
\Pr[Z_2>0]\le\mathbb E[Z_2]\le U_2<2^{-91}.
\]

The diagnostic value of -log2(U_2) is 91.5312913644. Adding the retained
occupation-one bound gives the exact checked inequality

\[
\Pr[Z_1+Z_2>0]\le U_1+U_2<2^{-49}.
\]

These statements cover only occupations one and two. They are not a bound
for messages with three or more nonzero rows.

## Weight-two epoch input

The selected map has distinct nonzero columns B(e_j). Therefore
B(e_j + e_k) is nonzero whenever j differs from k. No epoch of input weight
one or two can have a zero syndrome. This conclusion follows from the exact
map audit, not from an estimated kernel probability.

For j in {0,1,2}, let X_j have uniform support among the j-subsets of [t].
Retain M = 2^s - 1, kappa = M/(M-1), and the three classes Z,D,L from the
README. If a_w counts nonzero words of weight w in im(A), define

\[
m_j(z):=\frac{1}{M\binom tj}
\sum_w a_w\sum_{v=0}^j
\binom wv\binom{t-w}{j-v}z^{w+j-2v}.
\]

Terms with impossible binomial indices are zero. The overlap v counts the
ones of A(q) flipped by X_j. The complete audited A spectrum evaluates this
moment for a uniform nonzero q. In this configuration d_A = 48.

The existing transfers T_0 and T_1 remain unchanged. Add

\[
T_2(z):=\begin{pmatrix}
0&z^2&0\\
z^{d_A-2}/M&0&z^{d_A-2}\\
\kappa m_2(z)/M&0&\kappa m_2(z)
\end{pmatrix}.
\]

The first row follows from nonzero B(X_2): zero activates an arbitrary
nonzero state and emits exactly two bits. An arbitrary nonzero state emits
at least d_A - 2 bits. A near-uniform live law has moment at most kappa m_2.

For either nonzero entering class, condition on q and X_2. The fresh
multiplier makes alpha q uniform nonzero, independently of the emitted
weight. Termination has probability 1/M. Conditional on survival, the next
state is uniform outside zero and B(X_2); it belongs to class L. Mixtures
preserve the density bound, including after weighting by emitted weight.
The surviving entry uses the full moment, which conservatively duplicates
the termination mass. This proves the same weighted-measure invariant used
for occupation one.

## Exact region support law

Fix two nonzero outer words of weights a and b. Their independently permuted
supports S_1 and S_2 are independent uniform a- and b-subsets of [256].
Region i contains j_i = 1_{i in S_1} + 1_{i in S_2} ones.

When j_i = 2, the ones initially occupy different outer-row positions. The
region permutation sends them to a uniform two-subset of its L positions.
They cannot collide or cancel. Conditioned on the numbers of ones per epoch,
the within-epoch supports are independent and uniform. These facts justify
the following region envelopes, with E = L/t:

\[
R_j(z):=\frac{1}{\binom Lj}[u^j]
\left(T_0(z)+t uT_1(z)+\binom t2u^2T_2(z)\right)^E,
\qquad j\in\{0,1,2\}.
\]

Truncating the polynomial after degree two preserves these coefficients.
The coefficient of u^2 includes both possibilities: two ones in one epoch,
or one in each of two distinct epochs. The binomial factors count actual
subsets, without replacement.

## Two-row support count

Let e_Z = (1,0,0), and let 1 be the three-entry column of ones. Define

\[
F_{a,b}(z):=\frac{[x^ay^b]
e_Z\left(R_0+(x+y)R_1+xyR_2\right)^{256}\mathbf1}
{\binom{256}a\binom{256}b}.
\]

The four terms append one outer coordinate with support pattern 00,10,01,
or 11. Every ordered pair of supports appears once. Transfer multiplication
also retains the state across region boundaries. Thus

\[
\Pr[\operatorname{wt}(Y)\le H\mid a,b]
\le p_{a,b}(z):=\min\{1,z^{-H}F_{a,b}(z)\}.
\]

Write A_w(C) for the fixed BCH spectrum, and let U_w be its retained
deterministic shell caps. Choose two row positions and order them by index.
The two local nonzero messages may coincide; their values are otherwise
unrestricted. This gives

\[
\mathbb E[Z_2]
\le\binom L2\sum_{a,b}A_a(C)A_b(C)p_{a,b}(z)
\le\binom L2\sum_{a,b}U_aU_b p_{a,b}(z).
\]

Each sum covers weights 38,40,...,218 and 256. No LP objective from the old
inner is needed here; the deterministic BCH caps alone suffice.

## Arithmetic and checks

The certificate uses the fixed tilt z = exp(-exp(-75/10)). Epoch and region
entries are evaluated with 256-bit Arb and converted to binary64 upper
bounds. The outer support recurrence uses nonnegative sums and products,
rounding every operation upward. Positive underflows are replaced by the
smallest positive subnormal. Exact binomial reciprocals are enclosed by Arb.
All observed entries remain finite. The final aggregation treats every
binary64 upper bound as an exact rational number.

The certificate is `generated/t128_s15_q2_j-75.json`, with a small companion
NPZ file containing the region and pair bounds. Both are ignored by Git.
Their provenance includes the scripts, map manifest, arithmetic helper,
deterministic shell caps, and Johnson certificates used in aggregation.

Completed checks:

- 3,840 exact rational prefix comparisons against a full GF(16) kernel, with
  epoch weights zero, one, and two and all required entering-law classes.
- Explicit enumeration of zero-, one-, and two-position subsets in a small
  region, including same-epoch and different-epoch placements.
- All 25 outer weight pairs on a four-region toy instance, compared with
  exact rational enumeration, plus a positive-underflow test.
- Both zero-weight axes of the actual 257-by-257 pair table, compared with
  the independent occupation-one Arb recurrence.
- Bit-for-bit replay of the complete directed pair table and exact rational
  reaggregation of every supported BCH weight pair.
- An independent 512-bit Arb calculation for the dominant pair (38,38),
  using binary region powering and normalized support probabilities. Its
  value is below the stored bound, which exceeds it by about 2.01e-13
  relative. This independent high-precision check covers that pair, not the
  entire table.
- Exact comparison of U_1 + U_2 with 2^-49.

Reproduce from this worktree root:

```powershell
python -B workstreams/bch_rm2sub_bridge/test_occupation_two.py
python -B workstreams/bch_rm2sub_bridge/check_q2_dominant.py
```

The producer is `occupation_two.py`; it refuses to overwrite existing
certificate files. All changes remain in this integration directory, and
the other worktree has not been modified.

Next, derive occupation three with the same support law and state classes.
The audited t128_s15 kernel has minimum weight four, so weight-three epoch
inputs also cannot have zero syndrome. Kernel events first become possible
at occupation four and must be included there.
