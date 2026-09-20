# Bit transpose followed by one accumulator

## Conclusion

A pure bit transpose cannot replace the uniform interleaver in rate-half
Accumulator SPIN.  With probability tending to one, the resulting code has a
nonzero word of weight at most \(B_N\).  Therefore

\[
  B_N=O(\log N)
  \quad\Longrightarrow\quad
  \frac{d_{\min}}N\longrightarrow0
\]

in probability.  No choice of the logarithmic block constant repairs this
obstruction.

The obstruction applies to a deterministic transpose that preserves the same
two row positions in every column.  It does not apply automatically when each
transposed column receives an independent region permutation.

## Pure-transpose ensemble

Let \(B\) be even, let \(N=LB\), and set \(k:=B/2\).  For each
\(i\in\{1,\ldots,L\}\), sample an independent uniform linear injection

\[
  C_i:\mathbb F_2^k\longrightarrow\mathbb F_2^B.
\]

Write the outer word as an \(L\)-by-\(B\) binary matrix.  Row \(i\) is
\(C_i(x_i)\).  The deterministic transpose serializes the matrix in
column-major order:

\[
  (1,1),(2,1),\ldots,(L,1),(1,2),\ldots,(L,B).
\]

The inner map is the length-\(N\) accumulator.  The probability statements
below are over \((C_1,ldots,C_L)\).  The transpose and accumulator are
deterministic.

Let

\[
  U_i:=\operatorname{im}(C_i)\subseteq\mathbb F_2^B.
\]

Each \(U_i\) is an independent uniform \(k\)-dimensional subspace.

## Adjacent-row cancellation

**Lemma 1.** Suppose \(v\in U_i\cap U_{i+1}\) is nonzero.  Then the complete
code contains a nonzero word of weight \(\operatorname{wt}(v)\le B\).

*Proof.* Choose local messages whose encoded rows equal \(v\) in blocks
\(i\) and \(i+1\).  Set every other local message to zero.

Fix a column \(j\).  If \(v_j=0\), the column contains no accumulator input
ones.  If \(v_j=1\), the column contains two consecutive ones at rows \(i\)
and \(i+1\).  The first one changes the accumulator state from zero to one.
The second changes it back to zero.  Exactly one output bit is one.

Every column has even input parity, so the accumulator enters each new column
in state zero.  The total output weight is therefore

\[
  \sum_{j=1}^B v_j=\operatorname{wt}(v)\le B.
\]

The input message is nonzero because \(v\ne0\). \(\square\)

The proof uses only persistent adjacency.  The two rows must remain adjacent
in every serialized column.

## Intersection probability

Fix one \(k\)-dimensional subspace \(U\subseteq\mathbb F_2^{2k}\).  The
number of \(k\)-dimensional complements of \(U\) is \(2^{k^2}\).  The total
number of \(k\)-dimensional subspaces is the Gaussian binomial coefficient

\[
  {2k\brack k}_2.
\]

Hence two independent uniform half-dimensional subspaces have trivial
intersection with probability

\[
  p_k
  :=
  \Pr[U_1\cap U_2=\{0\}]
  =
  \frac{2^{k^2}}{{2k\brack k}_2}
  =
  \frac{\displaystyle\prod_{r=1}^k(1-2^{-r})}
       {\displaystyle\prod_{r=k+1}^{2k}(1-2^{-r})}.
  \tag{1}
\]

The sequence converges to

\[
  p_\infty
  =\prod_{r=1}^{\infty}(1-2^{-r})
  =0.2887880950\ldots.
  \tag{2}
\]

In particular, \(p_k\le2/3\) for every \(k\ge1\).  Thus each adjacent pair
intersects nontrivially with probability at least \(1/3\).

## Vanishing relative distance

Partition the outer blocks into the disjoint adjacent pairs

\[
  (1,2),(3,4),\ldots.
\]

Their intersection events are independent because the pairs use disjoint
outer-code samples.

**Theorem 1 (pure-transpose obstruction).** For the ensemble above,

\[
  \Pr[d_{\min}>B]
  \le
  p_{B/2}^{\lfloor L/2\rfloor}
  \le
  \left(\frac23\right)^{\lfloor L/2\rfloor}.
  \tag{3}
\]

Consequently, if \(L\to\infty\), then

\[
  \Pr[d_{\min}\le B]\longrightarrow1.
  \tag{4}
\]

For every schedule \(B_N=o(N)\), including \(B_N=c\log_2N\), the relative
minimum distance converges to zero in probability.

*Proof.* If a disjoint adjacent pair has a nontrivial intersection, Lemma 1
gives \(d_{\min}\le B\).  Therefore \(d_{\min}>B\) requires every disjoint
pair to have trivial intersection.  Independence and (1) give the first
bound in (3).  The bound \(p_k\le2/3\) gives the second.  Equations (4) and
the relative-distance statement follow. \(\square\)

## Lower rates

The same collision identifies the relevant lower-rate threshold.  If each
random block has dimension \(k=RB<B/2\), then

\[
  \mathbb E[|U_i\cap U_{i+1}|-1]
  =
  \frac{(2^k-1)^2}{2^B-1}
  =2^{-(1-2R)B+o(B)}.
  \tag{5}
\]

For \(B=c\log_2N\), the expected number of adjacent collisions changes scale
near

\[
  c(1-2R)=1.
\]

Equation (5) identifies the collision scale.  A complete lower-rate theorem
also needs a matching probability bound and control of non-collision profile
classes.  Those obligations remain open.

## Transpose plus region shuffles

The frozen Structured SPIN permutation is not the pure transpose analyzed
above.  It also samples coordinate assignments and independent permutations
inside the transposed regions.  Those permutations generally move a pair of
blocks to different relative positions in different columns.  The persistent
adjacency used in Lemma 1 then disappears.

That richer permutation is still computable, but total Hamming weight is not
a sufficient interface.  Condition on the vector of column weights.  Each
independently shuffled column has a two-state accumulator transfer matrix.
Multiplying the \(B\) matrices gives the conditional output enumerator.  The
remaining outer obligation is a joint envelope for the column-weight vector.

Thus the next exact computation has two inputs:

1. the mathematical distribution of every coordinate and region shuffle;
2. the joint column-profile law of the selected outer ensemble.

Without those two inputs, substituting the uniform-interleaver accumulator
bound would be invalid.

## Status

- **Proved:** adjacent equal block words give output weight at most \(B\).
- **Proved:** the exact half-rate intersection probability (1).
- **Proved:** pure transpose has vanishing relative distance for every
  \(B_N=o(N)\).
- **Not implied:** the same obstruction for transpose plus independent region
  shuffles.
- **Open:** the profile transfer and asymptotic constants for the richer
  factored permutation.
