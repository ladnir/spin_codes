# Goal 03: BCH moment gate for the two parity blocks

## Objective

The inner envelope depends only on the total number of active packets.  Goal
03 tests whether the ordinary BCH weight spectrum can control that packet
support after the two outer parity blocks are included.

The test compares two rigorous interfaces:

1. use the full spectrum on the data blocks and only deterministic activity
   information on the parity blocks;
2. use a second-moment inequality for two data blocks and both parity blocks.

The goal is diagnostic.  It evaluates low occupation shells before a complete
sum over all outer words.

## Local packet-support moment

For $x\in\mathbb F_{2^{64}}$, let $S_x$ be the number of nonzero packets
in the independently permuted BCH encoding of $x$.  For $0<t\le1$, define

\[
Q_x(t):=\mathbb E[t^{S_x}],
\]

where the expectation is over the permutation inside that BCH block.
The value $Q_x(t)$ depends only on the binary weight of the BCH encoding.
The exact BCH spectrum therefore determines

\[
M_1(t):=\sum_{x\ne0}Q_x(t),
\qquad
M_{2,*}(t):=\sum_{x\ne0}Q_x(t)^2,
\qquad
M_{2,\mathrm{all}}(t):=\sum_x Q_x(t)^2,
\qquad
M_3(t):=\sum_{x\ne0}Q_x(t)^3.
\]

## One active data block

Fix a data position with coefficient \(\alpha\ne0\).  Its three nonzero outer
symbols are $x,x,\alpha x$.  Hölder's inequality and the bijection
$x\mapsto\alpha x$ give

\[
\sum_{x\ne0}Q_x(t)^2Q_{\alpha x}(t)
\le M_3(t).
\]

This bound uses the ordinary spectrum.  It does not require the joint weight
profile of $x$ and $\alpha x$.

## At least two active data blocks

Fix an outer support of size $r\ge2$.  Choose two supported positions
$i\ne j$, and denote their values by $x,y$.  Fix the other $r-2$
nonzero values.  The two parity values are affine functions of $x,y$.
Their linear part is

\[
\begin{pmatrix}
1&1\\
\alpha_i&\alpha_j
\end{pmatrix}.
\]

The matrix is invertible because \(\alpha_i\ne\alpha_j\).  Hence the affine
map from $(x,y)$ to the parity pair is a bijection on \(\mathbb F^2\).

Apply Cauchy--Schwarz on the original domain \((\mathbb F^*)^2\).  The norm
of the two data factors is \(M_{2,*}(t)\).  The parity-pair image is a subset
of \(\mathbb F^2\), so its norm is at most \(M_{2,\mathrm{all}}(t)\).  Hence

\[
\sum_{x,y\ne0}
Q_x(t)Q_y(t)Q_{p_0}(t)Q_{p_1}(t)
\le M_{2,*}(t)M_{2,\mathrm{all}}(t).
\]

Summing the remaining data values gives the shell moment

\[
M_1(t)^{r-2}M_{2,*}(t)M_{2,\mathrm{all}}(t).
\]

## From moments to the inner envelope

For a packet cutoff $b$, Markov's inequality gives

\[
\mathbb E[\#\{\text{shell words with }H\le b\}]
\le
t^{-b}Z_r(t),
\]

where \(Z_1(t):=M_3(t)\) and
\(Z_r(t):=M_1(t)^{r-2}M_{2,*}(t)M_{2,\mathrm{all}}(t)\) for
\(r\ge2\).

Goal 02 supplies decreasing probability caps on consecutive support
intervals.  Summation by parts combines those caps with cumulative bounds at
the interval endpoints.  This avoids charging every cumulative count once
per interval.

## Acceptance criteria

Goal 03 is complete when the evaluator:

1. validates the exact BCH spectrum and every local packet-support law;
2. verifies the two parity coefficient matrix is invertible for distinct
   supported positions;
3. evaluates the Hölder, Cauchy, and parity-worst interfaces;
4. reports shell contributions for representative occupation sizes;
5. identifies whether a joint BCH-weight profile remains necessary.

The inequalities are exact.  The initial optimized numerical evaluation is
a floating-point diagnostic, not an outward-rounded certificate.
