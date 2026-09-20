# Goal 01: exact inner transfer bound

## Objective

For **Riffle RandomStepConv g=4 sigma=20**, derive and implement a rigorous
upper bound on

\[
p_h(D):=
\Pr\!\left[\operatorname{wt}(Y)\le D
\;\middle|\;
\operatorname{supp}(X)\text{ is a uniform }h\text{-subset of }[N]
\right].
\]

The probability is over the global packet permutation and the independent
setup matrices (M_1,\ldots,M_N).  Fix the nonzero values of the (h) active
input packets arbitrarily.  The random-map symmetry makes (p_h(D))
independent of those values.

This goal concerns only the inner encoder.  It does not use an outer spectrum
or claim an end-to-end distance.

## Exact transfer enumerator

Define

\[
b(z):=2^{-g}(1+z)^g,
\qquad
q:=2^{-\sigma}.
\]

Use state classes (Z:=\{0^\sigma\}) and
\(L:=\mathbb F_2^\sigma\setminus Z\).  With rows indexed by the current
state class and columns by the next class, define

\[
T_0(z):=
\begin{pmatrix}
1&0\\
q b(z)&(1-q)b(z)
\end{pmatrix},
\qquad
T_1(z):=
\begin{pmatrix}
q b(z)&(1-q)b(z)\\
q b(z)&(1-q)b(z)
\end{pmatrix}.
\]

The subscript records whether the current input packet is zero or nonzero.
Let (x) mark a nonzero input packet and set

\[
T(x,z):=T_0(z)+xT_1(z).
\]

Because the baseline discards the terminal state, define

\[
F_N(x,z):=(1,0)T(x,z)^N(1,1)^\mathsf T.
\]

The coefficient

\[
[x^h z^d]F_N(x,z)
\]

is the sum, over every (h)-subset of packet positions, of the probability
that the inner output has weight (d).  This coefficient averages over every
compatible sequence of random linear maps.

## Rank-one gap reduction

The active transition has rank one.  Write

\[
T_1(z)=\begin{pmatrix}1\\1\end{pmatrix}
\begin{pmatrix}qb(z)&(1-q)b(z)\end{pmatrix}.
\]

Fix (h\ge1\) active positions.  Let \(\ell_j\) be the number of inactive
positions after active position (j) and before the next active position.
For (j=h), let \(\ell_h\) count the inactive suffix.  The inactive prefix
does not affect the state or output.

Define

\[
d(z):=(1-q)b(z),
\qquad
G_z(r):=\frac{b(z)-d(z)r}{(1-r)(1-d(z)r)}.
\]

Define

\[
\alpha(z):=\frac{qb(z)}{1-d(z)}.
\]

For a gap of length \(\ell\), the exact scalar moment is

\[
f_\ell(z)=\alpha(z)+(1-\alpha(z))d(z)^{\ell+1}.
\]

The scalar (G_z(r)=\sum_{\ell\ge0}f_\ell(z)r^\ell) is the generating
function for one post-activation gap.  The exact fixed-support moment is

\[
[x^h]F_N(x,z)
=
[r^{N-h}]\frac{G_z(r)^h}{1-r}.
\]

This identity removes the zero-input coefficient that makes a direct
bivariate Chernoff bound ineffective for sparse support.

## Proof-facing coefficient bound

Expanding each factor (f_\ell(z)) by its two summands records how many live
episodes terminate.  For (0\le e\le h\), set (k:=h-e).  Then

\[
[x^h]F_N(x,z)
=
\sum_{e=0}^h
\binom he
\alpha(z)^e
\big((1-\alpha(z))d(z)\big)^k
[r^{N-h}]
\frac{1}{(1-r)^{e+1}(1-d(z)r)^k}.
\]

This positive decomposition permits a separate coefficient tilt for every
termination count.  For arbitrary radii \(0<r_e<1\), nonnegativity gives

\[
p_h(D)
\le
\frac{z^{-D}}{\binom Nh}
\sum_{e=0}^h
\binom he
\frac{\alpha(z)^e\big((1-\alpha(z))d(z)\big)^{h-e}}
{r_e^{N-h}(1-r_e)^{e+1}(1-d(z)r_e)^{h-e}}.
\]

The implementation minimizes every component analytically over (r_e), then
minimizes the sum over (z).  The target interface is the infimum of the
displayed bound over (z) and the component radii.

The outer BCH construction activates at least 18 packets.  Trivial bounds for
smaller (h) therefore do not affect the eventual outer integration.

## Acceptance criteria

Goal 01 is complete when all of the following conditions hold:

1. A proof derives both transfer matrices from the random-map experiment.
2. An exhaustive checker enumerates maps for small (g,\sigma,N) and matches
   the transfer coefficients exactly.
3. A target-size evaluator computes the tilted bound without coefficient
   expansion.
4. The evaluator reports the optimizing tilts and numerical safety margins.
5. A receipt distinguishes exact identities, outward bounds, and ordinary
   floating-point diagnostics.

The target-size run will initially scan representative values of (h).  It
will not yet sum the outer code spectrum.
