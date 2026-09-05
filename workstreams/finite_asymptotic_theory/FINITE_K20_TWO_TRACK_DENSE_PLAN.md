# Two-track finite dense proof at \(k=2^{20}\)

## Decision

Keep the \(B=240\) construction and its linear-time online map. Pursue two
finite arguments in the probability space of
`FINITE_K20_INDEPENDENT_SETUP.md`.

1. At relative distance 11%, permit one exact coefficient argument that sums
   all BA shells directly.
2. At relative distance 10.9%, develop a five-band certificate as the primary
   fallback.

The first track is time-boxed. A relaxation that introduces a factor per row
or independently conditions every band count does not qualify as the exact
11% attempt.

## Exact all-active identity at 11%

Fix \(Q=L\), a Chernoff witness \(z\in(0,1)\), and one row of the conditioned
BA ensemble. After the row-coordinate permutation, let

\[
  \nu(S)
\]

be the expected counting mass of support \(S\subseteq[B]\). Each accepted BA
code has dimension 120. Therefore

\[
  \sum_{S\subseteq[B]}\nu(S)=2^{120}-1.
  \tag{1}
\]

Equation (1) is exact. The pointwise estimate
\(\nu(S)\le \overline A(|S|)/(g\binom B{|S|})\) need not sum to the right-hand
side of (1). Multiplying that pointwise estimate across all rows can therefore
discard material finite slack.

For supports \(S_1,\ldots,S_L\), define the weight of region \(b\) by

\[
  K_b:=|\{j:b\in S_j\}|.
\]

Let \(R_k(z)\) be the finite nonnegative RM2Sub transfer for one length-\(L\)
region whose input is a uniform word of weight \(k\). Independent region
permutations give

\[
\begin{split}
 \mathbb E[Z_{d,L}]
 \le z^{-d}
 \sum_{S_1,\ldots,S_L}
 \left(\prod_{j=1}^L\nu(S_j)\right)
 e_0^{\mathsf T}
 \left(\prod_{b=1}^B R_{K_b}(z)\right)\mathbf1.
 \tag{2}
\end{split}
\]

The products in (2) follow the serialized region order. No commutativity is
assumed.

Define the multivariate polynomial

\[
  F_\nu(x_1,\ldots,x_B)
  :=\sum_{S\subseteq[B]}\nu(S)\prod_{b\in S}x_b.
\]

Then the sum in (2) equals

\[
  \sum_{\boldsymbol k\in\{0,\ldots,L\}^B}
  [x_1^{k_1}\cdots x_B^{k_B}]F_\nu(x)^L\,
  e_0^{\mathsf T}
  \left(\prod_{b=1}^B R_{k_b}(z)\right)\mathbf1.
  \tag{3}
\]

Equations (1)--(3) are the exact coefficient target. They sum shell mixtures
without row-conditioning or band-composition factors. Computing or bounding
(3) sharply remains open.

## Rejected 11% shortcut

`diagnose_finite_k20_qL_column_holder.py` replaces every \(R_k(z)\) by a
common positive weighted norm. It then applies Hölder's inequality across the
240 regions. This operation sums shell mixtures, but it destroys their
cross-region correlation.

At the tested tilts, the resulting upper bound fails by more than 1.7 million
bits. The failure is too large for witness refinement or outward rounding to
repair. The column-Hölder route is closed as a practical 11% option. The exact
coefficient target (3) remains open, but it requires materially stronger
structure than a common norm.

## Five-band fallback at 10.9%

Set

\[
 d_{109}:=231045=\lfloor0.109N\rfloor
\]

and partition the permitted weights into

\[
 [23,54],\ [55,94],\ [95,145],\ [146,185],\ [186,217].
 \tag{4}
\]

For band \(b\), fix \(y_b\in(0,1)\) and a pointwise majorant

\[
  \nu_b\le c_b\mu_{y_b},
  \tag{5}
\]

where \(\mu_y\) is the Bernoulli-\(y\) product law on \(\{0,1\}^{240}\).
Fix a composition \(\boldsymbol n=(n_1,\ldots,n_5)\) with
\(\sum_b n_b=L\). Let \(\boldsymbol p\) be a positive categorical
probability vector and set

\[
  q:=\sum_{b=1}^5p_by_b.
\]

In each region, condition iid categorical labels with probabilities
\(\boldsymbol p\) to have counts \(\boldsymbol n\). After removing that
conditioning, every inner input bit is iid Bernoulli-\(q\). Thus

\[
\begin{split}
 T(\boldsymbol n)
 :={}&\binom{L}{n_1,\ldots,n_5}
       \prod_{b=1}^5c_b^{n_b}\\
 &\quad\cdot\min\left\{1,
   \frac{M_{\rm in}(q,z)z^{-d_{109}}}
   {\Pr_{\boldsymbol p}[\operatorname{Mult}(L,\boldsymbol p)
      =\boldsymbol n]^B}
   \right\}
 \tag{6}
\end{split}
\]

is a valid upper bound for that composition, subject to the RM2Sub transfer
lemma stated below. Summing (6) over all
\(\binom{L+4}{4}\) compositions covers \(Q=L\).

## Convex composition cells

Fix \((\boldsymbol p,\boldsymbol y,z)\). Suppose the second argument of the
minimum in (6) is below one throughout a convex composition cell. Apart from
constants, the logarithm of (6) is

\[
  (1-B)\log\binom{L}{n_1,\ldots,n_5}
  +\sum_{b=1}^5n_b(\log c_b-B\log p_b).
  \tag{7}
\]

The function \(-\log\binom{L}{\boldsymbol n}\) is convex on the continuous
simplex because \(\log\Gamma(x+1)\) is convex. Therefore (7) is convex. The
unclipped reference exponent is also convex. A fixed witness is valid on a
simplex when every vertex has negative reference exponent. The maximum of
(7) on that simplex then occurs at a vertex.

This observation reduces the finite composition proof to an adaptive
barycentric simplex cover. Each accepted simplex records one positive
categorical witness and one Chernoff witness.

## Bernoulli densities above one half

The three-state source note introduces a marked-position probability \(p\)
and sets the actual bit probability to \(p/2\). The five-band calculation can
produce \(q>1/2\). The transfer identities themselves extend to every
\(q\in(0,1)\).

For iid Bernoulli-\(q\) input and \(z\in(0,1)\), set

\[
 u:=1-q+qz,\qquad
 v:=q+(1-q)z,\qquad
 r:=1-q-qz.
\]

Fourier inversion uses integer powers of \(r\). Hence the formulas for
\(k_{00}\), \(k_{01}\), \(\overline d\), and \(m\) remain algebraically
valid when \(r<0\). Their final values are nonnegative transition masses.
Nearest binary64 arithmetic is unsafe because the signed sums can cancel.
The diagnostic therefore reevaluates every selected transfer with 100-digit
arithmetic.

The extension is exact in form but is not yet outward-certified. The final
verifier must evaluate the signed Fourier sums with rational arithmetic or
directed intervals.

## Current evidence

`diagnose_finite_k20_band_pure.py` proves no claims, but its 100-digit final
evaluations show that all five pure compositions close at \(Q=L\). Their
10.9% margins are approximately

\[
  192233,\ 58996,\ 6016,\ 59921,\ 196148\text{ bits}.
\]

`diagnose_finite_k20_band_compositions_qL.py` tests 116 pure, pairwise, and
Dirichlet compositions. Every sampled term has more than 6,016 bits of
margin. The pure central composition is the sampled maximum. The equal
five-way composition has about 88,971 bits of margin.

The first axis-aligned box cover was inefficient because its cells did not
follow the simplex faces. A barycentric longest-edge implementation now
certifies cells, but the interrupted depth-12 diagnostic did not produce a
complete receipt. Boundary cells require targeted refinement.

## Status and next action

- **Proved algebraically:** the coefficient identity (2)--(3), the multiband
  reduction (6), and the fixed-witness convexity statement (7), subject to the
  existing RM2Sub envelope assumptions.
- **Diagnostic:** all numerical margins in this note.
- **Rejected:** common-norm column Hölder as the 11% coefficient shortcut.
- **Open:** a sharp evaluator for (3), a complete rational simplex cover for
  (6), every occupation below \(L\), and outward arithmetic.

The next computation should finish the \(Q=L\) barycentric cover at 10.9%.
It should treat lower-dimensional faces separately, then recurse only on the
interior simplices that fail their first witness.
