# Translating the random baseline to Structured SPIN

## Conclusion

The random analysis identifies the correct structured-outer interface.  With
an independent uniform coordinate permutation inside every active outer
block, followed by transpose and independent region permutations, the outer
code enters the first moment through its ordinary weight spectrum.  The
structured inner still enters through its own state transfer.  Thus a
random-like outer spectrum can be separated cleanly from the remaining
RM2Sub proof obligation.

The decisive asymptotic quantity is the spectral excess per active block.
If that excess is \(o(B)\), it does not change the limiting sparse exponent.
If it is \(O(1)\), a logarithmic block schedule has the same leading constant
as the random-outer benchmark.  In contrast, a repeated fixed-constituent BA
family has only polynomial sparse suppression under the currently available
weight-only interface and therefore does not obtain \(B=O(\log N)\) from
that interface.

This document proves the spectrum-to-profile translation and its comparison
lemma.  It does not prove the RM2Sub transfer theorem or authenticate the
modeled frozen outer spectrum.

## Probability space

Let \(N=LB\).  A local outer map has input dimension \(K_B\) and output
length \(B\).  For each block \(i\in[L]\), setup independently samples a
local outer map \(O_i\) from the declared outer ensemble.  It also samples an
independent uniform coordinate permutation \(\sigma_i\gets S_B\).

After applying the \(L\) coordinate permutations, arrange the outer output
as an \(L\)-by-\(B\) matrix and transpose it.  For every region
\(j\in[B]\), independently sample \(\pi_j\gets S_L\) and permute the \(L\)
bits in that region.  All local outer randomness, coordinate permutations,
region permutations, and inner randomness are mutually independent unless a
different joint law is stated explicitly.

For a realized local outer code, let \(A_{O,h}\) be the number of nonzero
messages whose output has Hamming weight \(h\).  Define

\[
  \overline A_{B,h}:=\mathbb E_O[A_{O,h}].
  \tag{1}
\]

The expectation is omitted for a fixed deterministic constituent.  If the
same sampled outer realization is shared across several blocks, the product
formula below must retain the joint expectation; one may not replace it by a
power of (1).

## Symmetrized local profile

Let \(e_h(y_1,\ldots,y_B)\) denote the elementary symmetric polynomial of
degree \(h\).  Define the expected symmetrized nonzero profile

\[
  S_B(\boldsymbol y)
  :=
  \sum_{h=1}^{B}
  \frac{\overline A_{B,h}}{\binom Bh}
  e_h(y_1,\ldots,y_B).
  \tag{2}
\]

The denominator in (2) has a direct meaning.  Conditional on a local output
having weight \(h\), the uniform coordinate permutation makes its support
uniform among the \(\binom Bh\) subsets of size \(h\).

Fix a set of \(q\) active outer blocks.  For

\[
  \boldsymbol a=(a_1,\ldots,a_B)\in\{0,\ldots,q\}^{B},
\]

the coefficient

\[
  [\boldsymbol y^{\boldsymbol a}]S_B(\boldsymbol y)^q
  \tag{3}
\]

is the expected number of messages on those \(q\) blocks that produce
exactly \(a_j\) ones in transposed region \(j\).  Equation (3) follows by
expanding one copy of (2) for each independently sampled active block.  The
complete occupation generating polynomial is

\[
  (1+S_B(\boldsymbol y))^L.
  \tag{4}
\]

Therefore the per-block coordinate permutation is not cosmetic in the
structured analysis.  It removes support geometry inside each local
weight class.  It was redundant in the random-outer calculation only because
a nonzero output of a uniform random injection was already uniform over all
nonzero binary words.

## Exact accumulator specialization

Consider the region-shuffled transpose followed by one accumulator.  Let
\(K_{L,a}(z)\) be the two-state region kernel defined in
`REGION_SHUFFLED_TRANSPOSE_ACCUMULATOR.md`, and let
\(e_0=(1,0)^{\mathsf T}\) and \(\boldsymbol1=(1,1)^{\mathsf T}\).
Let \(Z_{D,q}\) count the nonzero messages supported on exactly \(q\) outer
blocks whose encoded output has weight less than \(D\).  For an integer
threshold \(D\) and \(z\in(0,1)\), its first moment satisfies

\[
\begin{split}
  \mathbb E[Z_{D,q}]
  \le{}&
  \binom Lq z^{-D}
  \sum_{\boldsymbol a\in\{0,\ldots,q\}^{B}}
  [\boldsymbol y^{\boldsymbol a}]S_B(\boldsymbol y)^q
  \\
  &\mathrel{}\cdot
  e_0^{\mathsf T}
  \left(\prod_{j=1}^{B}K_{L,a_j}(z)\right)
  \boldsymbol1.
  \tag{5}
\end{split}
\]

To prove (5), first condition on the region-count vector
\(\boldsymbol a\).  Each region permutation makes its input uniform on the
weight-\(a_j\) slice, so its conditional accumulator transform is
\(K_{L,a_j}(z)\).  Multiplying the kernels preserves the accumulator boundary
state between consecutive regions.  Finally use

\[
  \boldsymbol1\{w<D\}\le z^{-D}z^w.
\]

All coefficients and all entries in (5) are nonnegative.  This observation
is what permits the comparison theorem below.

For RM2Sub, equation (5) is an interface rather than a completed formula.
The two-state accumulator kernel must be replaced by the authenticated
RM2Sub state-and-boundary kernel.  If lane occupancy or another invariant is
preserved, that invariant must remain in the coefficient type.  The ordinary
weight spectrum suffices only after the coordinate and region actions have
been proved transitive on every type used by that kernel.

## Random-outer reference

For a uniform random injective map from \(\mathbb F_2^{K_B}\) to
\(\mathbb F_2^B\), a fixed nonzero message is uniform on
\(\mathbb F_2^B\setminus\{0\}\).  Its expected local spectrum is

\[
  A^{\mathrm{rnd}}_{B,h}
  =
  (2^{K_B}-1)\frac{\binom Bh}{2^B-1}.
  \tag{6}
\]

Hence its symmetrized profile is

\[
  S_B^{\mathrm{rnd}}(\boldsymbol y)
  =
  \frac{2^{K_B}-1}{2^B-1}
  \left(\prod_{j=1}^{B}(1+y_j)-1\right).
  \tag{7}
\]

Equation (7) is the exact nonzero-row law.  Relaxing the row to a uniform
binary word gives the binomial region mixture used in the earlier random
calculation, together with its explicit conditioning factor.

## Spectrum comparison lemma

Define the pointwise spectral excess, in bits per active block, by

\[
  g_B
  :=
  \max_{1\le h\le B:\,\overline A_{B,h}>0}
  \log_2\frac{\overline A_{B,h}}
  {A^{\mathrm{rnd}}_{B,h}}.
  \tag{8}
\]

Missing structured weights have multiplicity zero and impose no constraint.

**Lemma 1 (nonnegative-transfer comparison).**  Suppose the local outer
setups used by distinct blocks are independent.  For every \(q\) and every
nonnegative functional of the transposed region-count profile, the structured
expected contribution is at most \(2^{g_Bq}\) times the corresponding exact
random-outer contribution.

*Proof.*  Equation (8) gives the coefficientwise inequality

\[
  S_B(\boldsymbol y)
  \preceq 2^{g_B}S_B^{\mathrm{rnd}}(\boldsymbol y).
\]

Both polynomials have nonnegative coefficients.  Taking the \(q\)-th power
preserves the inequality.  Applying any nonnegative profile functional and
summing its coefficients also preserves it. \(\square\)

Lemma 1 is deliberately stronger than necessary.  A transfer-weighted
comparison can succeed even when the pointwise spectrum comparison fails.
The latter is the correct next interface if a few structured weights have a
large excess but negligible RM2Sub bad-output probability.

## Consequence for logarithmic blocks

Suppose the exact random-outer version of the same structured-inner transfer
satisfies the sparse bound

\[
  \mathbb E[Z^{\mathrm{rnd}}_{D,q}]
  \le
  \binom Lq N^a2^{-\lambda_{\mathrm{rnd}}Bq}.
  \tag{9}
\]

Lemma 1 gives

\[
  \mathbb E[Z^{\mathrm{str}}_{D,q}]
  \le
  \binom Lq N^a
  2^{-(\lambda_{\mathrm{rnd}}-g_B/B)Bq}.
  \tag{10}
\]

Thus, if

\[
  \overline\epsilon
  :=\limsup_{B\to\infty}\frac{g_B}{B}
  <\lambda_{\mathrm{rnd}},
  \tag{11}
\]

then the structured sparse exponent may be chosen below
\(\lambda_{\mathrm{rnd}}-\overline\epsilon\).  A schedule

\[
  B_N=\lceil c\log_2N\rceil_{\mathcal B}
\]

closes the sparse union bound whenever

\[
  c>\frac{a+1}
  {\lambda_{\mathrm{rnd}}-\overline\epsilon}.
  \tag{12}
\]

If \(g_B=o(B)\), the structured outer has the same limiting block constant as
the random outer for this transfer proof.  If \(g_B=O(1)\), the conclusion is
stronger: each active block costs only a constant factor.  If
\(g_B=\epsilon B+o(B)\), the inner sparse exponent loses \(\epsilon\).

Equation (12) does not import the random-convolution theorem into RM2Sub.  It
transfers an outer spectrum only when the random and structured outers are
paired with the same inner transfer.  The random-convolution result supplies
a benchmark and a target exponent; a separate RM2Sub kernel bound is still
required.

## Frozen outer diagnostic

The frozen Structured SPIN proof status points to the adjacent exploration
receipt
`parityfanout31x33_outer256_d38_expected_spectrum.json`.  Comparing that
modeled expected spectrum with (6), at \(B=256\) and \(K_B=128\), gives

\[
  g_{256}=0.0065094697569\text{ bits},
  \qquad
  \frac{g_{256}}{256}=2.5427616\times10^{-5}.
  \tag{13}
\]

The maximum occurs at weight 241.  On the interval from weight 37 to
weight 219, the maximum excess is only
\(2.83681\times10^{-6}\) bits; from weights 64 through 192, the
listed spectrum agrees with the random reference to binary64 precision.
The receipt's Bernoulli-envelope cost \(128.0065094697569\) is exactly the
random dimension cost 128 plus the excess in (13), up to the negligible
nonzero-row normalization.

These numbers explain why the frozen first-moment diagnostic is close to the
random-outer calculation.  They are not a proof.  The receipt assumes a
modeled even-floor source spectrum and nearest-binary64 arithmetic, and the
frozen `MAIN_CODE_FREEZE.md` explicitly requires a certified spectrum or
sufficient envelope for the actual outer code.

## Implications for candidate scalable outers

The random baseline divides candidate outer families into two cases.

1. A family with \(g_B=o(B)\), or with an equally strong transfer-weighted
   comparison, preserves the random sparse exponent and can support
   \(B=\Theta(\log N)\), subject to the RM2Sub sparse and bulk bounds.
2. A family whose one-active weighted mass is only
   \(B^{-p+o(1)}\) cannot close the global union bound with logarithmic
   blocks.  The sufficient schedule becomes
   \(B_N=N^{\beta+o(1)}\) with \(\beta>1/(p+1)\).

The proved sparse path for a repeated fixed constituent followed by a fixed
number of accumulators places the currently analyzed BA construction in the
second case under its weight-only proof interface.  This does not rule out a
different scalable BA family or a sharper RM2Sub-weighted argument.  It does
show that matching the random-like exponential spectrum is the substantive
outer-code objective; the transpose and region shuffles do not by themselves
create that exponent.

## Remaining proof obligations

- Authenticate the frozen ParityFanout-31x33 expected spectrum or a sufficient
  transfer-weighted envelope.
- Define a scalable outer ensemble and prove either \(g_B=o(B)\) or the
  weaker RM2Sub-weighted comparison uniformly in \(B\).
- Authenticate the RM2Sub state kernel and prove that its type is complete
  after coordinate, transpose, region, lane, and boundary actions.
- Prove a sparse bound uniform from \(q=1\) through growing \(q=o(L)\).
- Prove the linear-occupation exponent gap.
- Specify a bounded-gap admissible family before applying the arbitrary-length
  wrapper.

The next numerical calculation should therefore keep the actual RM2Sub
kernel fixed and replace only the outer spectrum: evaluate the one-active and
bounded-\(q\) transfer once with (6) and once with the candidate structured
spectrum.  Their ratio measures the structured penalty directly and avoids
conflating it with the inner-code penalty.
