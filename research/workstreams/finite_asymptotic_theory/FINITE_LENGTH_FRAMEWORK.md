# Finite-length first-moment framework

## Purpose and status

This document states theorem targets for a SPIN encoder at one finite length.
The first-moment identities below are proved by elementary counting and
linearity.  A construction-specific distance claim still requires certified
outer-spectrum and inner-transfer bounds.

The framework distinguishes three assertions:

- an ensemble has a small expected number of bad nonzero codewords;
- a sampled code has the requested distance except with a bounded probability;
- at least one deterministic code with the requested parameters exists.

The first assertion implies the other two only through the explicit argument
given below.  A construction exists even when the displayed bound exceeds one;
that case simply gives no distance guarantee.

## Common notation

Fix integers \(K,N\ge 1\) and a target integer distance \(D\in\{1,\ldots,N\}\).
The intermediate and output alphabets are \(\mathbb F_2^N\).

A setup experiment samples

\[
  (O,\Pi,I)\gets \mathsf{Setup}_{K,N}.
\]

Here

- \(O:\mathbb F_2^K\to\mathbb F_2^N\) is the outer linear map;
- \(\Pi:\mathbb F_2^N\to\mathbb F_2^N\) is a coordinate permutation; and
- \(I:\mathbb F_2^N\to\mathbb F_2^N\) is the inner linear map.

The resulting encoder is

\[
  E:=I\circ\Pi\circ O.
\]

All probabilities and expectations in a theorem refer to the stated setup
experiment.  Separate sampling operations use fresh independent randomness
unless the experiment declares a dependence.

Define

\[
  Z_D:=\left|\left\{x\in\mathbb F_2^K\setminus\{0\}:
  \operatorname{wt}(E(x))<D\right\}\right|.
\]

The strict inequality is deliberate.  The event \(Z_D=0\) says both that
\(E\) is injective and that its image is a binary \([N,K,d_{\min}\ge D]\)
linear code.

## Exact theorem for a deterministic outer code

First fix an injective outer map \(O\).  Let

\[
  A_h^{\mathrm{out}}
  :=\left|\left\{x\in\mathbb F_2^K\setminus\{0\}:
  \operatorname{wt}(O(x))=h\right\}\right|.
\]

Assume that \(\Pi\) is uniform in the symmetric group \(S_N\).  Assume also
that \(I\) is independent of \(\Pi\).  The inner map may be deterministic.
For \(h\in\{1,\ldots,N\}\), sample \(U_h\) uniformly from the Hamming slice

\[
  \mathcal S_{N,h}:=\{u\in\mathbb F_2^N:\operatorname{wt}(u)=h\}
\]

independently of \(I\), and define

\[
  p_h^{\mathrm{in}}(D)
  :=\Pr_{I,U_h}\!\left[\operatorname{wt}(I(U_h))<D\right].
\]

**Theorem 1 (exact finite first moment).** Under the preceding experiment,

\[
  \mathbb E[Z_D]
  =\sum_{h=1}^{N} A_h^{\mathrm{out}}p_h^{\mathrm{in}}(D),
\]

and

\[
  \Pr\!\left[E\text{ is not injective or }d_{\min}(E)<D\right]
  \le \sum_{h=1}^{N} A_h^{\mathrm{out}}p_h^{\mathrm{in}}(D).
\]

*Proof.* For each nonzero \(x\), the word \(\Pi(O(x))\) is uniform in
\(\mathcal S_{N,h}\), where \(h=\operatorname{wt}(O(x))\).  Linearity of
expectation and grouping by \(h\) give the identity.  The bad event is
\(\{Z_D\ge 1\}\), so Markov's inequality gives the probability bound. \(\square\)

The theorem is exact before any spectrum or transfer inequality is inserted.
It returns a number for every finite parameter choice.  It does not promise
that the number is below one.

## Joint theorem for a random outer code

The outer code may be sampled.  Independence must then be stated before the
expected outer spectrum is multiplied by an averaged inner probability.

Sample \(O\gets\mathcal O_{K,N}\), where the ensemble is supported on
injective outer maps.  Independently sample a uniform
\(\Pi\gets S_N\) and \(I\gets\mathcal I_N\).  Define

\[
  A_h^{\mathrm{out}}
  :=\mathbb E_O\!\left[
  \left|\left\{x\ne 0:\operatorname{wt}(O(x))=h\right\}\right|
  \right].
\]

Use the same definition of \(p_h^{\mathrm{in}}(D)\) as above.

**Corollary 2 (independent random outer).** In this joint experiment,

\[
  \mathbb E[Z_D]
  =\sum_{h=1}^{N} A_h^{\mathrm{out}}p_h^{\mathrm{in}}(D),
\]

and the same sum upper-bounds the probability of noninjectivity or distance
below \(D\).

This factorization can fail if \(O\) and \(I\) share setup randomness.  In that
case one must condition on the shared setup or use the type theorem below.

## Envelope theorem with an explicit outer failure event

Some outer analyses prove a useful spectrum bound only on a good event.
Let \(G_{\mathrm{out}}\) be an event determined by \(O\), and suppose

\[
  \Pr[G_{\mathrm{out}}^c]\le\varepsilon_{\mathrm{out}}.
\]

For every \(h\), assume

\[
  \mathbb E\!\left[
    \left|\{x\ne0:\operatorname{wt}(O(x))=h\}\right|
    \mathbf 1_{G_{\mathrm{out}}}
  \right]\le \overline A_h
\]

and \(p_h^{\mathrm{in}}(D)\le\overline p_h\).  Independence of \(O\) from
\((\Pi,I)\) remains in force.

**Corollary 3 (finite envelope).** The bad-setup probability is at most

\[
  \varepsilon_{\mathrm{out}}
  +\sum_{h=1}^{N}\overline A_h\overline p_h.
\]

The indicator in the outer expectation avoids an invalid conditional-to-
unconditional step.  An alternative is to state a pointwise spectrum envelope
for every outer realization in \(G_{\mathrm{out}}\).

## Type-orbit theorem for structured interleavers

A factored Structured SPIN interleaver need not make a fixed weight-\(h\) word
uniform on the complete Hamming slice.  Its orbit can depend on the word's
block occupation, region counts, or another invariant.  The following theorem
is the exact replacement.

Let \(\mathcal T\) be a finite type set.  For every outer realization \(O\),
let

\[
  \tau_O:\mathbb F_2^K\setminus\{0\}\to\mathcal T
\]

assign a type to each nonzero message.  Define the conditional outer count

\[
  A_t(O):=|\{x\ne0:\tau_O(x)=t\}|.
\]

The setup may sample \(O\), \(\Pi\), and \(I\) jointly.  For every
\((O,t)\) with \(A_t(O)>0\), define

\[
  p_{O,t}(D)
  :=\frac{1}{A_t(O)}
  \sum_{\substack{x\ne0\\\tau_O(x)=t}}
  \Pr\!\left[\operatorname{wt}(I(\Pi(O(x))))<D\mid O\right].
\]

The conditional probability includes all remaining setup randomness.  This
definition does not assume transitivity inside a type.  If transitivity is
proved, the average can be replaced by the probability for one representative.

**Theorem 4 (exact type first moment).** For any joint setup distribution,

\[
  \mathbb E[Z_D]
  =\mathbb E_O\!\left[\sum_{t\in\mathcal T}A_t(O)p_{O,t}(D)\right].
\]

Consequently, the right-hand side upper-bounds the probability that \(E\) is
noninjective or has minimum distance below \(D\).

*Proof.* Expand \(Z_D\) as a sum of indicators over nonzero messages.  Take
conditional expectation given \(O\), group the indicators by \(t\), and then
average over \(O\).  Apply Markov's inequality. \(\square\)

For independent outer and inner randomness, a useful sufficient factorization
is

\[
  \mathbb E[Z_D]
  \le\sum_{t\in\mathcal T}\overline A_t\overline p_t,
\]

provided either

\[
  \mathbb E[A_t(O)]\le\overline A_t
  \quad\text{and}\quad
  p_{O,t}(D)\le\overline p_t\ \text{for every relevant }O,
\]

or an independently justified conditional product bound is available.  Two
separate averaged bounds do not, by themselves, bound the average product.

## Existence and margin statements

Let

\[
  \mu_D:=\mathbb E[Z_D].
\]

The exact conclusions are:

1. \(\Pr[Z_D\ge1]\le\mu_D\).
2. If \(\mu_D<1\), then at least one setup gives a binary
   \([N,K,d_{\min}\ge D]\) code.
3. If \(\mu_D\le2^{-\lambda}\), then the setup failure probability is at
   most \(2^{-\lambda}\).  The quantity \(\lambda\) is the first-moment
   margin in bits.
4. If \(\mu_D\ge1\), the construction still exists for every setup.  The
   first-moment calculation gives no distance certificate at \(D\).

A finite theorem should therefore report \((K,N,D,\mu_D)\), the exact setup
distribution, and every mathematical input used to upper-bound \(\mu_D\).
It should not promise one distance or margin uniformly over all small lengths.

## Status of the source materials

The following facts are already proved in the cited manuscript material:

- A uniform permutation sends a fixed weight-\(h\) word uniformly to the
  weight-\(h\) Hamming slice.
- The accumulator input-output enumerator is exact:
  \[
    A^{\mathrm{Acc}}_{h,j}
    =\binom{j-1}{\lceil h/2\rceil-1}
     \binom{N-j}{h-\lceil h/2\rceil}
  \]
  on its feasible range.
- The accumulator tail obeys
  \(p_h(D)\le(4eD/N)^{\lceil h/2\rceil}\) when
  \(D/N<1/4\) and \(h\le N/2\).
- The random dense inner material supplies explicit finite tail envelopes,
  subject to its stated setup distribution and parameter restrictions.

The following items are not yet unconditional theorems:

- The sampled-grid check behind the Random SPIN \(0.109\) checkpoint is not
  an interval proof over the full parameter continuum.
- The frozen Structured SPIN first-moment numbers use nearest binary64
  arithmetic.
- The frozen Structured SPIN outer spectrum is a modeled spectrum, not a
  proved spectrum of the actual structured outer constituent.
- A factored structured permutation must use Theorem 4 unless its required
  slice-uniformity or orbit-transitivity statement is proved.

## Certificate interface for one requested length

A finite certificate generator should accept a construction manifest and a
target \(D\).  It should emit:

- \(K,N,D\) and the admissibility checks;
- a complete description of the setup probability space;
- exact or outward-rounded outer bounds;
- exact or outward-rounded inner-transfer bounds;
- the complete summation partition;
- an outward upper bound on \(\mu_D\);
- source and input hashes; and
- a verifier result that recomputes the bound independently.

Binary64 discovery output may select partitions or witnesses.  It is evidence
for a candidate certificate, not evidence for the theorem itself.
