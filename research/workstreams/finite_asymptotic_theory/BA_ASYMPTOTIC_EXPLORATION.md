# BA outer codes: sparse asymptotics and block-size schedules

## Question and conclusion

Finite BA spectra can look random-like at linear weights while retaining a
polynomial sparse boundary.  Structured SPIN does not require the outer block
size to be logarithmic.  It requires the total one-active-block contribution
to vanish.

For the repeated-constituent BA ensemble defined below, a fixed number of
accumulators gives a provable polynomial lower bound on the weighted spectrum.
This bound rules out an exponential weighted-spectrum hypothesis.  It does
not rule out a sublinear outer block.

The authenticated EBCH128 diagnostics strongly suggest that the polynomial
lower bound is tight.  A matching upper bound remains open.

## BA ensemble

Fix a binary linear constituent

\[
  C_0\subseteq\mathbb F_2^b
\]

with dimension \(k_0\), minimum distance \(d_0\), and weight enumerator

\[
  W_0(z)=\sum_{h=0}^b A_h^{(0)}z^h.
\]

For every multiple \(B=mb\), define the pre-accumulator code

\[
  C_B^{[0]}:=C_0^{\oplus m}\subseteq\mathbb F_2^B.
\]

For \(j\in\{1,\ldots,\ell\}\), setup samples an independent permutation
\(\pi_j\gets S_B\) and defines

\[
  C_B^{[j]}:=\operatorname{Acc}_B(\pi_j(C_B^{[j-1]})).
\]

The BA outer is \(C_B^{[\ell]}\).  The expectation below is over
\((\pi_1,\ldots,\pi_\ell)\).

For \(0<\rho<1\), define the weighted nonzero spectrum

\[
  M_{B,\ell}(\rho)
  :=\mathbb E\!\left[
    \sum_{c\in C_B^{[\ell]}\setminus\{0\}}
    \rho^{\operatorname{wt}(c)}
  \right].
\]

When the structured inner admits the pointwise transfer bound

\[
  \Pr[\text{bad output}\mid c]\le \rho^{\operatorname{wt}(c)},
\]

the quantity \(M_{B,\ell}(\rho)\) upper-bounds the contribution of all
messages supported inside one BA block.  A profile-dependent transfer theorem
may use a different weighted sum.

## Exact finite recursion

Let \(A_{B,h}^{[j]}\) be the expected number of words of weight \(h\) after
\(j\) accumulators.  The initial spectrum is exact:

\[
  \sum_{h=0}^B A_{B,h}^{[0]}z^h=W_0(z)^m.
\]

For input weight \(h\) and output weight \(w\), define

\[
  T_B(h,w)
  :=
  \binom{w-1}{\lceil h/2\rceil-1}
  \binom{B-w}{\lfloor h/2\rfloor}.
\]

A uniform interleaver followed by one accumulator gives

\[
  A_{B,w}^{[j+1]}
  =\sum_{h=1}^B
  A_{B,h}^{[j]}
  \frac{T_B(h,w)}{\binom Bh}.
  \tag{1}
\]

Equation (1) is the exact finite recursion used by the existing BA spectrum
analyzer.  It preserves the total expected number of codewords.

## Proved sparse-path lower bound

Let

\[
  h_0:=d_0,
  \qquad
  h_j:=\left\lceil\frac{h_{j-1}}2\right\rceil
  \quad(1\le j\le\ell),
\]

and define

\[
  p_\ell(d_0):=\sum_{j=1}^{\ell}h_j-1.
  \tag{2}
\]

**Lemma 1 (constant-weight boundary path).** Fix \(C_0\), \(\ell\), and
\(\rho\in(0,1)\).  Along multiples \(B=mb\),

\[
  M_{B,\ell}(\rho)=\Omega\!\left(B^{-p_\ell(d_0)}\right).
  \tag{3}
\]

*Proof.* There are \(mA_{d_0}^{(0)}=\Theta(B)\) pre-accumulator words that
activate one constituent and have weight \(d_0\).  Restrict each accumulator
transition to output weight \(h_j\).  For fixed \(h\) and
\(r=\lceil h/2\rceil\),

\[
  \frac{T_B(h,r)}{\binom Bh}=\Theta(B^{-r}).
\]

The product of the \(\ell\) transition probabilities is
\(\Theta(B^{-\sum_jh_j})\).  Multiplying by the \(\Theta(B)\) initial words
and by the constant \(\rho^{h_\ell}\) proves (3). \(\square\)

Lemma 1 is an obstruction to one proof interface, not an impossibility theorem
for Structured SPIN.  It shows that a fixed-\(\rho\) weighted-spectrum bound
cannot decay faster than a power of \(B\).  A finer structured transfer law
could assign additional \(N\)-dependent suppression to these boundary types.

## Consequences for the outer block size

Let the complete Structured SPIN length be \(N\).  There are
\(L_N=N/B_N\) outer blocks.  Under the weight-only interface, the one-active
contribution is at most

\[
  \frac{N}{B_N}M_{B_N,\ell}(\rho).
  \tag{4}
\]

If one proves the matching upper bound

\[
  M_{B,\ell}(\rho)
  \le B^{-p+o(1)}+2^{-\Omega(B)},
  \tag{5}
\]

then every schedule

\[
  B_N=N^{\beta+o(1)},
  \qquad
  \beta>\frac1{p+1},
  \tag{6}
\]

makes (4) vanish.  Lemma 1 shows that this threshold cannot be improved
through the same weighted-spectrum interface when \(p=p_\ell(d_0)\).

More generally, the decay of \(M_B\) determines the schedule:

| Weighted-spectrum bound | Sufficient block schedule |
| --- | --- |
| \(M_B\le2^{-\kappa B}\) | \(B_N=\Theta(\log N)\) |
| \(M_B\le2^{-\kappa B^\alpha}\) | \(B_N=\Theta((\log N)^{1/\alpha})\) |
| \(M_B\le B^{-p+o(1)}\) | \(B_N=N^{\beta+o(1)}\), \(\beta>1/(p+1)\) |

Every displayed schedule is sublinear.  The arbitrary-length rounding loss is
therefore \(o(N)\).

## EBCH128 boundary prediction

For the fixed EBCH \([128,64,22]\) constituent, \(d_0=22\).  The predicted
exponents are:

| Accumulators | Boundary path | \(p_\ell(22)\) | Candidate schedule after a matching upper bound |
| ---: | --- | ---: | --- |
| 1 | \(22\to11\) | 10 | \(B_N=N^{1/11+\varepsilon}\) |
| 2 | \(22\to11\to6\) | 16 | \(B_N=N^{1/17+\varepsilon}\) |
| 3 | \(22\to11\to6\to3\) | 19 | \(B_N=N^{1/20+\varepsilon}\) |

The two-accumulator schedule is the current BA-3 case.  An additional
accumulator improves the sparse exponent by only three because the boundary
weight has already fallen to six.

## Finite diagnostics

The script `analyze_ba_asymptotics.py` composes the authenticated
`scripts/ebch128_64_spectrum.csv` with the existing binary64 implementation
of (1).  The calculations are diagnostics; they are not outward-rounded
certificates.

At \(\rho=0.1\) and two accumulators, the results are:

| \(B\) | \(\log_2 M_{B,2}(0.1)\) | Boundary lower bound | Gap (bits) | First weight with expected multiplicity at least one |
| ---: | ---: | ---: | ---: | ---: |
| 256 | -66.7171 | -73.6143 | 6.8972 | 30 |
| 512 | -83.7900 | -89.8153 | 6.0253 | 58 |
| 1024 | -100.2088 | -105.9148 | 5.7060 | 114 |
| 2048 | -116.3987 | -121.9644 | 5.5656 | 227 |
| 4096 | -132.4896 | -137.9891 | 5.4994 | 453 |

The fitted log-log slope over \(B\in\{1024,2048,4096\}\) is approximately
\(-16.14\).  The boundary lower bound has fitted slope approximately
\(-16.04\).  The gap decreases toward about 5.5 bits.

The same check at \(B\in\{512,1024,2048\}\) gives:

| Accumulators | Predicted slope | Fitted full-spectrum slope | Fitted boundary slope |
| ---: | ---: | ---: | ---: |
| 1 | -10 | -10.21 | -10.06 |
| 2 | -16 | -16.30 | -16.07 |
| 3 | -19 | -19.25 | -19.08 |

For two accumulators, changing \(\rho\) among 0.05, 0.1, and 0.2 leaves the
observed power close to \(-16\).  This stability is expected because \(\rho\)
changes the constant price of the final bounded weight, not the power of
\(B\).

These diagnostics support the conjecture

\[
  M_{B,2}(\rho)=\Theta_\rho(B^{-16})
  \tag{7}
\]

for the repeated EBCH128 ensemble and fixed \(\rho\) in the tested range.
They do not prove the upper bound in (7).

## Missing provenance for other finite diagnostics

The existing EBCH32 and shortened-XBCH64 BA receipts refer to

- `scripts/ebch32_16_delta8_spectrum.csv`; and
- `scripts/xbch64_32_philips_spectrum.csv`.

Those files are absent from this worktree.  This exploration did not
reconstruct them from derived receipts.  Their predicted sparse exponents are
still determined by their stated minimum distances:

\[
  p_2(8)=5,
  \qquad
  p_2(12)=8.
\]

Numerical confirmation requires the original authenticated spectra.

## Matching upper-bound target

The next theorem should prove, for fixed \(\rho\) in the transfer range,

\[
  M_{B,2}(\rho)\le C_\rho B^{-16}+2^{-\kappa_\rho B}
  \tag{8}
\]

for all multiples of 128.  A useful proof partition is:

1. one active EBCH128 constituent;
2. a fixed number \(q\ge2\) of active constituents;
3. \(q\to\infty\) with \(q=o(B)\); and
4. \(q=\Theta(B)\).

The first class should yield the \(B^{-16}\) leading term.  Each fixed
\(q\ge2\) should lose a strictly larger power.  The final class should follow
from a negative linear-weight saddle-point exponent.  The third class is the
main uniformity gap between the fixed-\(q\) and linear-weight analyses.

## Complexity implication

Repeating a fixed constituent, applying a fixed number of permutations, and
applying a fixed number of accumulators costs \(O(B)\) per BA block.  Across
\(N/B\) blocks, the outer layer costs \(O(N)\) for every sublinear schedule in
(6).  The factored transpose and region permutations also route \(O(N)\)
symbols.

The block-size exponent affects setup structure and proof constants, but it
does not by itself destroy linear encoding time.

## Status

- **Proved:** the exact recursion (1) and sparse-path lower bound (3).
- **Diagnostic:** the finite EBCH128 slopes and constant gaps.
- **Conjectured:** the matching upper bound (8).
- **Conditional consequence:** (8) implies the admissible schedule
  \(B_N=N^{1/17+\varepsilon}\).
- **Open:** the type-dependent replacement for \(M_B(\rho)\) under the actual
  Structured SPIN interleaver and RM2Sub transfer law.
