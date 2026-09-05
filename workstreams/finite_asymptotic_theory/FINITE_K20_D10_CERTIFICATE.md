# Finite 10% certificate for EBCH32--ParityFanout--BA

## Certified statement

Set

\[
 k=2^{20},\qquad N=2^{21},\qquad d=209715=\lfloor0.10N\rfloor.
\]

The bounded setup algorithm makes at most six independent row draws for each
of 8,192 rows. It tests a candidate row by enumerating all \(2^{128}\)
messages. It accepts the candidate exactly when every nonzero row word has
weight in \([24,232]\), and it selects the first accepted candidate. Setup
returns failure if any row has no accepted candidate after six draws.

Conditioned on setup success, this procedure samples the independent-row
conditional ensemble specified in
`FINITE_K20_EBCH32_PARITYFANOUT_SETUP.md`. For the resulting linear encoder

\[
 E:\mathbb F_2^{2^{20}}\longrightarrow\mathbb F_2^{2^{21}},
\]

define

\[
 Z_d(E):=
 \left|\left\{x\ne0:\operatorname{wt}(E(x))\le d\right\}\right|.
\]

The outward receipts prove

\[
 \mathbb E[Z_d(E)]<2^{-51}.
 \tag{1}
\]

Markov's inequality therefore gives, conditioned on setup success,

\[
 \Pr[d_{\min}(E)\le209715]<2^{-51}<2^{-40}.
 \tag{2}
\]

The setup receipt also proves that one unconditioned row draw is rejected with
probability at most

\[
 q:=\mathtt{0x1.362062013db23p-10}.
\]

A union bound places the probability that any row exhausts six draws at most
\(8192q^6\). Exact rational addition of this abort bound and the conditional
distance bound proves

\[
 \Pr[\text{setup fails or }d_{\min}(E)\le209715]<2^{-45}<2^{-40}.
 \tag{3}
\]

If setup succeeds and the event in (3) does not occur, then

\[
 d_{\min}(E)\ge209716,
 \qquad
 \frac{d_{\min}(E)}N\ge0.10000038146972656>0.10.
\]

The probability in (3) is over the complete code-setup experiment. Message
words are enumerated; they are not sampled.

## Outer setup

Each of the 8,192 outer rows starts with eight copies of the genuine
\([32,16,8]\) extended BCH code.  Its exact constituent enumerator is

\[
 1+620z^8+13888z^{12}+36518z^{16}
 +13888z^{20}+620z^{24}+z^{32}.
\]

The row then applies an independently sampled ParityFanout-31x33 map and two
independently interleaved terminated accumulators.  The setup conditions each
row on having no nonzero final word outside weights \([24,232]\).  The outward
setup receipt proves

\[
 \Pr[G_{256}]\ge0.998816961305865
\]

for one unconditioned row draw. Unbounded rejection sampling therefore needs
at most 1.0011844399324061 trials per accepted row in expectation. The
certified bounded algorithm uses six attempts. Conditional on finding an
accepted draw, first-success selection has exactly the same conditional row
law as unbounded rejection sampling. Independence across rows is preserved
when the complete setup succeeds.

All rows are independent.  The route permutations, region permutations, and
RM2Sub-S19 multipliers are independent of the rows and of one another, as
specified by the setup document.

## First-moment partition

Let \(Q\) be the number of active outer rows.  The receipts partition every
nonzero message as follows.

| Occupations | Bad-weight threshold | Certified expected-count margin |
|---:|---:|---:|
| \(Q=1\) | 230686 | 51.845632 bits |
| \(2\le Q\le64\) | 230686 | 95.572310 bits |
| \(65\le Q\le8192\) | 209715 | 468.225918 bits |

The first two receipts bound a larger bad-weight event.  They therefore remain
valid at threshold 209715.  The combined verifier adds their exact binary64
upper bounds to the conservative dense bound \(2^{-468}\).  The exact rational
sum is less than \(2^{-51}\), which proves (1).

## Dense one-band reduction

The 11% proof needs several outer-weight bands because the BA boundary shells
are too expensive under one Bernoulli reference.  At 10%, one band suffices.

For an accepted outer row, let \(A_h\) be the conditioned expected number of
words of weight \(h\).  For \(r>1\) and \(y\in(0,1)\), define

\[
 M_r(y):=
 \sum_{h=24}^{232}
 A_h^r
 \left[\binom{256}{h}y^h(1-y)^{256-h}\right]^{1-r}.
 \tag{4}
\]

Fix an occupation \(Q\), an active-label probability \(p\in(0,1)\), a value
probability \(y\), a Chernoff surprisal \(s>0\), and an order \(r>1\).  Let
\(I(py,s)\) be the logarithm of the finite RM2Sub moment for Bernoulli-\(py\)
input.  Define

\[
 R_Q:=I(py,s)+ds
 -256\log\!\left[
   \binom{8192}{Q}p^Q(1-p)^{8192-Q}
 \right].
 \tag{5}
\]

When \(R_Q<0\), the fixed-witness Hölder bound is

\[
 \log\mathbb E[Z_{d,Q}]
 \le
 \log\binom{8192}{Q}
 +\frac Qr\log M_r(y)
 +\frac{r-1}{r}R_Q.
 \tag{6}
\]

For fixed \((p,y,s,r)\), both \(R_Q\) and the right side of (5) are convex in
real \(Q\in[0,8192]\).  Indeed, \(\log\binom{8192}{Q}\), extended with the
gamma function, is concave.  Its coefficient in (5) is

\[
 1-256\frac{r-1}{r}<0
\]

for every saved witness.  All other terms are affine in \(Q\).  Therefore,
checking both endpoints of an interval proves (5)--(6) throughout that
interval.

The diagnostic selects 19 fixed witnesses whose integer intervals partition
\([65,8192]\).  The Arb verifier recomputes every spectrum moment, inner
transfer, endpoint constraint, endpoint exponent, and interval contribution.
It then outward-sums the 19 contributions.

## Arithmetic and artifacts

`certify_ebch32_parityfanout_ba_one_band_interval_cover_d10_outward.py` uses
256-bit Arb balls.  Binary64 witness values are interpreted as exact rational
numbers.  Optimization success flags and diagnostic objective values are not
trusted.

The final receipt is
`ebch32_parityfanout31x33_ba3_B256_combined_outward_d10.json`.  The certificate
manifest binds the proof note, generators, verifiers, receipts, outer spectrum,
and frozen RM2Sub selection by SHA-256 hash.

`certify_ebch32_parityfanout_ba_bounded_setup_d10.py` performs the exact
rational setup-abort calculation. Its receipt is
`ebch32_parityfanout31x33_ba3_B256_bounded_setup_combined_outward_d10.json`.
The combined failure margin is 45.323934 bits for display; the proof checks the
integer claim \(<2^{-45}\).

## Scope

This certificate proves a finite setup-and-distance statement. The setup
algorithm is exact and terminating, but its exhaustive acceptance test costs
at most

\[
 6\cdot8192\cdot2^{128}
\]

row-message evaluations. It is therefore a mathematical closure, not a
practical setup implementation. Two production obligations remain.

1. Setup needs an efficient exact or one-sided-safe test for \(G_{256}\), or
   an authenticated precomputed set of accepted rows.
2. The optimized implementation must be audited against the declared
   ParityFanout, accumulator, route, and RM2Sub probability space.

The online encoder has linear work because every constituent and inner step
has fixed size and every permutation makes one pass over its coordinates. The
theorem does not claim linear-time setup. The frozen implementation baseline
is the construction directory
`riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19`;
its source does not by itself close the implementation-equivalence obligation
for this EBCH32--BA row law.
