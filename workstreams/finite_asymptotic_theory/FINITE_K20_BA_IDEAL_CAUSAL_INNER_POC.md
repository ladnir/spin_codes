# Finite \(k=2^{20}\) BA plus ideal causal-inner experiment

## Status and conclusion

This document records the initial nearest-binary64 diagnostic. The later
prefix-rank route supersedes its incomplete occupation-one conclusion for a
fixed setup. `FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md` now gives a
complete outward all-message certificate for one frozen repeated-BA setup.

The experiment replaces RM2Sub-S19 by invertible random lower-triangular
Toeplitz convolution with packet width \(g=1\). Sample \(h_1,\ldots,h_{N-1}\)
independently and uniformly, set \(h_0=1\), and define

\[
 y_t=\sum_{i=0}^t h_i u_{t-i}.
\]

For each fixed nonzero input difference, the inner emits zero before its first
nonzero input bit. The activation output is one, and all later output bits are
independent and fair. The triangular suffix map from the fresh \(h_i\)
coefficients proves this law exactly.

For occupation \(Q=1\), the reuse-aware first moment has only

\[
 24.305345\text{ bits of margin}.
\]

It therefore does not reach the 40-bit target. Outer weights 2 through 5
dominate. Deleting every outer shell below weight 13 raises this same
occupation-one calculation to 40.58 bits. Thus, for the unconditioned
Golay--BA-3 ensemble, strengthening only the causal inner does not resolve the
finite low-occupation first-moment bound. A low-shell outer property is still
needed by this proof route.

This conclusion concerns the available first-moment proof. It does not prove
that a typical sampled code has a low-distance word, and it does not compare
the ideal inner pointwise with RM2Sub-S19 for every input.

## Construction used by the experiment

Select one length-720, dimension-360 Golay--BA-3 code and reuse it in every
outer row. The parent parameters are

\[
 B=720,\qquad L=2944,\qquad
 K_0=1{,}059{,}840,\qquad N=2{,}119{,}680.
\]

Zero-shortening \(K_0-2^{20}=11{,}264\) input coordinates gives the requested
message dimension without increasing the set of codewords. The diagnostic
uses the parent expected spectrum, so it is conservative for that shortened
subcode. The bad-weight threshold is

\[
 D=\lfloor0.11N\rfloor=233{,}164.
\]

The random choices are:

1. the two Golay--BA-3 accumulator interleavers, sampled once;
2. one local coordinate permutation for each outer row;
3. one permutation inside each transposed region; and
4. the ideal random linear inner.

Only the first item selects the repeated outer code. The remaining choices
belong to the routing and inner maps.

## Exact occupation-one calculation

Fix a nonzero outer word of weight \(w\). Its local coordinate permutation
makes its support a uniform \(w\)-subset of the \(B\) regions. If \(J\) is its
first occupied region, then

\[
 \Pr[J=j]
 =\frac{\binom{B-j}{w-1}}{\binom Bw}.
\]

The region permutation makes the activation position \(U\) uniform on
\(\{1,\ldots,L\}\). The active suffix length is

\[
 M=(B-J)L+(L-U+1).
\]

Under the ideal inner, the output weight on that suffix has distribution
\(1+\operatorname{Bin}(M-1,1/2)\). Hence the exact binary64 single-row
quantity is

\[
 p_w=
 \sum_j \Pr[J=j]\frac1L\sum_{u=1}^L
 \Pr[\operatorname{Bin}(M(j,u)-1,1/2)\le D-1].
\]

Let \(A_w(C)\) denote the number of weight-\(w\) words in the one sampled BA
code. Grouping all row placements of the same word before averaging gives
the reuse-aware union bound

\[
 \Pr[\text{a bad }Q=1\text{ word}]
 \le
 \sum_{w=1}^B \mathbb E[A_w(C)]\min\{1,Lp_w\}.
\]

This grouping is valid without assuming independence between different row
events. It is sharper than multiplying every expected shell by \(L\) before
capping the final sum.

## Numerical result

The direct row-by-row union bound has 20.244610 bits of margin. The
reuse-aware bound has 24.305345 bits. Its leading terms are:

| Outer weight \(w\) | \(\log_2\mathbb E[A_w]\) | \(\log_2 p_w\) | Forced-late \(\log_2\) probability | Reuse-aware contribution \(\log_2\) |
|---:|---:|---:|---:|---:|
| 2 | -28.1630 | -4.38 | -6.39 | -28.1630 |
| 3 | -27.0521 | -6.57 | -9.60 | -27.0521 |
| 4 | -26.3599 | -8.78 | -12.84 | -26.3599 |
| 5 | -25.8411 | -10.99 | -16.09 | -25.8411 |
| 6 | -25.4184 | -13.21 | -19.36 | -27.11 |

For weights 2 through 5, \(Lp_w\ge1\). The row-placement cap therefore
removes the factor \(L\), but the expected multiplicities of the rare outer
shells themselves remain only about 26 to 28 bits below one.

The forced-late probability counts activations whose remaining suffix has at
most \(D\) coordinates. Every zero-state causal linear inner fails on this
event because it cannot emit nonzero symbols before activation and the suffix
cannot contain more than \(D\) nonzero symbols. This identifies a genuinely
causal part of the obstruction rather than a weakness of RM2Sub-S19.

Deleting low-weight terms from the same first moment gives:

| Retained outer weights | Reuse-aware \(Q=1\) margin |
|---|---:|
| \(w\ge10\) | 34.43 bits |
| \(w\ge11\) | 36.46 bits |
| \(w\ge12\) | 38.51 bits |
| \(w\ge13\) | 40.58 bits |
| \(w\ge14\) | 42.66 bits |

This deletion is an algebraic diagnostic. It is not a theorem conditioned on
minimum distance at least 13; such a theorem needs a specified sampling or
verification rule and the corresponding conditional spectrum bound.

## Validation and remaining obligations

The evaluator computes SciPy binomial tails and all spectrum sums in nearest
binary64. A separate exhaustive check at \(B=5,L=4,D=5,w=2\) agrees with the
implemented first-activation formula to \(1.2\times10^{-16}\).

The receipt is ba3_B720_ideal_causal_g1_q1_d11.json, generated by
evaluate_ba_ideal_causal_inner_q1.py.

The subsequent prefix-rank analysis removed two apparent obligations. It
covers all occupations without multiplying BA spectra, and it does not
require a low-shell hypothesis for a fixed setup.

FINITE_K20_RANDOM_CONV_PREFIX_ANALYSIS.md records the result. Independent BA
rows retain 20.244610 bits over all messages. Four fixed repeated-BA setups
at \(B=720\) retain 190.521692 bits. Four more at \(B=240\) retain the same
margin.

The first two former obligations are closed. One obligation remains:

1. Benchmark Toeplitz multiplication or design a linear-time inner that
   retains enough of its suffix law.
