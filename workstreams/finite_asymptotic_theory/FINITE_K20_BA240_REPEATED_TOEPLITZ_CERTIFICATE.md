# Finite \(k=2^{20}\) repeated-BA Toeplitz certificate

## Certified claim

There is an explicit binary linear code construction with

\[
 k=2^{20},\qquad N=2{,}119{,}680
\]

such that, over one random lower-triangular Toeplitz kernel,

\[
 \Pr[d_{\min}\le 233{,}164] < 2^{-180}.
 \tag{1}
\]

Consequently, except with probability below \(2^{-180}\), the code has

\[
 d_{\min}\ge233{,}165,
 \qquad
 \frac{d_{\min}}N\ge0.1100000943>0.11.
\]

The outward computation gives 190.521691555 bits of margin. The theorem
claims the smaller integer margin of 180 bits.

This is a complete finite-distance certificate for the construction defined
below. It is not a linear-time certificate: full-length Toeplitz
multiplication still needs an implementation and complexity decision.

## Construction

Set

\[
 B=240,\qquad L=8832,\qquad N=BL,
 \qquad K_0=120L=1{,}059{,}840.
\]

The local base code is the direct sum of ten extended binary Golay
\([24,12,8]\) codes. In each Golay block, the verifier uses the cyclic
generator

\[
 g(x)=x^{11}+x^9+x^7+x^6+x^5+x+1
\]

and appends the parity coordinate at position 23. It checks the exact local
weight enumerator

\[
 1+759z^8+2576z^{12}+759z^{16}+z^{24}.
\]

Apply two permute-then-accumulate maps to the 240 generator columns. This
produces one \([240,120]\) Golay--BA-3 code. Reuse this same code in all
\(L\) rows. Each row has a fixed local coordinate order. After transpose,
each of the 240 regions has a fixed row order.

All these fixed permutations come from the seed

\[
 \mathtt{0x4241323430544f45}.
\]

The seed is not random in the theorem. It is a compact specification of one
deterministic setup. The receipt states SplitMix64, rejection sampling, the
descending Fisher--Yates convention, the order in which the stream is
consumed, and hashes of every resulting permutation stream.

To obtain dimension \(k=2^{20}\), retain the first \(2^{20}\) row-major
parent input coordinates and set the remaining 11,264 parent inputs to zero.
This zero-shortening selects a subcode and therefore cannot reduce minimum
distance.

Finally, sample independent uniform bits

\[
 h_1,\ldots,h_{N-1}\gets\mathbb F_2,
 \qquad h_0:=1,
\]

and map the routed outer word \(u\) to

\[
 y_t=\sum_{i=0}^t h_i u_{t-i},
 \qquad 0\le t<N.
 \tag{2}
\]

The same Toeplitz kernel is shared by every message.

## Probability space

The outer generator, both accumulator permutations, all row-local coordinate
orders, all region row orders, and the zero-shortening rule are fixed. The
only random variables are the \(N-1\) independent bits
\(h_1,\ldots,h_{N-1}\).

The failure event in (1) is the existence of a nonzero shortened message
whose output under (2) has Hamming weight at most 233,164. The verifier proves
the stronger statement for every nonzero message in the unshortened
\(K_0\)-dimensional parent code.

## Toeplitz suffix law

Fix a nonzero routed outer word whose first one is at coordinate \(\tau\).
Equation (2) gives zero before \(\tau\) and one at \(\tau\). For \(r\ge1\),

\[
 y_{\tau+r}
 =h_r+\sum_{i=0}^{r-1}h_i u_{\tau+r-i}.
\]

The map from the fresh kernel bits to the output suffix is triangular with
unit diagonal. Hence

\[
 \operatorname{wt}(y)
 \mathrel{\overset{\mathrm d}=}
 1+\operatorname{Bin}(N-\tau-1,1/2).
 \tag{3}
\]

No independence between different messages is asserted or required.

## Complete prefix identity

Let \(V\) be the fixed routed parent code. Define

\[
 \kappa_t=\dim\{v\in V:v_0=\cdots=v_{t-1}=0\},
 \qquad Z_t=2^{\kappa_t}-1.
\]

Thus \(Z_t-Z_{t+1}\) counts the nonzero words first activating at \(t\).
For \(D=233{,}164\), (3) and summation by parts give

\[
 \mathbb E_h[Z_{\mathrm{bad}}]
 =Z_0F_N+
 \sum_{t=1}^{N-D}
 Z_t\frac{\binom{N-t-1}{D-1}}{2^{N-t}},
 \tag{4}
\]

where

\[
 F_N=\Pr[\operatorname{Bin}(N-1,1/2)\le D-1].
\]

Equation (4) covers all nonzero parent messages. It does not enumerate
messages and does not split them by active-row occupation.

## Outward bound

The verifier rebuilds every row-local routed generator column and computes
its exact rank increments by Gaussian elimination over \(\mathbb F_2\).
Different rows use disjoint parent input coordinates, so summing these local
rank increments gives every global \(\kappa_t\) exactly. For a prefix interval
\([a,b]\), set

\[
 C=\max_{a\le t\le b}(\kappa_t+t).
\]

Since \(Z_t<2^{\kappa_t}\), the interval contribution to (4) is at most

\[
 2^{C-N}
 \left(
   \binom{N-a}{D}-\binom{N-b-1}{D}
 \right).
 \tag{5}
\]

The verifier uses one interval per transposed region. The quantity
\(\kappa_t+t\) is nondecreasing, so the region endpoint gives \(C\). The first
prefix-rank deficit occurs at \(t=10{,}081\), after the first region. The
kernel becomes zero at \(t=1{,}163{,}547\).

For the base term, the ratio of consecutive lower-tail binomial masses gives

\[
 F_N\le
 \frac{\binom{N-1}{D-1}}{2^{N-1}}
 \frac{N-D+1}{N-2D+2}.
 \tag{6}
\]

The verifier evaluates (5) and (6) with 256-bit python-flint Arb intervals.
It obtains

\[
 \log_2\mathbb E_h[Z_{\mathrm{bad}}]
 <-190.5216915553056652.
\]

Markov's inequality now gives (1). All integer thresholds are fixed before
the probability is evaluated.

## Reproduction and authenticated values

Run:

```powershell
python workstreams/finite_asymptotic_theory/certify_ba240_repeated_toeplitz_prefix_outward.py
```

The verifier writes
`ba3_B240_repeated_toeplitz_prefix_outward.json`. Important authenticated
values are:

- accumulator permutation hashes:
  `2d2cba5b7c03d751ff6cbc9ffa9dcbc52934e46c906225065510c9b3842628ae`
  and
  `42b4c26f55f0d2df7166077c11abf31753aac0712415d8e2b91783b05eccf366`;
- row-local permutation-stream hash:
  `985b15c90daaf15b8aaa048b5a2ddd96533a7f8a83da9659b2d4640180395647`;
- region permutation-stream hash:
  `d95b8118118c9b7b73dbcd93ebef14554682a418a32f5335e788eb5224a12f89`;
  and
- complete prefix-dimension hash:
  `aef6b0dc49c858a5b13053b741a1c6a50a309747233e7eacfcfb9b09167c9632`.

## Scope and remaining obligations

The following statements are proved for the specified finite construction:

1. the outer parent and shortened code dimensions;
2. the exact routed prefix dimensions;
3. the Toeplitz suffix law for every fixed nonzero input;
4. the all-message first-moment identity; and
5. the 180-bit finite 11% distance-failure bound.

The following statements are not proved here:

1. a probability of finding this setup by sampling a BA code and routes;
2. a linear-time bound for full-length Toeplitz multiplication;
3. competitive ordinary or transposed implementation performance; and
4. replacement of Toeplitz convolution by RM2Sub-S19 or another bounded-state
   inner while preserving this margin.

The certificate therefore resolves the finite proof-model question that
motivated the experiment. It also isolates the remaining issue: the current
outer is adequate once fixed, while the Toeplitz suffix law may be too costly
to realize directly.
