# Why random convolution closes for a typical repeated BA outer

## Conclusion

Random lower-triangular Toeplitz convolution closes the finite 11% distance
calculation for the sampled repeated-BA setups tested here. The calculation
covers every nonzero parent message at once. A subsequent verifier froze one
\(B=240\) setup and converted the calculation into an outward certificate
with 180 claimed bits of margin. See
`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md`.

At the current \(B=240\) design point, four fixed repeated-BA setups each
give 190.521692 binary64 bits of first-moment margin. Four \(B=720\) setups
give the same margin. These values are essentially the best possible margin
for the parent rate-one-half dimensions.

The corresponding unconditioned independent-row BA ensembles give only
10.623420 bits at \(B=240\) and 20.244610 bits at \(B=720\). Almost all of
that ensemble contribution comes from occupation one. Rare BA interleavers
create weight-2 through weight-5 outer words. Typical sampled BA codes do not
contain those words.

Indeed, Markov's inequality bounds the probability that one \(B=720\) BA
sample contains a word of weight 2 through 5 by \(2^{-24.621}\). At \(B=240\),
the analogous bound is \(2^{-16.285}\).

The random-outer comparison closes because a uniform random
\([720,360]\) code suppresses weight-2 through weight-5 words by hundreds of
additional bits. Agreement between the BA and random spectra near weight
\(B/2\) is irrelevant to this finite failure mode.

The new evidence changes the proof target. Increasing \(B\) is not necessary
for the tested fixed setups. A better route is to freeze or accept one
repeated BA code and its routing permutations after checking their prefix
ranks. The remaining construction question is whether Toeplitz convolution,
or a linear-time approximation to it, meets the performance target.

All numerical claims generated directly in this document use nearest
binary64. They are diagnostics. The separate frozen-setup certificate uses
256-bit outward Arb arithmetic.

## The random convolution

Fix the routed outer word

\[
 u=(u_0,\ldots,u_{N-1})\in\mathbb F_2^N.
\]

Setup samples independent uniform bits

\[
 h_1,\ldots,h_{N-1}\gets\mathbb F_2
\]

and sets \(h_0:=1\). The inner map is the causal Toeplitz convolution

\[
 y_t:=\sum_{i=0}^t h_i u_{t-i},
 \qquad 0\le t<N.
 \tag{1}
\]

The map is linear and invertible because its matrix is lower triangular with
unit diagonal.

Suppose \(u_\tau=1\) and \(u_t=0\) for \(t<\tau\). Then \(y_t=0\) for
\(t<\tau\), and \(y_\tau=1\). For \(r\ge1\),

\[
 y_{\tau+r}
 =h_r+\sum_{i=0}^{r-1}h_i u_{\tau+r-i}.
\]

Thus the map from \((h_1,\ldots,h_{N-\tau-1})\) to
\((y_{\tau+1},\ldots,y_{N-1})\) is triangular with unit diagonal. For this
fixed nonzero \(u\), the later output bits are independent and uniform.
Consequently,

\[
 \operatorname{wt}(y)
 \mathrel{\overset{\mathrm d}=}
 1+\operatorname{Bin}(N-\tau-1,1/2).
 \tag{2}
\]

The same convolution coefficients serve every message. The first-moment
argument below does not require independence between different messages.

## An all-message identity

Fix one outer code and all routing permutations. Let \(V\) be the resulting
parent code before convolution. For \(0\le t\le N\), define

\[
 \kappa_t:=
 \dim\{v\in V:v_0=\cdots=v_{t-1}=0\}
\]

and

\[
 Z_t:=2^{\kappa_t}-1.
\]

Thus \(Z_t\) counts the nonzero routed outer words whose first \(t\)
coordinates are zero. The number whose first nonzero coordinate is exactly
\(t\) equals \(Z_t-Z_{t+1}\).

Set

\[
 D:=\lfloor0.11N\rfloor
\]

and let \(Z_{\mathrm{bad}}\) count nonzero parent messages whose final weight
is at most \(D\). Equation (2) gives the exact identity

\[
 \mathbb E_h[Z_{\mathrm{bad}}]
 =
 \sum_{t=0}^{N-1}
 (Z_t-Z_{t+1})
 \Pr[\operatorname{Bin}(N-t-1,1/2)\le D-1].
 \tag{3}
\]

Summation by parts gives the computational form

\[
 \mathbb E_h[Z_{\mathrm{bad}}]
 =
 Z_0F_N+
 \sum_{t=1}^{N-D}
 Z_t
 \frac{\binom{N-t-1}{D-1}}{2^{N-t}},
 \tag{4}
\]

where

\[
 F_N:=
 \Pr[\operatorname{Bin}(N-1,1/2)\le D-1].
\]

Equation (4) counts all \(2^{K_0}-1\) nonzero parent messages. It does not
split them by occupation. Markov's inequality turns the right-hand side into
a distance-failure bound over the convolution coefficients.

The prefix dimensions \(\kappa_t\) are computable by Gaussian elimination on
the routed generator matrix. This calculation does not enumerate messages.

## Prefix counts for independent outer rows

The independent-row models admit a second exact calculation. For a local
code \(C\subseteq\mathbb F_2^B\), define

\[
 M_r(C):=
 \frac1{\binom Br}
 \sum_{\substack{S\subseteq[B]\\|S|=r}}
 \left|\{c\in C:c|_S=0\}\right|.
 \tag{5}
\]

If \(A_w(C)\) is the local weight spectrum, then

\[
 M_r(C)
 =
 \sum_{w=0}^{B-r}
 A_w(C)
 \frac{\binom{B-w}{r}}{\binom Br}.
 \tag{6}
\]

Let \(t=jL+u\), where \(0\le u<L\). A routed prefix of length \(t\) contains
\(j+1\) coordinates from \(u\) outer rows and \(j\) coordinates from the
other rows. If the local codes are sampled independently, then

\[
 \mathbb E[Z_t]
 =
 \overline M_{j+1}^{\,u}
 \overline M_j^{\,L-u}-1,
 \qquad
 \overline M_r:=\mathbb E_C[M_r(C)].
 \tag{7}
\]

Equation (7), inserted into (4), is exact for independent rows. For one
sampled code repeated in every row, the required quantity is instead

\[
 \mathbb E_C[
 M_{j+1}(C)^uM_j(C)^{L-u}].
 \tag{8}
\]

The expected spectrum determines \(\overline M_r\), but it does not
determine the high moment in (8). Replacing (8) by the powers in (7) is the
precise invalid step that a repeated-code proof must avoid.

For a fixed repeated code and fixed routing, equation (4) avoids (8)
entirely. It uses the actual prefix ranks.

## Results at \(B=720\)

The parent parameters are

\[
 B=720,\quad L=2944,\quad
 K_0=1{,}059{,}840,\quad
 N=2{,}119{,}680,\quad
 D=233{,}164.
\]

Zero-shortening 11,264 parent input coordinates gives \(k=2^{20}\). Using
the parent code in the calculation is conservative.

| Outer model | All-message margin | Dominant mechanism |
|---|---:|---|
| Independent uniform \([720,360]\) injections | 189.230985 bits | Full-length and earliest-prefix terms |
| Independent Golay--BA-3 rows | 20.244610 bits | Rare low-weight BA words at a late prefix |
| Ideal full-rank prefix profile | 190.521692 bits | Full-length and earliest-prefix terms |
| Fixed repeated BA, seed 20260902 | 190.521692 bits | Full-length and earliest-prefix terms |
| Fixed repeated BA, seed 20260903 | 190.521692 bits | Full-length and earliest-prefix terms |
| Fixed repeated BA, seed 20260904 | 190.521692 bits | Full-length and earliest-prefix terms |
| Fixed repeated BA, seed 20260905 | 190.521692 bits | Full-length and earliest-prefix terms |

The independent-BA occupation-one margin is 20.244610 bits. Subtracting its
contribution from the all-message result leaves 41.491 bits for every
occupation \(Q\ge2\) combined. Thus occupation one accounts for all but about
\(4.02\times10^{-7}\) of the unconditioned independent-BA first moment.

## The critical prefix

The independent-BA calculation is largest at

\[
 t=1{,}653{,}351
   =0.780000283\,N.
\]

The remaining suffix has length

\[
 N-t=466{,}329
     =0.219999717\,N.
\]

Half of this suffix is approximately \(0.11N\). Therefore the convolution
has about one-half probability of landing below the distance threshold once
an input first activates near this location.

More precisely, the conditional expected output weight is
\(1+(466{,}329-1)/2=233{,}165=D+1\).

At \(B=720\), the critical prefix contains 561 complete regions and part of
the next region. A one-row outer word survives the complete regions only if
its support lies in the remaining 159 coordinates.

For the BA ensemble, the expected number of nonzero local messages that
survive 561 uniform coordinate constraints is

\[
 2^{-31.753866}.
\]

Multiplying by \(L=2944\) row locations gives approximately
\(2^{-20.230304}\). This is the main finite loss.

For a uniform random \([720,360]\) injection, the corresponding local
expectation is exactly

\[
 \frac{(2^{360}-1)(2^{159}-1)}{2^{720}-1}
 \approx 2^{-201}.
\]

The local gap is 169.246 bits. This is also the approximate gap between the
random-outer and BA margins.

Weights 2 through 5 contribute 98.806% of the BA survivor expectation at
this prefix. Weight 2 alone contributes about 58.5%.

## Why the central spectrum comparison was misleading

At weight 2, the expected BA multiplicity is

\[
 \log_2\mathbb E[A_2^{\mathrm{BA}}]=-28.1630.
\]

For a uniform random \([720,360]\) injection, it is

\[
 \log_2\mathbb E[A_2^{\mathrm{rnd}}]=-342.0183.
\]

At weight 5, the corresponding values are \(-25.8411\) and \(-319.4677\).
The two ensembles can therefore look almost identical around weight 360 and
remain radically different for the finite distance proof.

The relevant statistic is not a bulk spectral norm. It is the low-fugacity
transform in (6), evaluated when roughly 78% of the coordinates are forced
to zero. That transform gives extreme low weights most of its mass.

The observed block-size sweep reinforces this explanation.

| \(B\) | \(L=N/B\) | Independent random margin | Independent BA margin | Reuse-aware \(Q=1\) margin |
|---:|---:|---:|---:|---:|
| 240 | 8832 | 54.044 | 10.623 | 15.456 |
| 360 | 5888 | 88.213 | 14.187 | 18.814 |
| 480 | 4416 | 122.214 | 16.706 | 21.127 |
| 720 | 2944 | 189.231 | 20.245 | 24.305 |
| 960 | 2208 | 190.522 | 22.750 | 26.513 |
| 1440 | 1472 | 190.522 | 26.275 | 29.720 |

Doubling \(B\) adds about six bits to the unconditioned BA margin in this
range. This behavior is consistent with polynomial suppression of the
lowest BA shells. A random linear code suppresses the same shells
exponentially in \(B\).

## Results at the current \(B=240\) design point

The current design uses \(B=240\) and \(L=8832\), with the same values of
\(K_0\), \(N\), and \(D\). Four independently seeded fixed repeated-BA
setups again give 190.521692 bits.

The independent-row random outer gives 54.043528 bits. The independent-row
BA ensemble gives 10.623420 bits. Thus the current block size already closes
the random-convolution model with more than ten ensemble-average bits under
the independent-row interface, and a typical fixed repeated setup behaves
much better.

Increasing \(B\) from 240 to 720 improves the probability of avoiding rare BA
interleavers. It does not improve the four tested fixed-code prefix bounds.

## Correct proof interface for one repeated code

The repeated-code result should use two levels.

1. Select or freeze one BA generator and the routing permutations.
2. Compute the complete prefix-rank sequence \((\kappa_t)_{t=0}^N\).
3. Reject the setup unless an outward evaluation of (4) meets the target.
4. Sample one Toeplitz convolution and pay the conditional first moment in
   (4).

This setup test uses linear algebra on a \(360\)-by-\(720\) local generator
and its routed copies. It never enumerates \(2^{360}\), \(2^{128}\), or
\(2^{K_0}\) messages.

The four PCG64 seeds are only diagnostics. The follow-up certificate uses a
different, stable SplitMix64/Fisher--Yates setup specification and binds the
selected generator, permutations, prefix ranks, and arithmetic receipt by
hashes.

The deterministic frozen setup avoids a setup-probability claim. A randomized
setup theorem would still need an acceptance-probability bound.

## Independent validation

The evaluator reconstructs the extended Golay constituent from its cyclic
generator. Exhaustive enumeration of its 4096 messages gives

\[
 1+759z^8+2576z^{12}+759z^{16}+z^{24}.
\]

A separate \(N=9\) enumeration checked the Toeplitz suffix law. For a fixed
input activating at coordinate 3, all 32 permitted output suffixes occur
exactly eight times as the 256 convolution tails vary.

A separate small-code enumeration checked equation (3). The direct
all-message expectation and the prefix formula agree to
\(8.9\times10^{-16}\).

Finally, an 80-decimal-digit calculation gives 192.706121979414 bits for the
full-length term. The binary64 evaluator gives 192.706121982403 bits, a
difference of \(3.0\times10^{-9}\) bits.

## What is proved, conditional, and open

The following mathematical facts are proved.

1. Random Toeplitz convolution has the suffix law in (2).
2. Equations (3) and (4) are exact for every fixed routed outer code.
3. Equations (6) and (7) are exact for independently sampled outer rows.
4. A fixed prefix-rank sequence can be checked without message enumeration.
5. For the frozen \(B=240\) setup, the all-message expected bad count is below
   \(2^{-180}\) in outward arithmetic.

The following numerical facts are diagnostics.

1. Independent BA rows give 10.623420 bits at \(B=240\).
2. Independent BA rows give 20.244610 bits at \(B=720\).
3. Each of eight tested fixed repeated-BA setups gives 190.521692 bits.
4. Occupations \(Q\ge2\) contribute only \(2^{-41.491}\) in the independent
   \(B=720\) model.

The following obligations remain open.

1. Benchmark fast binary Toeplitz multiplication at \(N=2{,}119{,}680\).
2. If strict linear time is required, design and analyze a bounded-state
   convolution whose all-message transfer retains enough of the Toeplitz
   margin.

`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md` records the frozen setup,
outward proof, and exact shortening rule. The next task is implementation:
benchmark one forward and one transposed multiplication by a fixed binary
Toeplitz matrix.
