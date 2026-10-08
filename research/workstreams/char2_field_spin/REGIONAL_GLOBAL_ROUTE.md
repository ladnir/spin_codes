# A proved dense class for the regional ensemble

The whole-code proof can be split at inverse-input weight. The dense class admits an upper bound that does not analyze the inner state or count source words. At the current sizes, the bound below gives more than forty bits for every independently chosen invertible inner.

This is not a distance certificate. The sparse and intermediate classes remain open for the revised candidates. In particular, this result does not rescue the rejected GRS[32,24], state-four regional construction: its window witness lies outside the class proved here.

## Ensemble and the event being bounded

Let \(C\subseteq\mathbb F_q^b\) be a fixed linear MDS code of dimension \(k\), and put \(p=b-k\). The outer code consists of \(L\) independent rows of \(C\). Its length is \(N=bL\).

For each outer coordinate, sample an independent permutation of the \(L\) row labels. Together these permutations form the single regional routing stage \(P\). Independently sample a nonzero field multiplier for every coordinate; their diagonal matrix is \(D\). Let \(B\) be an invertible inner map, independent of both \(P\) and \(D\). The encoded code is

\[
 BPD(C^{\oplus L}).
\]

All maps are linear over \(\mathbb F_q\). The inner may be random, provided conditioning on it preserves the specified laws of \(P,D\).

Fix an output vector \(y\ne0\), and write \(x=B^{-1}y\). Its inverse-input density is \(\alpha=\operatorname{wt}(x)/N\). The proof first bounds the probability that this fixed \(x\) lies in \(PD(C^{\oplus L})\). It then sums over low-weight output lines whose inverse inputs meet a density cutoff.

For the GF256 constructions, proving distance over GF256 also proves it over every compatible scalar extension. A nonzero component of an extension-field word is a GF256 word whose support is contained in the original support. No transitivity of GF256 diagonals on larger-field nonzero symbols is assumed.

## Membership bound, uniform over regional weight profiles

Let \(u_j\) be the weight of \(x\) in region \(j\), and set \(\alpha_j=u_j/L\). Thus \(b^{-1}\sum_j\alpha_j=\alpha\).

For a prescribed row support \(S\subseteq[b]\), let \(A_S\) count the outer words with exactly that support. Independent nonzero diagonals make the values on \(S\) uniform over \((\mathbb F_q^*)^S\). The row-membership probability is therefore

\[
 f(S)=\frac{A_S}{(q-1)^{|S|}}.
\]

We have \(f(\varnothing)=1\). MDS shortening gives \(A_S=0\) when \(1\le |S|\le p\), and \(A_S\le q^{|S|-p}\) otherwise. Consequently, for every nonempty \(S\),

\[
 f(S)\le q^{-p}\left(\frac q{q-1}\right)^{|S|}.
\]

The inverse route gives, in each region, an independent uniform \(u_j\)-subset of row labels. We compare this law with independent Bernoulli entries. For a positive tilt \(\lambda\), use column probabilities

\[
 \theta_j=\frac{\lambda\alpha_j}{1-\alpha_j+\lambda\alpha_j}.
\]

Conditioning each column on weight \(u_j\) recovers the required uniform subsets. At Bernoulli parameter \(\alpha_j\), that count is a binomial mode and has probability at least \(1/(L+1)\). Comparing the tilted count probability with the modal count probability yields

\[
 \Pr[x\in PD(C^{\oplus L})]
 \le (L+1)^b\lambda^{-\alpha N}
 \left[
   \prod_j(1-\alpha_j)
   +q^{-p}\prod_j\left(1-\alpha_j+
                \frac q{q-1}\lambda\alpha_j\right)
 \right]^L.
\]

The second product deliberately also counts the empty support; this only increases the bound. Applying arithmetic--geometric mean to each nonnegative product proves

\[
 \boxed{
 \Pr[x\in PD(C^{\oplus L})]
 \le (L+1)^b\inf_{\lambda>0}
 \lambda^{-\alpha N}
 \left[(1-\alpha)^b+
 q^{-p}\left(1-\alpha+\frac q{q-1}\lambda\alpha\right)^b
 \right]^L.}
\]

There is no union over regional profiles. The displayed inequality holds separately for every profile with mean density \(\alpha\).

An independent audit checked the tilted conditioning, the mode-probability loss, both applications of arithmetic--geometric mean, and the independence requirement on the inner.

## Finite, exact dense-class bounds

Take \(q=256\), target distance \(10\%\), and bad-weight threshold \(e=\lceil N/10\rceil-1\). The number of nonzero projective output vectors of weight at most \(e\) is at most \(2^{1.27N}\). One exact Chernoff verification is

\[
 \sum_{j=0}^{e}\binom Nj255^j
 \le (10/9)^N2295^{N/10},
 \qquad (10/9)^{100}2295^{10}<2^{127}.
\]

Dropping the projective normalization only weakens this bound.

The replay proves a uniform membership bound

\[
 \Pr[x\in PD(C^{\oplus L})]\le (L+1)^b2^{-1.28N}
\]

throughout each dense interval below.

| Outer | Rate | Covered inverse-input densities |
|---|---:|---:|
| GRS[32,24] | 3/4 | \(\alpha\ge0.30\) |
| GRS[33,24] | 8/11 | \(\alpha\ge0.28\) |
| GRS[32,23] | 23/32 | \(\alpha\ge0.27\) |
| GRS[32,22] | 11/16 | \(\alpha\ge0.24\) |
| GRS[34,24] | 12/17 | \(\alpha\ge0.25\) |
| GRS[64,48] | 3/4 | \(\alpha\ge0.30\) |
| GRS[36,24] | 2/3 | \(\alpha\ge0.25\) |

Hence the probability that **any** low-weight output has an inverse input in the stated dense class is at most

\[
 \boxed{(L+1)^b2^{-N/100}.}
\]

This is a bound on the existence of a bad codeword in that class. It is not an expected source-word count averaged over a rare route event.
The argument is still a union bound. It succeeds after removing the low-inverse-weight class; it does not repair the unrestricted first moment.

For an interval \([j/200,(j+1)/200]\), the replay chooses a rational \(\lambda\ge1\). It upper-bounds the decreasing first term at the left endpoint, the increasing second term at the right endpoint, and \(\lambda^{-\alpha N}\) at the left endpoint. Raising to the power 200 turns the claimed exponent into an exact rational inequality. Floating optimization only proposes the rational tilt; it cannot cause an inequality to pass.

With \(K\) rounded to a multiple of \(64k\), conservative integer margin lower bounds are:

| Outer | K16 | K18 | K20 |
|---|---:|---:|---:|
| GRS[32,24] | 496 | 3054 | 13475 |
| GRS[33,24] | 512 | 3149 | 13896 |
| GRS[32,23] | 537 | 3217 | 14090 |
| GRS[32,22] | 578 | 3381 | 14745 |
| GRS[34,24] | 527 | 3244 | 14318 |
| GRS[64,48] | 197 | 2690 | 13048 |
| GRS[36,24] | 558 | 3435 | 15160 |

These are dense-class margins only. The script computes the rounded lengths and verifies each displayed class has at least forty bits.

## What is still needed for the complete code

The remaining event has \(0<\operatorname{wt}(B^{-1}y)<\alpha_0N\), where \(\alpha_0\) is the chosen cutoff. The previous window obstruction shows why this part cannot use an unconditional source-word first moment: a rare route can create a high-dimensional bad subspace.

A valid hybrid decomposition is

\[
 \Pr[\mathrm{bad}]
 \le \Pr_P[P\notin\mathcal G]
 +\mathbb E_{P,B}\!\left[
  \mathbf1_{P\in\mathcal G}
  \min\{1,F_{\mathrm{sparse}}(P,B)\}\right]
 +(L+1)^b2^{-N/100},
\]

where \(\mathcal G\) is a specified good-route event and \(F_{\mathrm{sparse}}\) is a valid conditional union bound over the remaining output lines, with only \(D\) still random. Clipping must occur before the expectation over routes. A complete proof must define \(\mathcal G\), cover every sparse bad word, and evaluate the first two terms. Checking a finite family of rectangular windows is insufficient.

For identity-feedthrough inners, increasing the injection period does not add constraint capacity. Increasing state dimension or adding outer parity changes that capacity, but neither change alone supplies the missing covering theorem.

### A more useful local transfer

A superregular map on \(t+s\) input/state coordinates has branch distance \(t+s+1\). If a nonzero branch has \(w\) nonzero input symbols and \(v\) nonzero output symbols, then

\[
 v\ge\max\{0,t-s+1-w\}.
\]

This removes low-input, low-output local branches deterministically. It is a genuine change from identity feedthrough. However, local support inequalities alone do not control dense packet inputs: the formal pattern of all input and state coordinates nonzero, with few output coordinates nonzero, satisfies those inequalities. Such a pattern is not asserted to be a codeword; it demonstrates the incompleteness of support-only screening.

With an independent uniform \(\mathrm{GL}(s,q)\) refresh before each packet and fresh nonzero input diagonals, a superregular local map also gives a tractable **fixed-word** transfer. For \(w\) nonzero inputs and \(r\) forced-zero output/state coordinates, the number of assignments with nonzero incoming state is

\[
 \sum_{a=0}^{w}(-1)^a\binom wa
 \left[q^{\max(0,w-a+s-r)}-q^{\max(0,w-a-r)}\right].
\]

Normalize by \((q-1)^w(q^s-1)\). For zero incoming state, use only the second power without subtraction and normalize by \((q-1)^w\). Forcing all next-state coordinates separates the zero-state transition. Inclusion--exclusion over the surviving output coordinates then gives exact output weights.

This formula was independently audited as a fixed-word law. It does not make the action on a surviving message subspace uniformly random. Its most credible use is a hybrid: deterministic sparse-branch exclusions and conditional route clipping for the low-weight class, this exact transfer where a value union is affordable, and the proved dense bound above. No whole-code conclusion follows until the low-weight contribution is evaluated.

## Replay

Run the mathematical reference only:

```sh
python -B research/workstreams/char2_field_spin/theory/regional_global_route.py --self-test
```

The reference checks 256 exact GF4 profile/tilt cases, then verifies every density interval using rational arithmetic. `--show-tilts` emits all rational witnesses. No native implementation, remote execution, or benchmark is included.
