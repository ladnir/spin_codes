# Dense groups with mostly weight-two columns

This restricted certificate covers a family excluded by the typical-histogram
argument. It does **not** establish distance for the full grouped SPIN code.
The target remains K=2^20, N=2^21, and output weight at most 209715.

Fix four BCH codewords a,b,c,d in an active group. Require a+b+c+d=0,
where addition is XOR. Every column then has weight zero, two, or four.
Call a column exceptional when its weight is not two. The certificate
restricts every active group to at most e exceptional columns.

## Counting this family

Let C be the implemented binary [256,128] code. Set f=a+b and h=1+a+c,
where 1 is the all-ones vector. The ordered tuple is uniquely determined by
(a,f,h), with a,f in C and h in the coset 1+C. A column is exceptional
exactly when f is zero and h is one there. Consequently its exception
count is the weight of h outside the support of f.

Write A_v for a coefficientwise upper bound on the number of weight-v
codewords. Let d_v bound the dimension of C shortened to any v coordinates.
For fixed f, restriction of 1+C outside its support has fibers of size
either zero or 2^dim(C_support(f)). Hence the number of eligible tuples is
at most

    C_e = min(2^384,
              2^128 sum_v A_v 2^d_v sum_{j=0}^{min(e,256-v)} binom(256-v,j)).

This is an upper count, not an exact enumerator. When e<256, none of these
tuples is the all-zero group.

## A direct lane-shuffle bound

Fix the coordinate and group permutations and the inner setup. Keep the
independent lane shuffles random. At any weight-two column, its four input
bits X are uniform among the six weight-two patterns. For every fixed
four-bit offset a and every 0<=z<=1,

    E[z^wt(a+X)] <= m_2(z) := (1+4z^2+z^4)/6.

The possible moments are z^2, (z+z^3)/2, and m_2(z), according to the
offset's weight. Their differences from m_2 are respectively
(1-z^2)^2/6 and (1-z)^2(z^2-z+1)/6, both nonnegative.

Condition on the past before each inner epoch. Its state is then fixed,
and its output is x+E(state), before feedback updates the state. Independent
current lane shuffles therefore give one factor m_2 per weight-two window.
All other windows contribute at most one. Iterating over epochs proves
the bound without assuming that the state is uniform. In particular it
applies to the target two-update inner.

There are 2048 possible active groups. For q active groups satisfying the
exception restriction, the union contribution is at most

    U_q(z) = binom(2048,q) C_e^q z^(-209715) m_2(z)^(q(256-e)).

Each q uses its own positive rational z<=1. A floating calculation proposes
z only; Arb evaluates the bound and sums occupancies outward.

## Complementary shortening bounds

For any coordinate set S of a binary [n,k] code,

    dim C_S = |S|-(n-k)+dim (C^perp)_(S^c).

Indeed, puncturing the dual to S has rank |S|-dim C_S, and its kernel is
the dual shortened to S^c. Thus a uniform dual shortening bound improves
d_v to min(d_v, v-128+d_dual(256-v)).

`dual_shortening.py` replays the complete Fourier-case certificate giving
dual distance at least 30. It also reads the production generator and
checks its rank, Q<=C<=P, and all-ones containment. The last check makes
C^perp even, allowing the existing exact even-code LP witnesses. Rejected
numerical proposals contribute no bound. Hamming and Griesmer bounds
remain available at every length. Exhaustive small-code tests check 432
complementary-support identities, including codes without all ones.

Combined with the existing positive-polynomial refinement, selected
shortening dimensions improve as follows:

| Support size | Previous cap | With dual complement |
|---|---:|---:|
| 128 | 47 | 47 |
| 160 | 77 | 66 |
| 192 | 107 | 72 |
| 200 | 115 | 75 |

## Restricted certificate results

The following use both positive-polynomial and dual-complement shortening.
Each row sums all occupancies in its stated range; every active group must
satisfy that row's exception restriction. The classes overlap, so their
margins must not be added.

| Exceptions per group | log2 tuple-count upper | Active groups covered | Summed margin, bits |
|---|---:|---|---:|
| 0 | 315.6150 | 1587--2048 | 174.9647 |
| 4 | 338.0151 | 1815--2048 | 95.6561 |
| 8 | 354.0938 | 2037--2048 | 150.1702 |

Independent 192- and 384-bit runs verify these bounds. The program restores the requested
Arb precision after importing legacy modules, which can change its global
context. The new proof components also pass 134 exact lane-moment and
small-code count tests. Continuous validity of the lane bound follows
from the polynomial identities above, not just the sampled test values.

Reproduce from the repository root:

    python -B research/workstreams/permutation_locality/even_column_certificate.py --dual-shortening --positive-shortening --precision 384 --exceptions 0 4 8
    python -B research/workstreams/permutation_locality/even_column_certificate.py --dual-shortening --positive-shortening --precision 192 --exceptions 0 4 8
    python -B research/workstreams/permutation_locality/shortening_moments.py --dual-shortening --positive-polynomials

The complementary dimensions also improve the general four-row support
counts through containment moments. Log2 upper counts change from
350.9023 to 348.9427 at support 128, from 431.2912 to 396.1708 at 160,
and from 461.9518 to 419.8266 at 176. These comparisons use the original
basis-count bounds as the baseline. They are exact integer count bounds;
the displayed logarithms are rounded diagnostics.

The general intermediate regime and the remaining atypical dense
histograms remain open. Next, use the new dimensions in the occupancy
cover and measure their effect before adding further state refinements.
Production and the paper are unchanged.
