# Blockwise SPIN: capped support composition

## Result and remaining gap

The blockwise local operator now gives a complete restricted-class bound.
For GRS[128,64], \(K=2^{20}\), packet width64, and state16 over
\(\mathrm{GF}(65536)\), the probability of a word below 40% distance whose
message occupies exactly one outer row is less than \(2^{-135}\).
This includes every choice of that row, every coordinate support, and every
field message in the class. It is not a whole-code distance certificate.

The subsequent [quiet-subcode argument](HD_BLOCKWISE_QUIET_SUBCODE.md)
deterministically rejects the width64/state4 and width64/state8 controls at
40%, for both outer lengths. Their singleton bounds below remain correct,
but those constructions are not viable whole-code candidates. Fix any
\(a\) source rows. The \(ak\)-dimensional message space has a nonzero vector
with every preterminal packet syndrome zero whenever
\(ak>s(N/t-1)\). The state stays zero, so \(Y=X\) and
\(\operatorname{wt}Y/N\le a/R\). At the listed K20 geometries, this gives
12.5% for width64/state4, 25% for width64/state8, and 50% for
width64/state16. These are upper bounds on actual code distance for every
setup, not failures of a proof technique.

The method caps the conditional failure probability of each exact source
support before averaging over the regional route. A fractional-power
relaxation then composes with the positive local operator. This is different
from counting all bad codewords in an uncapped first moment.

Two limitations are also quantified. A route-only inverse-word bound loses
more than 2.29 million bits against the 40% output ball. Moreover, one common
fractional exponent still fails on the known multi-row routing traps.
The next many-row bound must retain additional source or routing constraints.

## Conditional cap for one exact source support

Use the exact maps and three-component state measure from
[HD_BLOCKWISE_MAPS.md](HD_BLOCKWISE_MAPS.md). Setup consists of independent
uniform regional permutations \(\Pi\) and independent packet mixers \(M\).
Each setup is fixed for all messages. There are no independent input
diagonals.

Let the outer be a direct sum of \(R\) copies of an MDS \([b,k]\) code over
\(\mathbb F_q\). Write \(p=b-k\). For one prescribed exact coordinate
support of size \(w\), the number of local codewords is

\[
 m_w=
 \begin{cases}
 1,&w=0,\\
 0,&1\le w\le p,\\
 (q-1)\displaystyle\sum_{j=0}^{w-p-1}
 (-1)^j{w-1\choose j}q^{w-p-1-j},&w>p.
 \end{cases}
\]

Fix a nonempty exact support pattern \(P\) across all outer rows, and let
\(w_a\) be its weight in row \(a\). Its number of source words is
\(B_P=\prod_a m_{w_a}\). Skip patterns for which \(B_P=0\). Global scalar
multiplication partitions the remaining words into \(B_P/(q-1)\) lines.

Fix the route \(\Pi\). It sends \(P\) to a sequence of packet weights
\(w_i(\Pi,P)\). Let \(T_w(\theta)\) be either valid local matrix from the
maps note, with the entrywise maximum over feasible syndrome flags. For
\(0<\theta<1\), every fixed source word with support \(P\) satisfies

\[
 \mathbb E_M\theta^{\operatorname{wt}(Y)}
 \le G_P(\Pi;\theta):=
 e_0\left(\prod_iT_{w_i(\Pi,P)}(\theta)\right)\mathbf1.
\]

The same matrix bound applies to all amplitudes on \(P\); the actual states
and syndromes need not agree. For an integer output threshold \(e\), the
conditional probability of at least one bad line on this support is at most

\[
 \min\left\{1,\frac{B_P}{q-1}\theta^{-e}G_P(\Pi;\theta)\right\}.
\]

This cap occurs before route averaging. It charges a rare route only once,
even when its setup admits many bad source words with the same support.

## Fractional-power local operator

Fix \(0<\rho\le1\). For \(x\ge0\), \(\min(1,x)\le x^\rho\). For a
finite family of nonnegative terms, \((\sum_jx_j)^\rho\le\sum_jx_j^\rho\).
Apply the second inequality to the state paths in the matrix product. With
entrywise powers \(\widehat T_w=T_w^{\circ\rho}\),

\[
 \Pr_{\Pi,M}[\text{a bad line has exact source support }P]
 \le
 \left(\frac{B_P}{q-1}\right)^\rho\theta^{-\rho e}
 \mathbb E_\Pi\left[
 e_0\left(\prod_i\widehat T_{w_i(\Pi,P)}\right)\mathbf1
 \right].
\]

No independence between different messages is used. The field-codeword
count receives the power \(\rho\); regional route probabilities do not.
The number of distinct coordinate supports also remains unpowered. In
particular, \(m_w^\rho\) is the prescribed-support coefficient. Replacing
\({b\choose w}m_w^\rho\) by \(({b\choose w}m_w)^\rho\) would be invalid.

The same derivation permits different witnesses \((\rho,\theta)\) for
different support classes. It does not justify choosing them after observing
the random route and then reusing an unconditional route average.

### Exact regional composition

Suppose a fixed set of \(h\) source positions is active in one region of
length \(R\), partitioned into \(B=R/t\) packets. Its routed support is a
uniform \(h\)-subset. Consequently its averaged operator is

\[
 Q_h=
 \frac1{{R\choose h}}[z^h]
 \left(\sum_{w=0}^t{t\choose w}z^w\widehat T_w\right)^B.
\]

Matrix multiplication preserves packet chronology. The coefficient counts
all support placements, including several active positions in the same
packet. This formula is valid without treating packet occupancies as
independent Bernoulli variables.

For a fixed set of \(A\) active source rows, one region contributes the
multivariate polynomial

\[
 K_A(z_1,\ldots,z_A)=
 \sum_{J\subseteq\{1,\ldots,A\}}Q_{|J|}\prod_{a\in J}z_a.
\]

The coefficient of \(z_1^{w_1}\cdots z_A^{w_A}\) in
\(e_0K_A^b\mathbf1\) averages all exact coordinate supports with those row
weights. Weighting it by \(\prod_a m_{w_a}^\rho\), then multiplying by
\({R\choose A}(q-1)^{-\rho}\theta^{-\rho e}\), bounds that active-row class.
This is a valid finite many-row interface. Its coefficient state is still
large; it is not presented as an efficient whole-code algorithm.

## Directed one-active-row certificate

For \(A=1\), the regional operators reduce to

\[
 Q_0=\widehat T_0^B,\qquad
 Q_1=\frac1B\sum_{j=0}^{B-1}
 \widehat T_0^j\widehat T_1\widehat T_0^{B-1-j}.
\]

Each of the \(R\) positions is equally likely. Every packet has \(t\)
positions, giving the unpowered factor \(1/B\). The local weight-one bound
already covers every position and nonzero amplitude inside that packet.

Let \(g_w=[z^w]e_0(Q_0+zQ_1)^b\mathbf1\). The complete class bound is

\[
 \Pr[\text{bad word with exactly one active outer row}]
 \le R(q-1)^{-\rho}\theta^{-\rho e}
 \sum_{w=p+1}^b m_w^\rho g_w.
\]

The new helper evaluates this expression with directed192-bit Arb arithmetic.
The local matrices are exact rational expressions before taking fractional
powers. It sums every coefficient, not a sampled range of row weights.
All receipts below use \(q=65536\), \(N=2^{21}\), one refresh
round per packet, the density allocation for the local uniform component,
and \(e=838860\), corresponding to strict output weight below \(0.4N\).

| Outer | Rows \(R\) | Width/state | \(\theta\) | \(\rho\) | Singleton failure exponent, directed lower |
|---|---:|---:|---:|---:|---:|
| GRS[128,64] | 16384 | 64/16 | 9993/10000 | 3/4 | 135.7782932202 |
| GRS[128,64] | 16384 | 64/8 | 9993/10000 | 3/4 | 135.8310085705 |
| GRS[256,128] | 8192 | 64/4 | 9987/10000 | 3/4 | 222.6995272210 |
| GRS[256,128] | 8192 | 64/8 | 9987/10000 | 3/4 | 222.6866669800 |
| GRS[128,64] | 16384 | 32/8 | 9993/10000 | 3/4 | 143.5027165137 |
| GRS[256,128] | 8192 | 16/4 | 9987/10000 | 3/4 | 223.5258335159 |

The table's decimal truncations are conservative lower bounds. The machine
receipts also provide integer guarantees of135,143,222, or223 bits. These results do
not repair the multi-row small-state failures found by the routing analysis.
In particular, a strong singleton result is compatible with a bad code caused
by partially shortened collections of many rows.

The smaller-packet controls use the same exact map family: \(U=\{0,\ldots,t-1\}\),
\(V=\{128,\ldots,128+s-1\}\), residue interpolation \(A\), and \(C=A^T\).
The degree proof of \(CA=0\) applies because \(2s-2<t-1\); direct finite-field
checks also pass. Their quiet-subcode upper distance is 50%, so that particular
obstruction does not reject 40%. No conclusion about the remaining source
classes follows from this fact or from their singleton receipts.

## Why one common fractional exponent is insufficient

Let \(Z_P(\Pi)=\mathbb E_M[\text{number of bad lines on support }P]\).
The fractional relaxation, even before replacing it by local upper bounds,
contains \(\sum_P\mathbb E_\Pi Z_P(\Pi)^\rho\).

Use the fixed64-row route event from
[HD_BLOCKWISE_GLOBAL.md](HD_BLOCKWISE_GLOBAL.md). Every mixer setup has a
bad subspace of dimension at least2064 on that event. Thus
\(\sum_PZ_P(\Pi)\ge q^{2062}\) there, using the same conservative projective
count bound as the global note. Since
\(\sum_PZ_P^\rho\ge(\sum_PZ_P)^\rho\), at \(\rho=3/4\),

\[
 \log_2\left(\sum_P\mathbb E_\Pi Z_P(\Pi)^\rho\right)
 \ge13880.2090017988.
\]

Therefore this common-exponent bound cannot close the whole code, regardless
of tighter local matrices. This is a lower bound on a proof expression, not
a lower bound on distance-failure probability. The exact conditional cap
remains valid; replacing every cap by the same power is the relaxation that
fails. Additional source-class caps or a route-good event must be retained
before the many-row calculation.

## Inverse-word routing: valid bound and quantified loss

The inner is an involution because \(CA=0\). For a fixed proposed low-weight
output \(Y\), the inverse word \(X=E_M(Y)\) depends on the mixers but not the
route. This is useful, but it does not supply random nonzero amplitudes for
the outer membership test.

For region \(j\), let \(m_j(v)\) count occurrences of field value \(v\) in
the fixed inverse word. A uniform permutation maps that multiset to any one
feasible ordered column with probability

\[
 \gamma_j=\frac{\prod_v m_j(v)!}{R!}.
\]

Fix any \(k\) outer coordinate regions as an information set. Once their
permutations are exposed, the MDS code uniquely prescribes the remaining
columns. Each remaining region either has the wrong multiset, contributing
probability zero, or matches with probability \(\gamma_j\). Independent
regional permutations therefore give an outer-membership bound
\(\prod_{j\notin I}\gamma_j\). The information set may be chosen from the
fixed inverse word before exposing the route.

This type-only bound is too weak for a union over the whole 40% output ball.
Even all-distinct regions give only

\[
 (b-k)\log_2(R!)=13167821.4531505847
\]

bits of outer-membership cost. A directed lower bound on one projective
40%-ball shell exceeds that by more than2290113 bits. This calculation already
divides by \(q-1\) for global scalar lines. Thus any successful inverse proof
must retain multiset consistency or its coupling to the inner-state
randomness, rather than only the displayed type probabilities.

The old support-only MDS membership probability is unavailable here. A
constant column is unchanged by its permutation, and a repeated outer
codeword can have membership probability one. Independent input diagonals
would change that probability space; they are absent from this construction.

## Validation and next step

The new [helper](theory/hd_blockwise_composition.py) and
[seven tests](theory/hd_blockwise_composition_tests.py) check the
noncommuting single-insertion identity, fractional state-path inequality,
unpowered support multiplicities, explicit small route enumeration, type
probabilities, directed singleton receipts, and the common-exponent
obstruction. An independent review found no issue in the cap placement or
the regional chronology.

The next useful step is a many-row source-weight-class cap combined with
the established route exclusions. The exact regional \(Q_h\) formula is
available for that work. Brute-force bivariate and higher-support sweeps
would not address the identified common-exponent obstruction, so they were
not started. No native run, production change, or whole-code claim accompanies
this checkpoint.
