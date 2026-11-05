# Whole-word moment fibers of the G16 outer

The failed regional envelopes do not establish that the G16 construction
has insufficient distance. They discard dependencies between a fixed RS
row's sixteen symbols. This note gives a different, whole-word interface.
It proves a small exceptional-setup bound, not a new distance certificate.
No native implementation or frozen distance receipt changes.

## The kernel is a genuine subcode

Work over $F=\operatorname{GF}(256)$, with $q=256$ and
$C=\operatorname{RS}[16,8]\subseteq F^{16}$.
Fix a set $A$ of $a$ distinct affine mixer labels $u_i$.
An ambient message in this row set is $c=(c_i)_{i\in A}\in C^A$;
this space includes zero constituent rows.
Let $d_{ij}$ be the independent nonzero input-diagonal coefficients.
Define the two-moment map

$$
S_A(D):C^A\longrightarrow F^{32},\qquad
c\longmapsto
\left(\sum_{i\in A}d_{ij}c_{ij},
      \sum_{i\in A}u_i d_{ij}c_{ij}\right)_{j=1}^{16}.
$$

The affine two-sum mixer is $H=I+\mathbf1\mathbf1^T+uu^T$.
If $S_A(D)c=0$, its input and output agree in every column, before the
output diagonal. Consequently this subcode remains supported on rows $A$.
For every setup,

$$
\dim\ker S_A(D)\geq8a-32.
$$

Thus arities eight and nine have kernel dimension at least 32 and 40.
This is not a rare cancellation that a proof may discard.
It is also not a low-distance counterexample to the complete routed inner
code. A lower bound on kernel dimension is not a weight-enumerator upper
bound.

## A simultaneous rank guarantee

Fix $A$ before sampling $D$.
A nonzero left-null candidate is a vector
$\lambda=(\lambda_{0j},\lambda_{1j})_{j=1}^{16}\in F^{32}\setminus\{0\}$.
For each row, put

$$
v_i(j)=\lambda_{0j}+u_i\lambda_{1j}.
$$

The candidate annihilates $S_A(D)$ exactly when
$\operatorname{diag}(d_i)v_i\in C^\perp$ for every $i\in A$.
Distinct labels imply that at most one $v_i$ is identically zero.
If a nonzero $v_i$ has support size $w$, its diagonal randomization is
uniform on the nonzero values of that fixed support.
The shortened dual RS code has dimension $w-8$ when $w\geq9$.
An information set therefore bounds its full-support words by
$255^{w-8}$. The membership probability is at most $255^{-8}$;
it is zero for $1\leq w\leq8$.
These membership events use independent diagonals in different rows.

There are $a(q^{16}-1)$ nonzero candidates with one identically zero row.
The remaining $q^{32}-1-a(q^{16}-1)$ have no such row.
Scalar multiples describe the same left-null event, so a projective union
bound gives

$$
\Pr_D[\operatorname{rank}S_A(D)<32]
\leq b_a:={
 (q^{32}-1-a(q^{16}-1))255^{-8a}
 +a(q^{16}-1)255^{-8(a-1)}\over255}.
$$

Let $E$ be the event that every one of the 8192 source groups has full
syndrome rank for every row set of size at least five. No independence
between overlapping row sets is needed:

$$
\Pr[\neg E]\leq
8192\sum_{a=5}^{16}\binom{16}{a}b_a=:\varepsilon_{\rm rank}.
$$

The new exact-rational helper `theory/hd_two_sum_syndrome_rank.py` evaluates
this expression. It lies strictly between $2^{-47}$ and $2^{-46}$.
Its logarithm at 256 and 384 working bits is enclosed around

$$
\log_2\varepsilon_{\rm rank}
=-46.6757337702933229908836468566516754416544641825\ldots.
$$

The respective reported radii are below $4.26\cdot10^{-75}$ and
$7.32\cdot10^{-114}$. These values concern only the rank event.
The initial six exact small tests passed in 0.027 seconds. An additional
GF8 test now checks the whole-row marginals and row-Hölder inequality.
The seven rank tests and four local-count tests pass together in
0.038 seconds. They also check GF4 projective classes, shortened-dual
atoms, equal-sized fibers, and gauge invariance.

The probability argument uses ideal independent diagonals before native
normalization. Multiplying every input coefficient in one column by the
same nonzero scalar rescales that column's two syndrome coordinates.
It preserves rank. This proves transfer through the existing column gauge;
it does not assert independence of the normalized coefficients.

## What a good setup supplies

Fix a realization in $E$, a row set $A$ of size $a\geq5$, and a syndrome
$\sigma\in F^{32}$. Its fiber has exactly

$$
|S_A(D)^{-1}(\sigma)|=q^{8a-32}
$$

ambient messages. Uniform ambient messages consequently give independent
uniform moment pairs across all sixteen columns.
After conditioning on a syndrome, the remaining row supports need not be
independent between columns.

There is a stronger row-marginal fact. For $a\geq6$, put $k=a-5$.
Fix any $k$ rows and prescribe their complete RS words arbitrarily.
The remaining five-row syndrome map is surjective, so exactly $q^8$
completions produce the requested syndrome.
Hence every $k$-row marginal of the uniform fiber is uniform on $C^k$.
This statement includes zero rows. Conditioning all rows to be nonzero
would change those marginals.

In particular, every complete row is uniform on the $q^8$ RS words.
A union bound over zero rows therefore gives

$$
\#\{c\in S_A(D)^{-1}(\sigma):c_i\ne0\text{ for every }i\in A\}
\geq q^{8a-32}(1-a/q^8),\qquad a\geq6.
$$

Thus almost all words of each arity-eight or arity-nine fiber have exactly
that active row set; the nuisance kernel is not explained by smaller-row
subcodes. This is an **unweighted** count. It does not justify dropping
nonzero-row indicators from a weighted sum or claiming independent rows
after conditioning them all to be nonzero.

For nonnegative row functions $f_i:C\to\mathbb R_{\geq0}$, generalized
Hölder over all $k$-row subsets now gives the useful bound

$$
\sum_{c\in S_A(D)^{-1}(\sigma)}\prod_{i\in A}f_i(c_i)
\leq q^{8a-32}
       \prod_{i\in A}\left(q^{-8}\sum_{c\in C}
                    f_i(c)^{a/k}\right)^{k/a}.
$$

To see the exponents, each row belongs to $\binom{a-1}{k-1}$ of the
$\binom ak$ subsets. Each Hölder factor therefore raises a row function
to $a/k$. The known $k$-row marginal then factors its expectation.
For exact-active row sets, put the nonzero-row indicator inside each
$f_i$. The sum on the right is then over $C\setminus\{0\}$.
This retains a whole RS word, including its nine-coordinate support floor,
rather than selecting a new row set or support independently per column.

For the zero syndrome, the inside output is $d_i c_i$.
Nonzero diagonal scaling preserves its support. If
$f_i(c)=t_i^{\operatorname{wt}(c)}\mathbf1_{\{c\ne0\}}$, the right-hand
side uses the exact ordinary RS weight enumerator at $t_i^{a/k}$.
For other syndromes the inside row has a fixed affine offset; its marked
coset sum must be bounded separately. It is not the same enumerator.

## Where the fractional exponent may enter

Fix the active source set $G$, its exact row sets $A_g$, all input
diagonals, and the complete route. Let $a_g\geq5$ and
$d=\sum_{g\in G}(8a_g-32)$.
Let $B_{G,\mathbf A}$ be the event that a message with exactly these
active groups and row sets produces output weight at most $L=734003$.
Fix a weight marker $0<z\leq1$.
For a syndrome array $\boldsymbol\sigma$, let
$\Psi_{\boldsymbol\sigma}$ be the mean of
$z^{\operatorname{wt}(Y)}\prod_{g,i}\mathbf1_{\{c_{g,i}\ne0\}}$
under the product of uniform ambient fibers, the independent output
diagonals, and the inner refresh keys.
Thus the row-set indicators are retained in the expectation, rather than
silently conditioning the fiber distribution.

The conditional first moment and
$\min(1,x)\leq x^\alpha$, followed by subadditivity over syndromes, give

$$
\Pr[B_{G,\mathbf A}\mid D,\Pi]
\leq\left({z^{-L}\over255}\right)^\alpha q^{\alpha d}
       \sum_{\boldsymbol\sigma}
                \Psi_{\boldsymbol\sigma}^{\alpha},
\qquad 0<\alpha\leq1.
$$

This formula is used only for realizations in $E$.
Average its right-hand side over those input diagonals and routes, then
union over the fixed source and row sets. Charge $\Pr[\neg E]$ once.
The factor $255$ removes projectively equivalent nonzero messages.
Source-set and row-set unions remain outside $\alpha$.

The fiber cardinality is therefore legitimately charged under $\alpha$.
However,

$$
(\mathbb E_{\rm fiber}\Phi)^\alpha
\not\leq\mathbb E_{\rm fiber}[\Phi^\alpha]
$$

in general; Jensen gives the opposite direction. Feeding fiber samples
directly into the old packetwise fractional provider would make precisely
this invalid step. A sound implementation must first bound the unpowered
fiber expectation, for example by a positive sum of the row-factor bounds
above. It may then raise that complete bound to $\alpha$.

## Cheapest next obstruction check

The full-rank fact alone is not the missing gain. The previous uniform
column comparison already assigns about $q^{8a-32}$ mass to the all-zero
moment branch. Replacing that mass by its exact fiber cardinality saves
only its small diagonal-density inflation. Any useful gain must come from
retaining the whole-row support correlations inside the fiber.

Before implementing a general fiber operator, test one positive term:
all groups have arity eight or nine, all syndromes vanish, and the regional
boundary path remains in state zero. Use an **unpowered** regional bound
$R_1(m,f)_{00}\leq A_{00}B_{00}^m$ covering all feasible pairs.
Write

$$
P_C(t)=\sum_{c\in C\setminus\{0\}}t^{\operatorname{wt}(c)},\qquad k=a-5.
$$

The corresponding term in the proposed whole-word majorant is

$$
\binom{8192}{r}\binom{16}{a}^{r}
\left[
 {z^{-L}\over255}A_{00}^{256}
 \left\{q^{8a-32}
       \left(q^{-8}P_C(B_{00}^{a/k})\right)^k\right\}^{r}
\right]^\alpha.
$$

This is one term of an upper-bound construction, not a lower bound on
actual failure probability. If it already exceeds the target, that proposed
majorant cannot close without another refinement. If it passes, the other
syndromes and state paths still require analysis.
The computation needs one degree-16 RS polynomial, not an operator grid.
No unpowered-provider control has run yet. A saved fractional regional bank cannot be
silently converted to $R_1$ by taking its $1/\alpha$ power.

This next check deliberately gives up route fractionalization until a
valid fiber-aware local interface is available. It isolates whether exact
RS row correlations offer enough improvement to justify that larger step.

### Six local counting comparisons

The isolated helper `theory/hd_two_sum_fiber_activity.py` ran six
counting comparisons in 0.076 seconds, using arities eight and nine and
the three previously frozen 0-to-0 activity slopes. It rebuilt no regional
provider and fitted no parameter. All integer support counts are exact;
the subsequent positive polynomial evaluations are FLOAT, not certificates.

Let $b$ be a saved log2 activity marker and retain $\alpha=1/5$.
The fiber expression is evaluated at $t=2^{b/\alpha}$ and then raised
to $\alpha$. This keeps the final activity monomial $2^{b\operatorname{wt}}$.
It is a counting comparison, not a way to recover an unpowered regional
envelope from a fractional one.

Two old orderings provide controls. Both keep one fixed labelled row set,
exclude rows that never appear, and require total input weight at least
$9a$. They still relax each row's individual nine-column floor.
The mean-first control sums the sixteen-column quiet law before taking
$\alpha$. The column-first control powers the exact quiet-mask counts
inside each column instead. A positive seen-row recurrence performs both
sums without signed FLOAT cancellation.

| Frozen anchor | Arity | New fiber then alpha | Old mean then alpha | Old column alpha |
|---:|---:|---:|---:|---:|
| 18 | 8 | -33.436412 | -51.693390 | +30.066918 |
| 18 | 9 | -36.261598 | -51.755142 | +40.678815 |
| 24 | 8 | -40.187456 | -62.732239 | +22.449129 |
| 24 | 9 | -43.997203 | -64.174635 | +32.184757 |
| 32 | 8 | -45.906847 | -71.879915 | +16.075733 |
| 32 | 9 | -50.490501 | -74.471154 | +25.081406 |

These are local log2 weighted counts, without source counts, inner path
intercepts, or the global forbidden-weight prefix.
Postponing mask fractionalization recovers about 62--77 bits per group
relative to the column-first comparison. The new fiber Hölder bound is
nevertheless 15.5--26.0 bits worse than the old mean-first control.
The measured opportunity is the order of fractionalization, not a rank-only
entropy saving. The new expression does retain every row as an actual
nonzero RS word; no wholly empty row creates its apparent gain.

The artifact is `theory/hd_two_sum_fiber_activity_anchors_float.json`,
SHA-256 `410d975a6becba7a20717ace183717863db7cbe5a6a94f23ff715a392e4930e7`.
Its source hash is
`489d1ce418f04467934ac8442cdabfbb6361ee5d5ba3a954c2e08ecdd308a649`.

## A column transfer before mask fractionalization

The following alternative does not condition on the rank-good event.
It returns to the expected outer counting measure and retains each whole
column's mask sum before applying $\alpha$.

Fix $G$ and the same exact row sets $A_g$ across all sixteen columns.
Write $a_* = \sum_g |A_g|$.
In one column, a comparison input is $x=(x_g)_{g\in G}$ with
$x_g\in F^{A_g}$. These inputs may contain zeros.
Put $k(x)=\sum_g\operatorname{wt}(x_g)$ and
$h(x)=\#\{g:x_g\ne0\}$.
Let $\Pi_j$ be the sixteen regional permutations and swap bits belonging
to column $j$. All these route blocks are independent and identically
distributed: the same sixteen-region alternating matching pattern repeats
in each column. The fixed row sets are not resampled between them.

For one fixed route block and comparison input, let
$\mathcal R_z(\Pi_j,Hx)$ be the existing unpowered two-state upper matrix
through its sixteen regions, in their original chronological order.
Here the route is still fixed, so the matrix uses the corresponding
packet-placement caps, not an already route-averaged regional table.
The input amplitudes are averaged over the independent output diagonal,
and the inner refresh keys are already averaged in the local caps.
The two coordinates bound zero-state mass and the established nonzero-state
density. Positive mixtures remain valid upper operators on this cone.
In particular, $\mathcal R_z$ is not claimed to bound every fixed refresh
realization and is not a stochastic matrix.

Define the positive matrix polynomial

$$
T_{G,\mathbf A}(\Pi_j;X,Y)=255^{-a_*/2}
 \sum_{x\in\prod_g F^{A_g}}
 X^{k(x)}Y^{h(x)}\mathcal R_z(\Pi_j,Hx).
$$

The scalar normalizer distributes the prescribed-support comparison
$255^{-8a_*}$ across sixteen columns. Exact joint input/output support
counts may implement the sum; field vectors need not be enumerated.
Crucially, all masks contributing to one coefficient are summed before
that coefficient is powered.

Let $\mathcal L$ sum coefficients with total $X$ degree at least $9a_*$
and total $Y$ degree at least $9r$, where $r=|G|$.
These are necessary aggregate constraints for the actual RS messages.
They do not assert each individual row floor or exclude every synthetic
empty group in the relaxed comparison.
At the complete fixed route, the conditional first moment is at most

$$
{z^{-L}\over255}\,
\mathcal L\left[e_0^T
 T(\Pi_1;X,Y)\cdots T(\Pi_{16};X,Y)\mathbf1\right].
$$

Expand only its column-degree allocations and boundary-state paths before
fractionalization. With entrywise powers, define

$$
U_{k,h}:=\mathbb E_{\Pi_j}
 \left([X^kY^h]T(\Pi_j;X,Y)\right)^{\circ\alpha},
\qquad U(X,Y)=\sum_{k,h}U_{k,h}X^kY^h.
$$

Independence of the route blocks gives the valid class bound

$$
\Pr[B_{G,\mathbf A}]
\leq\left({z^{-L}\over255}\right)^\alpha
 \mathcal L\left[e_0^T U(X,Y)^{16}\mathbf1\right].
$$

The fixed $\mathbf A$ remains an argument of every coefficient of $U$.
Any subsequent row-set Hölder bound must use these new coefficients.
The old per-mask arity moments cannot be substituted unchanged.
Source and row-set unions remain outside $\alpha$, and the global
projective factor occurs once.

### Computation and a common route partition

The coarsest safe shortcut is Jensen in its upper direction:

$$
U_{k,h}\leq
 \left(\mathbb E_{\Pi_j}[X^kY^h]T\right)^{\circ\alpha}.
$$

The mean on the right may use unpowered regional operators, since those
sixteen routes are independent. This shortcut retains the mask sum but
loses route fractionalization. The exact local counting controls above
do not establish that its complete bound will pass.

Each coefficient has four scalar entries, with $0\leq k\leq a_*$ and
$0\leq h\leq r$. The final two-degree positive matrix power is finite.
The hard part is evaluating the column coefficients: without a further
reduction, sixteen regional occupancy coordinates give $(r+1)^{16}$
possibilities, and source-pair geometry may require labelled subsets.
An unpowered entrywise exponential envelope permits a cheap path expansion,
but the complete column mask sum must still precede the power.

A sharper shortcut conditions on a finite event $E_\ell$ of the route:

$$
U_{k,h}\leq\sum_\ell\Pr[E_\ell]
 \left(\mathbb E[[X^kY^h]T\mid E_\ell]\right)^{\circ\alpha}.
$$

Every event must concern the placement of the entire fixed source set $G$
and must be common to all masks in $T$.
The simplest proposed event records collisions among those potential
positions in the four 16-coordinate emission blocks of each packet.
Call a placement good when each block contains at most four positions
from $G$. Every subset chosen by an outer mask then has at most four
positions per block. Any four columns of the RS$[16,4]$ feedback block
are independent. Therefore a nonempty forced-nonzero packet input cannot
map state zero back to state zero on a good placement.
This is only a **per-packet** assertion. It does not make the regional
0-to-0 entry zero: one packet can activate the state and a later packet
can return it to zero. In particular, it does not by itself remove the
observed two-input regional activation/return contribution.

The complementary placements use the existing unpowered packet cap.
Their probabilities stay outside $\alpha$ in conditional Jensen; a rare
collision probability $p$ is not automatically weakened to $p^\alpha$.
A histogram that forgets group labels requires an exact conditional
average or a universal bound over its labelled placements. It does not
authorize independent relabelling for each mask or between regions.

The next gate should compare this common-source good/bad partition against
the unconditional Jensen column bound for one frozen class. It must first
derive the corresponding conditional mean operators. No conditional
provider, column grid, Monte Carlo certificate, or new distance receipt
has been computed in this note.

One optional refinement concerns two consecutive packets, not a general
region. Each input coordinate affects at most four feedback coordinates.
If the first packet starts at zero and activates the state, with input
weight $m_1$, its pre-refresh state has weight at most $4m_1$.
The current MDS refresh has branch number 17, so its output state has
weight at least $17-4m_1$. A following input of weight $m_2$ cannot cancel
that state when $4(m_1+m_2)<17$.
This only uses an upper bound on a feedback column's support; some columns
have weight one, and none are assumed to have exactly four nonzeros.
An intervening empty-packet refresh changes the state and invalidates this
two-consecutive-packet deduction. A useful state-weight refinement would
therefore need an additional proof interface. It is not part of the main
column-before-alpha proposal or of any recorded control.
