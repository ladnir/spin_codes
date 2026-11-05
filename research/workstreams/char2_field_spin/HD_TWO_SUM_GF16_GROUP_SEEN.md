# Positive group-activity constraints: a bounded negative control

This control adds two positive constraints to the existing fractional,
fixed-type state-path bound. One excludes every wholly empty source group.
The other retains the necessary aggregate floor of $9r$ nonempty
group-column incidences. The constraints are evaluated independently.
Their minimum remains an upper bound, but it does not impose their intersection.

Six saved-bank FLOAT controls completed in 0.651 seconds. The aggregate
floor improves the $r=18$, $F=1$ estimate by 307.952216 bits, but every
control still fails. Neither constraint improves the retained estimates
at $r=24$ or $r=32$. No regional operator, setup sample, or witness was rebuilt
or fitted; no directed certificate is claimed.

## Fixed objects and the two necessary constraints

The construction is the existing g16 affine two-sum outer, alternating
independently permuted and fairly swapped pairs, and h4 inner. Fix an indexed
source set $G$ of size $r$. Its stratum is
$F=\max(a(G),b(G))=1$, where $a(G)$ and $b(G)$ count full source pairs in
the two adjacent matchings. Every regional subset $U\subseteq G$ obeys
$f(U)\le\min(F,\lfloor|U|/2\rfloor)$. The saved regional envelopes use only
this upper constraint. The exact indexed-source count $N(r,F=1)$ is charged
once, outside the fractional exponent.

Every group in $G$ contains a nonzero RS$[16,8]$ row. That row is nonzero
in at least nine of the sixteen outer columns. The input diagonals preserve
these supports, and each column mixer and output diagonal are invertible.
Consequently each group has a nonempty output mask in at least nine columns.
For its activity bits $e_{g,j}\in\{0,1\}$, this implies both

$$
\forall g:\quad\sum_{j=1}^{16}e_{g,j}\ge1,
\qquad
\sum_{g=1}^{r}\sum_{j=1}^{16}e_{g,j}\ge9r.
$$

The two implemented bounds impose these respective consequences, not the
individual nine-column floor.

## Activity patterns remain fixed through Hölder

For a chronological local state path $p$, let $f_{a,e}(p)$ be its
nonnegative column contribution for fixed group type $a\in\{1,\ldots,16\}$
and activity $e$. The fixed-row-set moments, regional mask markers, and
optional support-floor witnesses are already included. The empty branch
uses only the zero mask; the nonempty branch sums the labelled positive-weight
masks. Type weights remain fixed across paths and outer columns.

Define

$$U_e(p):=\left(\sum_{a=1}^{16} f_{a,e}(p)^{16}\right)^{1/16}.$$

Fix the entire labelled activity-bit array before applying the existing
sixteen-column Hölder inequality over the group-type vector. The resulting
bound is a product of the corresponding $U_{e_{g,j}}(p_j)$. Only then sum
over allowed activity arrays. This chronology does not resample group types
or constituent row sets. Splitting activity branches can introduce slack:
$U_0(p)+U_1(p)$ may exceed the old norm of $f_{a,0}(p)+f_{a,1}(p)$.

Let $c_p$ be the positive path coefficient, including its multiplicity and
any existing empty-region cut coefficient. For $0\le\ell\le r$, define the
ordered two-state matrix

$$
(M_\ell)_{ij}:=
\sum_{p:i\to j}c_p\,U_0(p)^{r-\ell}U_1(p)^\ell.
$$

No factor $\binom r\ell$ is included in $M_\ell$.

## Two positive evaluations

For exact nonzero-group exclusion, let $c$ count groups already seen in
previous columns. A transition that first sees $t$ new groups has block

$$
T_{c,c+t}=\binom{r-c}{t}\sum_{j=0}^{c}\binom cj M_{t+j}.
$$

The binomial factors select newly active unseen groups and currently active
previously seen groups. Compose sixteen column transitions in order. Start
at inner state zero and $c=0$; sum both terminal inner states only at $c=r$.
There are $2(r+1)$ states, at most 66 in this control.

For the aggregate floor, form the positive matrix polynomial

$$K(Y)=\sum_{\ell=0}^{r}\binom r\ell M_\ell Y^\ell.$$

With $e_0^T=(1,0)$ and $\mathbf1=(1,1)^T$, retain coefficients of degree at
least $9r$ in $e_0^T K(Y)^{16}\mathbf1$. The maximum degree is $16r\le512$.
Both evaluations preserve chronological matrix multiplication and use
positive sums only.

An exact per-group floor cannot simply move inside a columnwise Hölder
norm. For example, with $f_0(X)=1$ and $f_1(X)=X$,
$[X^{\ge1}](f_0f_1)=1$, whereas the geometric mean of the same tails in
$f_0^2$ and $f_1^2$ is zero. The activity-history and aggregate-polynomial
evaluations avoid that invalid interchange.

## Six frozen results

The source bank is `hd_two_sum_gf16_source_strata_anchors_float.json`, SHA-256
`5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf`.
Only its $F=1$ stratum is used. The witnesses retain $\alpha=1/5$,
$z=691/700,172/175,171/175$ at the three respective source weights, all
numeric slopes, and the two stored type-marker templates. Every old baseline
was reproduced before its new evaluation.

| $r$ | Stored template | Old log bound | Seen-all log bound | Aggregate-$9r$ log bound |
|---:|---|---:|---:|---:|
| 18 | old-all-zero | 2757.796905 | 2732.763375 | 2449.844689 |
| 18 | joint-empty-floor1 | 3787.468689 | 3787.468891 | 3787.468891 |
| 24 | old-all-zero | 3785.532458 | 3785.953415 | 3785.953415 |
| 24 | joint-empty-floor1 | 5448.493390 | 5448.493503 | 5448.493503 |
| 32 | old-all-zero | 6035.768213 | 6036.455559 | 6036.455559 |
| 32 | joint-empty-floor1 | 8041.831804 | 8041.831812 | 8041.831812 |

The old estimate is retained whenever it is smaller. These are positive
logarithmic majorants, not successful failure-probability certificates.
The result covers neither all source strata nor any low-total-arity restriction.

The helper and driver pass 14 exact or mocked tests in 0.021 seconds. Tests
enumerate labelled activity words with noncommuting endpoint matrices,
fixed-type Hölder examples, exact zero branches, binomial multiplicities,
baseline identity failures, source scope, and bounded partial results.

Provenance:

- `theory/hd_two_sum_gf16_group_seen.py`:
  `9f049719e50f04a663ef3c4415fd34e43f81e28cd062312aac847d87dcb0e733`
- `theory/hd_two_sum_gf16_group_seen_control.py`:
  `36701fa670e8c1606e33be0e0a4648089b458b77132a92ab2ee8c744dccc95bb`
- `theory/hd_two_sum_gf16_group_seen_f1_anchors_float.json`:
  `0e1b81676a8c52de24f61fee0da141a03cc2e6bf500a21eff96e0d143348d074`

The current evidence supports retaining the positive interfaces, but not
expanding to a large age-histogram DP without a new diagnostic. In particular,
removing empty-group patterns helps the first anchor but does not explain
the larger-anchor gap. The next proof step should isolate the remaining
live-state or fixed-type relaxation, while preserving the separate directed
high-total-arity receipts. No construction change, performance claim, or
whole-code certificate follows from this control.
