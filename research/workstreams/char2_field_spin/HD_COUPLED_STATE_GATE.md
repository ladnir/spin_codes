# A common-potential gate for the coupled five-cone state

This bounded test compares the proposed coupled recurrence against one
fixed scalar moment law. It does not fit a new weight mark or moment
shape. Its purpose is to determine whether all five boundary-state cones
can retain the successful scalar bound with modest additional cost.

The recurrence and the normalized references are defined in
[the five-cone interface](HD_COUPLED_TWO_STREAM_CONES.md). State order is
$Z,L,R,B,P$. In particular, $P$ represents arbitrary probability mass,
not density relative to a uniform nonzero state. Its potential need not
be less than one.

## Fixed scope and scalar comparison

Fix $r=18$ active groups, eight nonzero RS constituents per group, and
$F_{\max}(G)=1$. Keep each constituent row set fixed through all sixteen
outer columns. The mark is $z=691/700$; all local and regional operators
are unpowered first moments.

Let $R_{m,f}$ be the new five-cone regional cap for a fixed input support
with $m$ nonzero coordinates and $f$ full pairs. The provider averages
the paired placement and fair swaps, raises its packet polynomial to
128, and applies one global 4096-slot normalizer. There is no independent
half normalization, tensor product, or reset of either state.

The frozen scalar target is

$$
M_0=1,\qquad M_m=\eta B^m\quad(m>0),
$$

where $\log_2\eta=-24.529436473994224$ and
$\log_2B=-8.907630419774677$. These are exact binary64 literals from the
pinned earlier F1 zero-state envelope. Using them as a target does not
transfer a receipt from the old construction. The new matrices must
satisfy their own comparison below.

Write a positive potential as $v_i=2^{u_i}$, with $u_Z=0$. Compute

$$
\lambda=\max\left(1,
 \max_{m,f,i}\frac{(R_{m,f}v)_i}{M_m v_i}\right),\qquad
e=\max(0,-\min_i u_i).
$$

The maximum includes every $0\le m\le18$ and
$0\le f\le\min(1,\lfloor m/2\rfloor)$. For every chronological regional
word, positivity gives

$$
e_Z^T\left(\prod_{j=1}^{256}R_{m_j,f_j}\right)\mathbf1
\le 2^e\lambda^{256}\prod_{j=1}^{256}M_{m_j}.
$$

Here $e_Z$ selects the zero initial state and $\mathbf1$ is the terminal
all-ones vector. The terminal cost occurs once. All references, including
$P$, have total mass one, so there is no separate terminal density factor.

## A log-domain initializer, not a numerical certificate

The helper proposes $u$ by a small linear program. For a row with $k$
nonzero entries, impose one constraint per such entry:

$$
\log_2R_{m,f}[i,j]+u_j-u_i-\log_2M_m+\log_2k\le t.
$$

The equal allocation $1/k$ makes these constraints sufficient for the
full weighted row sum to be at most $2^tM_m v_i$. This allocation is
deliberately conservative. Impose $-512\le u_i\le512$, $u_Z=0$,
$t\ge0$, $s\ge0$, and $s\ge-u_i$. Minimize $t+s/256$.

All coefficients are logarithms of moderate magnitude. The LP never
receives tiny potential values such as $2^{-128}$ as constraints. After
solving, the helper recomputes every complete row sum with log-sum-exp
and derives $\lambda$ and $e$ from those sums. It does not trust the LP's
feasibility tolerance, objective, or endpoint variable as a bound. This
remains a FLOAT screen, not directed arithmetic.

Positive $u_P$ is allowed. For example, the three-edge cycle with weights
$Z\to L:1$, $L\to P:2^{-128}$, $P\to Z:2^{-8}$ has balanced log
potentials $(0,-136/3,112/3)$. Restricting all nonzero-state potentials
to values below one would discard this balance.

## Whole-class count and normalization

Use the same raw fixed-row-set weights $W(w)$ and
$\kappa(w)=\min(w,8)$ as in
[the matched two-stream gate](HD_TWO_STREAM_MEAN_GATE.md). Define

$$
Q(X)=\sum_{j=0}^{16}\binom{16}j\eta^j(1-\eta)^{16-j}
 \left(\sum_{w=0}^{j}\binom jw W(w)B^wX^{\kappa(w)}\right)^{18}.
$$

Let $N_1$ be the exact number of source sets in the declared F1 stratum.
The new class first moment is at most

$$
\frac{N_1z^{-734003}}{255\,255^{1152}}\,
2^e\lambda^{256}\sum_{k\ge1296}[X^k]Q(X)^{16}.
$$

The raw $W$ weights already count fixed row sets through the sixteen-column
Hölder bound. There is no additional row-set factor. The whole-word
surrogate floor includes all outer masks, and the regional comparison
includes all five-cone state paths. This covers the declared homogeneous
class only, not all source weights or cohorts.

## Source and test status

`theory/hd_coupled_state_gate.py` is frozen at SHA256
`974db96f33705f1d58f83fece7da5e1bd408e082946115bd9e865565eaa4efe0`.
It pins the packet caps, regional provider, scalar artifact, and reused
positive-polynomial helper. The output retains all 36 regional log
matrices, with null encoding exactly zero entries, not a numeric floor.

Six small tests pass in 0.015 seconds. They cover a tiny-weight LP with
positive $u_P$, full-row recomputation, noncommuting word domination,
complete histogram serialization, the fixed scalar shape, and guards.
No production bank or gate has been run by this helper's author.

A failed gate would reject this fixed scalar-shape comparison, not the
coupled recurrence itself. Any later moment refinement should reuse the
saved bank and identify which row or cone forces the surplus before
requesting another provider or native experiment.

## Frozen-shape result and a different state comparison

Root's bounded run completed in 1.180 seconds. The scalar comparison is
$-88.130205160649$, but the recomputed potential penalty is
$11904.443985543905$, giving log2 bound $+11816.313780383256$.
The saved artifact is `theory/hd_coupled_state_gate_float.json`, SHA256
`0c385a1a2bb3a04d3991cdefb6a312fe59e13add4b2b2afecfcf60a8173b5da0`.

This failure is largely forced by the common scalar shape. At occupancy
18 the live $B\to B$ entry is about $2^{-140.67}$, whereas
$M_{18}$ is about $2^{-184.87}$. A diagonal ratio does not depend on the
potential. Thus no alternative potential can remove that approximately
44.2-bit per-region mismatch while keeping the same scalar target.
This rejects that comparison, not the coupled construction.

The next comparison retains separate zero and live coefficients. Choose
one strictly positive vector $(w_L,w_R,w_B,w_P)$ with sum one, and define
the live family

$$
\nu_w=w_L\nu_L+w_R\nu_R+w_B\nu_B+w_P\nu_P,
$$

where $\nu_L,\nu_R,\nu_B$ are the fixed normalized references from the
five-cone interface, and $\nu_P$ is any probability measure allowed by
cone $P$. The weights remain fixed across every histogram and region.
The member $\nu_P$ may depend on the past, just as it does in the original
five-cone comparison. Positive mixtures of possible $P$ outputs are still
members of that family.

For a five-cone coefficient matrix $T$, define the two-cone cap $S$ by

$$
S_{00}=T_{ZZ},\qquad
S_{01}=\max_{j\in\{L,R,B,P\}}\frac{T_{Zj}}{w_j},
$$
$$
S_{10}=\sum_iw_iT_{iZ},\qquad
S_{11}=\max_j\frac{\sum_iw_iT_{ij}}{w_j},
$$

where both sums range over the four live cones. To verify the direction,
let $E$ be the two-by-five coefficient embedding with rows
$(1,0,0,0,0)$ and $(0,w_L,w_R,w_B,w_P)$. Then
$ET\le SE$ entrywise. Applying this inequality repeatedly preserves
chronological products, even when the matrices do not commute. Since
$E\mathbf1=\mathbf1$, the reduced terminal vector has two ones.
There is no additional endpoint factor.

The saved-bank helper `theory/hd_coupled_state_compress.py` uses the fixed
log weights $(-64,-128,0,-172)$ and normalizes them once. It compresses
each actual $(m,f)$ matrix before taking any entrywise maximum over $f$.
Its source SHA256 is
`3cd4d2de91fa0bb1fa13b267d0157a1bd81f66c86f0f194bf2acd5f66e56d94f`.
Six tiny tests pass in 0.012 seconds, including exact coefficient
domination, noncommuting products, normalized terminal mass, and arbitrary
$P$ mixtures. An independent read-only review found no material issue.
The helper reuses the saved bank; it performs no provider build or bound
evaluation. The following separate-entry gate supplies the numerical
comparison without changing that bank.

## Separate-entry moment result

`theory/hd_coupled_state_moment_gate.py` keeps the normalized live weights
fixed. It compresses each histogram before taking the entrywise maximum
over its allowed full-pair counts. The zero-to-zero entry uses the
prescribed occupancy-5/18 secant; the two cross entries use their own
occupancy-2/18 secants. Each intercept covers every occupancy. The live
diagonal uses its maximum over all occupancies. The zero-state law is
accepted only when its nonzero-atom probability is at most one.

The positive path expansion retains distinct activation and return atoms.
Its cache records both directed transition counts. The resulting column
matrix polynomial is raised to sixteen, retaining every coefficient of
degree at least 1296. There is no common-potential endpoint penalty: the
normalized two-cone reference ends against the all-ones vector.

Root's saved-bank run completed in 2.390 seconds:

| Control | FLOAT log2 result |
|---|---:|
| Matched old shared-cross comparison | +3497.359683472526 |
| Coupled, distinct entry moments | +2906.6054663517416 |

The old result reproduced exactly. The coupled comparison improves it by
590.754217 bits, but still fails the required negative bound by a large
margin. The comparison includes every outer mask and boundary-state path
in the selected homogeneous arity-eight, r18, F1 class. It is not merely
the quiet component. Its mathematical comparison is valid, but these
FLOAT values are not directed receipts or a whole-code certificate.

The source SHA256 is
`d2258e48f19bf3c8b7abea6a3d95024c10ce4d7c7286f49d009f52d8a4d7e89a`;
the seven-test file SHA256 is
`e6c9e375f009779cefc12f96bdf29c081460529983f8571874c8e3ca8df7cf78`.
The result is `theory/hd_coupled_state_moment_gate_float.json`, SHA256
`9f4fc00fdee8d2bd72f2411401a1e8b7897f31d915b1fbac8e388053cac3a5ae`.
Seven tiny tests passed, and independent reviews checked the separate
atoms, noncommuting chronology, fixed-row-set counting, and terminal mass.
No new regional bank or native experiment was used for this comparison.

## Three-packet return: moving the XOR after refresh

The next algebraic control changes the recurrence to
$(U',V')=(V,M(U+c)+V)$, where $c=CX$. This is a different construction,
not a reinterpretation of the saved five-cone matrices. Its local return
event has a useful exact distinction.

Fix nonzero feedbacks $c_1,c_2,c_3$ independently of three fresh invertible
refreshes $M_1,M_2,M_3$, and start at $(0,0)$. Write
$K(x,y)=\Pr[Mx=y]$ for their common refresh kernel. Before-refresh XOR
has states $(0,a)$ and $(a,b)$ after the first two packets, where
$a=M_1c_1$ and $b=M_2(a+c_2)$. Thus its third state is zero exactly when
$c_2=c_3$ and $M_1c_1=c_3$. After-refresh XOR instead has
$b=M_2c_2+a$, so its third state is zero exactly when
$M_1c_1=M_2c_2=c_3$. Consequently,

$$
\Pr[\text{before returns at 3}]
=1_{c_2=c_3}K(c_1,c_3),\qquad
\Pr[\text{after returns at 3}]
=K(c_1,c_3)K(c_2,c_3).
$$

If $D=q^h-1$ and every nonzero refresh atom is at most $\beta/D$, these
are bounded by $\beta/D$ and $\beta^2/D^2$, respectively. This statement
conditions on fixed feedbacks and the zero initial state. It is not a
two-hit assertion from an arbitrary incoming state or a complete marked
packet bound.

`theory/hd_coupled_state_return.py` checks this with the actual toy law
$D_2HD_1$ over GF(4), state dimension two, and
$H=\left(\begin{smallmatrix}1&1\\1&\theta\end{smallmatrix}\right)$,
where $\theta^2+\theta+1=0$. It retains all 81 equally likely nonzero
diagonal-key tuples, including duplicate maps. Exact state counting
covers all 3375 fixed nonzero feedback triples and all three fresh keys.
The worst return probabilities are exactly $1/9$ before refresh and
$1/81$ after refresh, attaining both stated bounds with $\beta=5/3$.
Five tiny tests pass in 0.762 seconds.

The distinction is conditional concentration, not a universal average
improvement. If all three feedbacks are independent uniform nonzero
vectors, both return probabilities average to $1/D^2$; the toy check
gives $1/225$ for each. A useful next proof must retain the actual
feedback law, emission weights, and all persistent boundary states.
