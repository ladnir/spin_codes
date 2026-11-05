# Post-refresh coupling: two-hit returns and a scoped mean gate

This note studies a proposed inner recurrence. It changes neither the outer
code nor the paired-route law. The first saved-bank FLOAT gate improves the
preceding coupling control but still fails: its log2 upper bound is
`+2282.476801465376` for the selected homogeneous class. This is not a distance
certificate, a native implementation, or a timing result.

## Recurrence and setup randomness

Work over $\mathbb F_q$, with $q=256$. Let $A:\mathbb F_q^{16}\to
\mathbb F_q^{64}$ be the established four-block RS emission map, and let
$C:\mathbb F_q^{64}\to\mathbb F_q^{16}$ satisfy $CA=0$.
The state is $(U,V)\in(\mathbb F_q^{16})^2$, initially zero. At packet position
$t$, emit and update

\[
Y_t=X_t+AU_{t-1},\qquad
(U_t,V_t)=\bigl(V_{t-1},M_t(U_{t-1}+CX_t)+V_{t-1}\bigr).
\]

All $M_t$ are sampled independently **once in setup** and fixed for every
message and encoding. In a local proof, a fresh refresh means an unexposed
independent setup key at a new packet position. It does not mean per-call
randomness. Each $M_t$ is invertible and obeys the established atom bound

\[
\Pr[M_t a=b]\le\varepsilon:=255^{-16}
\quad\text{for fixed }a,b\ne0.
\]

The probability is over that setup key. Set $d=q^{16}-1$ and
$\beta=d\varepsilon$. The reference $\mu$ is uniform on the $d$ nonzero state
vectors. Marked measures use $z^{\operatorname{wt}(Y)}$, with $0<z\le1$.

The preceding coupling used $M_t(U+V+CX)$ instead. Moving the extra addition
after $M_t$ changes the sampled linear map; old distance receipts do not
transfer merely because the arithmetic count is similar.

## Exact two-packet identity

Fix the old state $(u,v)$ and two input packets independently of their two
refresh keys. Write $c_i=CX_i$ and

\[
d_1=M_1(u+c_1),\qquad d_2=M_2(v+c_2).
\]

After two packets the state is

\[
(U_2,V_2)=(v+d_1,\ v+d_1+d_2).
\]

The two emitted packets are $X_1+Au$ and $X_2+Av$. Their product mark is
therefore independent of both refresh keys under this conditioning. The
innovations $d_1,d_2$ are independent, but either may be identically zero.
If both refresh inputs are nonzero, every joint output atom is at most
$\varepsilon^2$. The two final coordinates then differ, so neither a zero
state nor a diagonal state is possible.

The zero-innovation cases cannot be discarded. For example, if $v=c_2\ne0$,
then $d_2=0$. Returning to zero can require only the one hit $d_1=v$.
Thus there is no universal two-hit return bound for arbitrary old states.

## First return from zero

Fix a complete input word independently of all inner refresh keys. Remove
leading packets with zero feedback, and assume $c_1=CX_1\ne0$.
With $a=M_1c_1$ and $b=M_2c_2$, the state after two packets is $(a,a+b)$.
Return at packet three occurs exactly when

\[
M_1c_1=M_2c_2=c_3\ne0.
\]

Its probability is at most $\varepsilon^2$. For the preceding coupling,
the corresponding condition is only $M_1c_1=c_2=c_3\ne0$.

There is also a bounded first-return statement for a longer fixed word.
Define $V_{-1}=V_0=0$. The recurrence becomes

\[
V_t=V_{t-1}+M_t(V_{t-2}+c_t),\qquad U_t=V_{t-1}.
\]

For every fixed $h\ne0$, induction gives
$\Pr[V_t=h]\le t\varepsilon$. If the refresh input is zero, the previous
value persists. Otherwise, conditioning on the past leaves a fixed-target
refresh hit of probability at most $\varepsilon$.

Let $\tau$ be the first return to $(0,0)$ after activation. If $\tau=n\ge3$,
then $V_{n-1}=0$ and $V_{n-2}=c_n=:h\ne0$. The condition $h\ne0$ follows
because $h=0$ would give a return one packet earlier. Packet $n-1$ must hit

\[
M_{n-1}(V_{n-3}+c_{n-1})=h.
\]

Conditioning through packet $n-2$ and then using the preceding marginal bound
proves

\[
\Pr[\tau=n]\le(n-2)\varepsilon^2,\qquad
\Pr[\tau\le T]\le\frac{(T-2)(T-1)}2\varepsilon^2.
\]

These are unweighted fixed-word bounds, not a union over codewords. Marks
at most one may be discarded to apply them. They cannot be multiplied by
an independently averaged emission mark without a further conditioning
argument. In particular, diagonal point states can persist when feedback
cancels their refresh input.

## Eight positive comparison families

The packet operator uses the following normalized references or families.
They dominate weighted histories; they do not assert actual state independence.

| Name | Reference or family on $(U,V)$ |
|---|---|
| $Z$ | Point mass at $(0,0)$ |
| $L$ | $\delta_0\otimes\mu$ |
| $R$ | $\mu\otimes\delta_0$ |
| $B$ | $\mu\otimes\mu$ |
| $D$ | $(u,u)$ with $u\gets\mu$ |
| $P$ | An arbitrary probability measure on nonzero diagonal points |
| $F$ | Arbitrary nonzero first coordinate, with conditional second-coordinate density dominated by $\mu$ |
| $G$ | An arbitrary nonzero first coordinate and zero second coordinate |

An $F$ coefficient bounds the mass of a dominating conditional-flat family,
not merely the mass of an unrestricted pair distribution. Positive sums of
$P$, $F$, or $G$ contributions remain in their respective families.

The implementation retains two important distinctions. A return from $D$
can create a diagonal point spike, so it goes to $P$, not $D$. A nonreturn
from $P$ reaches $G$ only after a fixed-target refresh hit, paying $\beta/d$.
The following $G\to Z$ transition may be cheap; its incoming coefficient
already contains this extra factor.

For $B$, the marked variable $W=M(U+CX)$ is independent of the reference $V$.
Consequently the old flat-state cap $T_{11}$ gives $B\to B$ coefficient
$T_{11}$ and $B\to R$ coefficient $T_{11}/d$. This reuse is valid after
refresh; it was not the corresponding argument for the previous coupling.
The source fixes one $D$-row decomposition and does not entrywise-min it
with a different decomposition into incompatible families.

## One live reference, chosen before the bank

In the live order $(L,R,B,D,P,F,G)$, fix unnormalized weights

\[
\widetilde w=(1,2^{-128},1,1,2^{-128},2^{-128},2^{-256}),
\qquad w=\widetilde w/\sum_i\widetilde w_i.
\]

This is a prescribed trial reference, not an optimized law. Its powers of
two do not replace the exact factors $1/d$ in packet caps. The rule leaves
common births into $L$ and $D$ unpenalized before normalization. It gives
$R$ and $P$ one inverse-state-scale weight, and $G$ two. The $F$ weight
inherits its $P$ birth scale. This ancestry motivates the trial; only the
explicit matrix domination below establishes its validity.

Let $T$ be any eight-state regional coefficient matrix. With the same $w$
for every region and histogram, define

\[
S_{00}=T_{ZZ},\quad S_{01}=\max_j T_{Zj}/w_j,\quad
S_{10}=\sum_iw_iT_{iZ},\quad
S_{11}=\max_j\frac{\sum_iw_iT_{ij}}{w_j},
\]

where $i,j$ range over live families. For the coefficient embedding
$E=(e_Z;(0,w))$, these definitions give $ET\le SE$ entrywise. Positive
matrix multiplication preserves this inequality in the original order.
Since $E\mathbf1_8=\mathbf1_2$, the terminal vector is exactly
$\mathbf1_2$. There is no repeated normalization or extra endpoint factor.
The family representatives may change between steps, but the weights cannot.

## Saved-bank gate and outcome

The single raw bank uses $z=691/700$, all occupancies $0\le m\le18$, and
full-pair counts $0\le f\le\min(1,\lfloor m/2\rfloor)$. It contains 36
eight-state matrices. The paired-route normalization is applied once after
the ordered packet power. No source geometry is resampled across regions.

After compression, the gate takes the entrywise maximum over $f$. It uses
the prescribed $00$ secant $(5,18)$, separate $01$ and $10$ secants $(2,18)$,
and a constant $11$ bound. Every intercept covers all $m$, including zero.
The $00$ law must be a valid two-atom probability law; failure is reported
as incomplete rather than clipped. Off-diagonal zeros at $m=0$ are relaxed.

The positive matrix-polynomial engine includes every mask and boundary path.
It preserves each outer row set across all 16 columns through the existing
fixed-row-set Holder comparison. The surrogate input degree is
$\kappa(w)=\min(w,8)$, and the complete word retains degrees at least 1296.
Raw outer weights omit the RS denominator, so the final expression is

\[
\frac{N_1z^{-734003}}{255\,255^{1152}}
\sum_{k\ge1296}[X^k]e_0^TQ(X)^{16}\mathbf1_2,
\qquad e_0=(1,0)^T.
\]

Here $N_1$ counts source sets of size 18 with matching maximum exactly one.
It is charged once. The fixed-row-set weights already contain all choices
of eight rows; there is no additional $\binom{16}{8}^{18}$ factor.

Root ran the one frozen gate in 1.223 seconds. It returned

| Quantity | FLOAT log2 value |
|---|---:|
| Positive complete-word tail | $-2379.9793367247244$ |
| Source/mark/projective/RS prefix | $4662.4561381901$ |
| Selected-class upper bound | $2282.476801465376$ |

The preceding coupled-state gate returned $2906.6054663517416$.
The new gate remains positive and does not close this class. It is not a
counterexample to the construction or to the exact two-hit return lemma.
It also makes no claim about other arities, source strata, or the whole code.

The frozen moment laws expose the next diagnostic. A two-region cross
excursion at equal occupancy $m$ costs, per region, approximately
$-126.6323814-3.4960240m$ bits. The constant live-to-live cap is
$-137.901186$. Thus sustained live transitions are larger for $m\ge4$ in
this comparator. A saved-polynomial path diagnostic should distinguish
that constant-envelope loss from a genuine limitation of the regional cap.
No further parameter fit or bank rebuild is justified by this observation alone.

Root's saved-polynomial diagnostic identified a more specific loss. Its
dominant word has nine degree-144 columns and seven empty columns. A
degree-144 transition from zero to live uses eight empty regions, then one
activation and seven live regions. The following live-to-zero column uses
the reverse pattern: seven live regions, one return, then eight empty regions.
In the quiet weight-eight branch, the output mask equals the fixed eight-row
set. Consequently these alternating early/late masks cannot arise from one
fixed row-set vector across columns. The fixed-row-set Holder comparison is
valid, but its columnwise factors permit this incompatible dominant pattern.
This diagnoses slack in the comparator, not an actual bad codeword.

## Whole-word group Minkowski at a fixed marker

The following comparison keeps each row set fixed through the complete
word. It moves a norm across groups only after retaining all columns.
Fix one scalar Chernoff marker $x\ge1$. This section makes no
coefficientwise claim about fractional roots of polynomials.

Let $\mathcal A=\{J\subseteq[16]:|J|=8\}$. For each $J$, let $N_J(T)$
be the exact single-column outer comparison count for output mask
$T\subseteq[16]$. The count omits the whole-word RS denominator. For a
vector of nonnegative regional atoms $t=(t_1,\ldots,t_{16})$, define

\[
K_J(t,x)=\sum_{T\subseteq[16]}N_J(T)
  \left(\prod_{i\in T}t_i\right)x^{\kappa(|T|)}.
\]

For the quiet-only comparison, restrict the sum to quiet masks. For the
complete comparison, use $\kappa(w)=\min(w,8)$ as above.

The matrix-valued moment laws induce a positive finite measure $\Lambda$
on all 256 chronological state transitions and their atom choices. The
initial state is zero; either final state is summed. Its mass need not be
one. It is common to every group and every row set. Let
$\mathcal L(h)=\int h\,d\Lambda$, and write $t_c(\omega)$ for the 16
atoms in outer column $c$. Define

\[
f_J(\omega)=\prod_{c=1}^{16}K_J(t_c(\omega),x).
\]

Expanding the row-set sum for the 18 labelled groups gives the common
positive numerator

\[
\mathcal L\left[\left(\sum_{J\in\mathcal A}f_J\right)^{18}\right].
\]

The $L^{18}(\Lambda)$ triangle inequality implies

\[
\mathcal L\left[\left(\sum_Jf_J\right)^{18}\right]
\le
\left(\sum_J\mathcal L[f_J^{18}]^{1/18}\right)^{18}.
\]

Each term on the right keeps one fixed $J$ through all 16 columns.
Define $Q_J(x)_{ij}$ by summing the positive atom/path weights for one
16-region column from state $i$ to state $j$, multiplied by $K_J(t,x)^{18}$.
Then chronological factorization gives

\[
\mathcal L[f_J^{18}]=e_0^TQ_J(x)^{16}\mathbf1_2.
\]

No state matrices commute in this argument. The factorization retains
their original order, and the same 16-region law repeats each column.
The selected-class first-moment comparison is therefore bounded by

\[
\frac{N_1z^{-734003}x^{-1296}}{255\,255^{1152}}
\left(\sum_{J\in\mathcal A}
  (e_0^TQ_J(x)^{16}\mathbf1_2)^{1/18}\right)^{18}.
\]

The sum over $J$ already charges all row-set choices. There is no additional
$\binom{16}{8}^{18}$ factor. The one marker $x$ is common to all terms and
is not selected separately inside the sum. The expression is a proposed
next interface; no numerical evaluation is reported here.

The common measure requirement also fixes the treatment of empty regions.
Outside a homogeneous row set, use the moment-law matrix $S(0)$, including
its relaxed off-diagonal entries. Replacing it by the smaller exact empty
matrix only in those homogeneous terms is not justified by Minkowski.
For example, two groups with row sets $\{1\}$ and $\{2\}$ can jointly
cover a two-region activation/return path. Either homogeneous assignment
leaves one required region empty and would vanish under exact empty cross
entries. A valid refinement must preserve that mixed-group coverage by a
separate positive argument; it cannot simply decrease the right-hand norms.

## Positive isolation of the nonquiet remainder

For arity eight, quiet masks have weight at most eight and every nonquiet
mask has weight at least nine. Thus a nonquiet branch can be identified
from its output mask without knowing $J$.

A safe split may use coherent row sets for the all-quiet contribution and
the previous Holder comparison for the remainder. In each group-column
polynomial, separate weights $w\le8$ from weights $w\ge9$. Attach a Boolean
flag to the latter. Propagate the flag by logical OR through the 18-group
power, ordered column matrices, and 16-column word power. Retain the same
degree tail. Flag one is a positive upper bound for the entire contribution
containing at least one nonquiet group-column.

Add that bound to the coherent all-quiet bound. Do not subtract separately
bounded totals. Earlier full/quiet FLOAT equality also does not establish
that the nonquiet remainder is below the certificate target; an explicit
positive flagged calculation is required.

That positive calculation is now complete. The saved-only OR-flag run took
4.084 seconds and gave nonquiet log2 upper bound `+2226.121421174038`.
Recombining the two flags reproduced `+2282.476801465376`, the complete
Holder comparator. Thus reusing that remainder cannot close the selected
class, even if the all-quiet contribution is improved. No provider or
parameter fit was used in this control.

## Full fixed-row-set comparison by three positive branches

A second comparison retains the same fixed row set for nonquiet columns.
Fix $A_0\in\mathcal A$ and write $O=[16]\setminus A_0$. For inside weight
$b$, let $Z_1(b)$ and $Z_2(b)$ count nonzero assignments satisfying one or
two independent Vandermonde checks, respectively. The established exact
column kernel has nonnegative branch counts

\[
q_b=Z_2(b),\qquad h_b=Z_1(b)-Z_2(b),\qquad
f_b=255^b-q_b-8h_b.
\]

The outside supports have sizes zero, seven, and eight, so these branches
are disjoint. Let $e_b(t_{A_0})$ denote the elementary symmetric polynomial
on the eight inside atoms. The complete kernel with the safe surrogate
degree has the exact decomposition

\[
\begin{aligned}
K_q&=\sum_{b=0}^8q_b e_b(t_{A_0})x^b,\\
K_f&=x^8\left(\sum_{b=0}^8 f_b e_b(t_{A_0})\right)\prod_{i\in O}t_i,\\
K_h&=x^8\left(\sum_{b=0}^8 h_b e_b(t_{A_0})\right)
        \sum_{j\in O}\prod_{i\in O\setminus\{j\}}t_i,\\
K_{A_0}&=K_q+K_f+K_h.
\end{aligned}
\]

For positive $\lambda_q+\lambda_f+\lambda_h=1$, convexity gives the
pointwise inequality

\[
K_{A_0}^{18}\le
\lambda_q^{-17}K_q^{18}+
\lambda_f^{-17}K_f^{18}+
\lambda_h^{-17}K_h^{18}.
\]

The declared trial is $(\lambda_q,\lambda_f,\lambda_h)=(99/100,1/200,1/200)$.
It is common to all atoms and row sets. The hole branch also satisfies

\[
\left(\sum_{j\in O}\prod_{i\in O\setminus\{j\}}t_i\right)^{18}
\le8^{17}\sum_{j\in O}\prod_{i\in O\setminus\{j\}}t_i^{18}.
\]

Apply these inequalities inside the same one-column positive measure.
Let $Q_{q,A_0}$ integrate $K_q^{18}$. Let $Q_{f,A_0}$ integrate
$(\sum_b f_b e_b)^{18}\prod_{i\in O}t_i^{18}$, without the scalar factor
$x^{144}$. Define $Q_{h,A_0,j}$ similarly using $h_b$ and omitting outside
position $j$. Then the complete column matrix is bounded entrywise by

\[
\widetilde Q_{A_0}=
\lambda_q^{-17}Q_{q,A_0}+
x^{144}\lambda_f^{-17}Q_{f,A_0}+
x^{144}\lambda_h^{-17}8^{17}\sum_{j\in O}Q_{h,A_0,j}.
\]

The sum over eight omitted positions is already present; do not add a
further factor of eight. Let $\mathsf C$ be the common positive $2\times2$
matrix-valued measure on local atoms that generates the regional envelope.
At full outside positions, the atom integral is
the common moment matrix $E_{18}=\int t^{18}\,d\mathsf C(t)$. At an omitted
position it is $E_0=\int1\,d\mathsf C(t)$. Neither may be replaced by a
smaller raw regional matrix after applying the common-measure argument.
The eight marked inside positions retain their ordered path calculation.

For a gap of $n$ outside positions, the exactly-one-hole contribution is

\[
H_n=\sum_{k=0}^{n-1}E_{18}^{k}E_0E_{18}^{n-1-k}.
\]

An exact degree-one flag across all outside gaps sums one omitted position
per column. It must not select one hole independently in every gap.

Finally, sum the three positive matrix branches **before** taking
$\widetilde Q_{A_0}^{16}$. Apply the preceding whole-word group Minkowski
bound with these matrices. This keeps one $A_0$ across columns, permits
their outer branches to vary, and charges no additional row-set count.
The construction yields a full-mask comparison interface. No numerical
outcome for this repair is asserted here.

## Reproducibility and scope

The twelve new compression/adapter tests pass. They include exact
noncommuting coefficient embeddings, conditional-flat mixture closure,
the two-packet identity, fixed-word first-return checks, and fail-closed
scope/provenance tests. The independently reviewed packet suite has seven
additional exact/light tests. All new evidence is proof-side FLOAT or exact
toy algebra; there was no native edit or benchmark.

| File | SHA256 |
|---|---|
| `hd_post_refresh_packet.py` | `b71bc026e32fb673478958df2fc0ec243c704024f796988d773d30c4fb607fbd` |
| `hd_post_refresh_regions.py` | `5942384325e38e77cf50db1a02b4b8ca30e50002256d9ceac25359158b2615ee` |
| `hd_post_refresh_compress.py` | `411271e8a507021ae81e030d31a1f20ea611da3c96d976543b8c389e1fb20fcc` |
| `hd_post_refresh_moment_gate.py` | `1aa9f0bf621ccae4dce7f7522b9e2db2afaa4d2edf09a14206443303e4ed5b70` |
| `hd_post_refresh_regions_float.json` | `40fcaa450f5d38447148d40f275b122c74c29d7abd341ab8c0e7899bd2223671` |
| `hd_post_refresh_moment_gate_float.json` | `81f8d1457d1900dacbc311ef1c55e12f446beaa8a9408e20c57f0173b2fa8b04` |
| `hd_post_refresh_diagnostic_float.json` | `22d24e9eda3e7ddc2911d1f84317429e6fb314f88536f58ae41d271b254dd21b` |
| `hd_post_refresh_nonquiet_gate_float.json` | `621e09c1faa97a98521b6173330d4492820d1b13ebede52a71a6fab4347b8d4b` |

The gate accepts the raw artifact and its SHA256 explicitly. It also checks
the reviewed source pins and complete histogram domain. The next bounded
step is the full coherent-row-set comparison using the three positive
branches. A separate two-packet cap can retain the independent innovation
cases if that control identifies a remaining local state loss.
