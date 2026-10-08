# Sparse multirow rank cover

The positive zero-set expansion improves the sparse existence bound without
counting all field-valued messages. Exact regional averaging then preserves
packet occupancy. These changes give useful bounds for several small active-row
classes, but they do not complete the rate-11/16, 10% distance proof.

This note concerns the GRS[32,22] outer over GF256 and the Cauchy inner from
[the diagonal-refresh proof](DIAGONAL_REFRESH_PROOF.md). It does not change the
construction, routing stage, or previously proved controls. Every probability
below refers to independently sampled ideal setup.

## Setup and the remaining class

There are \(L\) outer rows and \(N=32L\) output symbols. Each outer-coordinate
region independently permutes its rows; these permutations form the single
regional routing stage. The packet width is sixteen. The four-symbol state
starts at zero, and the final state is discarded.

Let \(C\) be the packet injection matrix and \(A\) its unscaled Cauchy emission
matrix. The packet equations are

\[
 X=D_{\mathrm{in}}\ell,\qquad
 R=D_{\mathrm{s}}S,\qquad
 Y=X+\eta AR,\qquad S'=R+CX,
 \quad\eta=2.
\]

Every diagonal entry is independently sampled from \(\mathbb F_{256}^{*}\).
The fixed maps vary through their prescribed four phases. In particular,
\(CA=0\), every four injection columns are independent, and every
\(I+\eta C_JA_J\) is invertible. The emission matrix is superregular; the
injection matrix is not.

Fix an exact routed input support \(I\). Write \(W=|I|\), and let \(a\) be
the number of nonempty outer rows before routing. Each nonempty row has at
least eleven support coordinates. Its shortened MDS dimension is its support
size minus ten. Consequently, the full shortened source dimension is

\[
 D=W-10a,\qquad 11a\le W\le32a.
\]

The existing directed class certificate covers \(W\ge\lceil .11N\rceil\).
Separate exact source counts cover one and two active rows. The unresolved
class has at least three active rows and \(0<W<.11N\).

## Positive zero-set expansion inside the rank cover

Let \(\Omega\) be the common span of all homogeneous equations already imposed
on the shortened source. A branch is alive only if its required nonzero input
and state forms remain realizable. Fix \(1\le\tau\le255\) and \(0<y\le1\).
Use the potential \(\tau^{-\dim\Omega}\).

For an actual packet output \(Y\), the positive identity

\[
 y^{\operatorname{wt}(Y)}
 =\sum_{Z\subseteq\{j:Y_j=0\}}
   y^{16-|Z|}(1-y)^{|Z|}
\tag{1}
\]

avoids enumerating exact output supports while dropping their nonzero tests.
The selected set \(Z\) imposes only its zero equations.

For every setup containing a bad nonzero word, fix one such word. Follow its
actual zero/nonzero state path and sum all selected-zero subsets from (1).
Every resulting branch has \(\dim\Omega\le D-1\). If the bad-output threshold
is \(e=\lceil N/10\rceil-1\), its total potential is therefore at least
\(y^e\tau^{-(D-1)}\). Markov's inequality gives

\[
 \Pr[\text{a bad word on }I]
 \le \tau^{D-1}y^{-e}
      \mathbb E[\text{total alive selected-zero potential}].
\tag{2}
\]

This is an existence bound over equation spaces, not a first moment over
field-valued words. Different branches may overlap; their contributions are
nonnegative. The multiplier is \(\tau^{D-1}\), not the number of messages.

### Packet transfer

Write \(w\) for the packet input weight. A selected zero set contains \(i\)
active-input positions and \(z\) empty-input positions. Let
\(\Gamma(z,\tau)\) be the arbitrary-state-image diagonal rank bound proved in
the diagonal-refresh note. On an incoming nonzero-state branch, the marginal
continuation bound is

\[
 M_w=\sum_{i=0}^{w}\sum_{z=0}^{\min(3,16-w)}
 \binom wi\binom{16-w}z
 y^{16-i-z}(1-y)^{i+z}\Gamma(z,\tau)\tau^{-i}.
\tag{3}
\]

The state diagonal supplies the root factor first. The \(i\) input-only zero
equations then consume distinct fresh input diagonals. A return to zero state
uses up to four remaining input pivots. Thus its initial bound is the same sum
with \(i<w\) and an extra factor \(\tau^{-\min(w-i,4)}\). Denote it by \(R_w\).
The omitted case \(i=w\) is impossible by invertibility of the residual state
map. In particular, \(R_0=0\).

For one active input, a stronger joint exposure gives

\[
 R_1\le y^{16}\tau^{-4}.
\tag{4}
\]

A return forces every incoming state coordinate to be nonzero because the
injection column has full support. Reveal the input diagonal first, then impose
four alive pivots using the fresh state diagonals. The returned output has all
sixteen coordinates nonzero by the exact cubic pairing. Equation (4) replaces
the earlier whole-packet return bound; it is not an additional independent
factor.

The resulting two-state matrix for a fixed input support of weight \(w\) is

\[
 T_w=\begin{pmatrix}
  \mathbf1_{w=0}+\mathbf1_{w\ge5}y^w\tau^{-4}
    &\mathbf1_{w>0}y^w\\
  R_w&M_w
 \end{pmatrix},\qquad
 T(x)=\sum_{w=0}^{16}\binom{16}w x^wT_w.
\tag{5}
\]

Both continuation entries are uniform over the common equation history. The
marginal nonzero entry may count return cases too; that relaxation increases
the bound. No exact law on formal state ranks is assumed.

## Exact averaging over active-row occupancy

Fix a set of \(a\) outer rows. If a routed support has \(h_j\) coordinates in
region \(j\), the probability that all its coordinates come from these rows is

\[
 \prod_{j=1}^{32}\frac{\binom a{h_j}}{\binom L{h_j}}.
\tag{6}
\]

The regional permutations are independent. Formula (6) preserves sampling
without replacement. It does not sample freshness flags independently.

Put \(G=L/16\), and define the matrix

\[
 U_a:=\sum_{h=0}^{a}
    \tau^h\frac{\binom ah}{\binom Lh}[x^h]T(x)^G.
\]

Dropping infeasible outer support patterns only adds positive terms. Equations
(2), (5), and (6) therefore yield the finite whole-class bound

\[
 \Pr[\text{a bad word with exactly }a\text{ active rows}]
 \le \binom La\tau^{-10a-1}y^{-e}
       e_0U_a^{32}\mathbf1.
\tag{7}
\]

This bound also relaxes the restriction \(W<.11N\). Enforcing a weight mark can
only improve it. The factor \(\tau^{-10a}\) is justified on the original exact
\(a\)-row class before infeasible supports are added; it does not assert the
same dimension for those added patterns.

## Joint return audit

The independent return implementation is
[`cauchy_joint_return.py`](theory/cauchy_joint_return.py). Its selected-zero
rank law is compatible with (1)--(7).

Before revealing any current diagonal, fix the common equation space and an
incoming state support \(H\subseteq\{1,2,3,4\}\). Impose the zero state
coordinates outside \(H\). Discard any branch on which a required input form
or state coordinate in \(H\) becomes identically zero. Omitting this initial
rank charge is conservative.

For a fixed input placement \(I\), let
\(B_I=E_I+\eta AC_I\) be the returned-output map. Given selected output zeros
\(Z\), form

\[
 K_{I,H,Z}:=\begin{pmatrix}C_{H^c,I}\\(B_I)_Z\end{pmatrix}.
\]

Choose independent pivot columns of \(K_{I,H,Z}\). Reveal all nonpivot input
diagonals first. A left inverse gives one alive equation for each pivot input
diagonal, hence a factor \(\tau^{-\operatorname{rank}K_{I,H,Z}}\). Finally,
the \(|H|\) fresh state diagonals impose the return equations and contribute
\(\tau^{-|H|}\). These are disjoint sets of random variables.

The implementation retains exactly the feasible pairs \((H,Z)\), sums their
positive marks, and **sums** over \(H\). A maximum over \(H\) would be invalid
for an existential witness that can choose its state support. Projective
arrangement counts are used only to determine feasibility, not as message
multiplicities. Averaging input placements is valid when the surrounding
coefficient transfer literally sums those placements uniformly.

The alternate coarse return law also passes the audit. After imposing a
return, the input-form image has dimension at most
\(d=w-\operatorname{rank}C_{H^c,I}\). An exact output-zero pattern is determined
by a canonical basis for the span of its zero forms. For basis rank \(k<d\),
there are at most \(\binom{16}k\) choices, each costing \(\tau^{-k}\). The
minimum feasible output weight supplies the remaining \(y\)-factor. This is
a replacement exact-support cover with uniform continuation, not pointwise
domination of every selected-zero branch.

## Executed screens and remaining loss

At requested size \(2^{20}\), padding gives \(L=47680\), \(N=1525760\), and
\(e=152575\). The authoritative regional screen is
[`cauchy_regional_rank_screen.py`](theory/cauchy_regional_rank_screen.py).
The following are floating evaluations of upper bounds, not directed
certificates. Negative margins mean that the tested upper bound exceeds one;
they are not lower bounds on construction failure.

| Active rows | Initial zero-set return law | Joint two/three-input return law |
|---:|---:|---:|
| 2 | 70.42 bits | not recorded here |
| 3 | not recorded here | 102.21 bits |
| 4 | 36.66 bits | 71.98 bits |
| 8 | -34.82 bits | -19.29 bits |
| 16 | -163.22 bits | -143.75 bits |
| 32 | not recorded here | -363.71 bits |
| 64 | -842.71 bits | -752.85 bits |

The initial entries for two, four, sixteen, and sixty-four rows were
independently replayed. The joint-return entries were executed by the parent
agent after the mathematical audit above.

Changing packet width or outer length does not presently remove the gap.
The parent agent's width-eight screen, with only the older one-input joint
return law, gave 43.70, -28.65, and -154.20 bits at four, eight, and sixteen
rows. GRS[64,44] gave 80.56 bits at eight rows but -67.02 at sixteen.
GRS[128,88] gave 64.38 bits at sixteen rows but -275.46 at thirty-two.
These exploratory comparisons do not select a new construction.

## Next finite refinement: actual state-support masks

The two-state return law forgets how a sparse incoming state was created.
Tracking its actual coordinate support can retain this cost. Empty-input
packets preserve the support exactly under diagonal refresh.

Fix an incoming support \(H\), an intended next support \(K\), an input
placement \(I\), and selected output zeros \(Z\). The current independent
diagonals multiply required nonzero forms. Their numerical joint coefficient
matrix is

\[
 \mathcal B_{H,K,I,Z}=
 \begin{pmatrix}
  \eta A_{Z,H}&E_{Z,I}\\
  E_{K^c,H}&C_{K^c,I}
 \end{pmatrix}.
\tag{8}
\]

Independent pivot columns and the alive diagonal lemma give
\(\tau^{-\operatorname{rank}\mathcal B_{H,K,I,Z}}\), uniformly over all
dependencies among the underlying source forms. Feasibility requires that its
kernel does not force any current diagonal variable to zero. It must also
leave each intended next-state coordinate in \(K\) nonzero. Since GF256 is
larger than the number of required nonzero forms, individual nonvanishing
tests suffice for simultaneous feasibility by a union bound over proper
hyperplanes.

For a return, (8) reduces to the audited exponent
\(|H|+\operatorname{rank}K_{I,H,Z}\). For no input, only \(K=H\) is possible.
For a one-input birth, \(K\) is the full four-coordinate support. This gives
a finite, directly testable mask transfer. The bounded reference
[`mask_state_rank.py`](theory/mask_state_rank.py) implements weights zero,
one, and two. Weight three and above explicitly raise an unsupported-case exception.
No evaluated whole-class improvement from this refinement is claimed here.

The interface `local_mask_transfer(weight,y,tau,phase)` returns a positive
sixteen-by-sixteen matrix. Rational marks give exact rational entries. Each
matrix averages all \(\binom{16}{w}\) input placements; a surrounding support
generating function must still multiply it by \(\binom{16}{w}\). The phase is
explicit, and there is no assumed distribution on incoming masks.

The companion `input_profiles(weight,phase)` exposes complete integer
\((r,z,c)\) lists for every matrix entry. Each term evaluates to
\(c\,y^{16-z}(1-y)^z\tau^{-r}/\binom{16}{w}\). This convention includes
zero-state births: their inactive output coordinates contribute the positive
binomial identity, which collapses to \(y^w\). The fixed-placement interface
`placement_profiles(phase,incoming,inputs)` omits placement normalization.

For an empty packet, the only nonzero entry in incoming row \(H\ne\varnothing\)
is

\[
 T_{0}[H,H]=\sum_{z=0}^{|H|-1}\binom{16}z
       y^{16-z}(1-y)^z\tau^{-z}.
\]

For one or two inputs, the reference enumerates selected-zero sets with fewer than
\(|H|\) empty-input roots. It computes the rank and feasibility of (8) by
batched Gaussian elimination over the exact GF256 multiplication table.
This uses no floating-point linear algebra. For one input, zero-state births
have full next-state support. A one-input return is possible only from full
incoming support, and its entry is exactly \(y^{16}\tau^{-4}\).

For two inputs, the enumeration includes cancellation outside the incoming
support. No one-input support restriction is reused. A zero-state packet
cannot return because its two injection columns are independent. Its
nonzero next masks are checked independently by enumerating all 255
projective ratios of its two nonzero inputs. Two-input returns are checked
against the independent projective-arrangement reference after summing the
incoming masks.

Simultaneously translating every packet label and residue node by an element
of the embedded GF4 preserves their pairwise differences. In the chosen
ordering, translation permutes input indices and state coordinates by XOR
with the same two-bit index. The implementation transports both incoming and
next masks under this permutation. It reduces the 120 input pairs to 36
orbit representatives, retaining every distinct translate exactly once.
Pairs with a nontrivial stabilizer are not multiplied by four. The phase
remains fixed throughout this reduction; phases are not averaged.

A conservative fallback must not divide a nonzero-state upper bound among
unknown masks. Broadcasting it to every mask is valid but duplicates mass.
An additional pooled nonzero state can instead forget the mask and retain the
original two-state bound until the next return to zero. Neither fallback is
implemented in this bounded checkpoint.

## Replay and scope

The new local reference is
[`sparse_multirow_zero_set.py`](theory/sparse_multirow_zero_set.py). Its three
tests pass: the positive identity, exact rational support summation, and
twenty-four exhaustive GF16 common-quotient cases. Those cases vary the mark
and rank tilt and expose every current state/input diagonal.

The nine mask tests also pass. They compare batched ranks and feasibility
against independent scalar elimination, check exact empty-packet masks and
one-input births/returns, and verify two-input births and returns in all four
phases. They compare the orbit reconstruction with direct translated
placements and preserve the four original one-input profile hashes. They
include ninety-six exhaustive GF4 comparisons of the joint diagonal rank
charge with arbitrarily dependent source forms. The high-weight rejection is
tested explicitly.

```sh
python -B research/workstreams/char2_field_spin/theory/sparse_multirow_zero_set_tests.py
python -B research/workstreams/char2_field_spin/theory/test_cauchy_joint_return.py
python -B research/workstreams/char2_field_spin/theory/mask_state_rank_tests.py
```

No native implementation, benchmark, remote run, or whole-code certificate
was produced. The next useful gate adds a justified pooled fallback for
weights at least three, then composes the tested mask tables with exact regional occupancy.
A directed certificate is appropriate only if the full sparse-class bound
closes.

## Independent audit: dependent source forms and the pooled composition

The numerical matrix rank in (8) does not assume that the underlying source
forms are independent. The alive predicate is essential to this statement.
The following argument establishes the required conditional bound directly.

Fix a common equation space `Omega`, and work in its quotient. Let
`f_1,...,f_m` be arbitrary source forms, each required nonzero on a witnessing
source. They may coincide or satisfy other linear relations. Let `B` be a
fixed numerical coefficient matrix of rank `r`, and let the entries of `D`
be fresh independent nonzero GF256 scalars. Append the equations `BDf=0`.
Write `Delta` for their additional rank. On a branch where a required form
vanishes identically on the resulting kernel, define the alive potential to
be zero. Then, for `1 <= tau <= 255`,

\[
 \mathbb E_D[\tau^{-\Delta}\mathbf1_{\rm alive}]
 \le\tau^{-r}.
\tag{9}
\]

Choose `r` independent columns of `B`, and apply row operations so these
columns form the identity. Expose every nonpivot diagonal first. Each of the
`r` remaining equations now has the form `d f + g`, using its own fresh
pivot scalar. All earlier equations are independent of that scalar. In the
quotient by those earlier equations, discard the branch if `f` is zero.
Otherwise there are two cases:

- If `g` is independent of `f`, every scalar choice adds one rank and keeps
  `f` nonzero on the new kernel. The potential contracts by `tau^-1`.
- If `g` lies in the span of `f`, at most one nonzero scalar makes the new
  equation vanish. Every other scalar forces `f=0` and kills the branch.
  The surviving zero-rank mass is at most `1/255 <= tau^-1`.

Discarding branches that kill other required forms only decreases the mass.
Iterating the conditional bound proves (9). This proof is uniform over the
common history and every dependency among the source forms.

Two examples locate the precise assumption. If `B=I_2` and `f_1=f_2=f`, the
new formal rank is only one. Omitting the alive predicate would falsely
compare `tau^-1` with `tau^-2`. The correct alive mass is zero, because the
equations force the required form `f` to vanish. If instead `B=(1,1)` and
`f_1=f_2=f`, the branch survives exactly when the two diagonals agree. Its
zero-rank mass is `1/255`, which satisfies the one-pivot bound.

For each intended successor mask `K`, apply (9) to its own joint matrix (8).
Impose zero successor coordinates outside `K`, and retain the nonzero tests
inside `K`. Summing the resulting upper bounds over `K` is a positive
existence cover. It need not be a stochastic transition matrix: one source
kernel may contain witnesses with several different successor masks. The
sum may overcount such a kernel, but it cannot omit its witness.

The table's numerical-kernel tests provide a necessary feasibility filter.
They do not establish that every common-history source kernel realizes that
coefficient pattern. Additional dependencies can kill the actual branch;
the uniform alive bound already permits that loss. No distribution on
formal state ranks or on unknown coordinate masks is used.

An independent first-principles check used three deliberately dependent
source histories over GF4. For each history it enumerated all 81 current
nonzero four-coordinate diagonals and all 64 source vectors. It tested
actual nonzero input/state values, selected output zeros, and the actual
successor support. Across every two-coordinate zero set, successor mask,
and tilt `tau` in `{1,2,3}`, all 144 comparisons satisfied the summed rank
upper action. This check tests witness support, not merely agreement
between two matrix-rank implementations.

### Seventeen-state regional implementation

The read-only audit of
[`mask_regional_screen.py`](theory/mask_regional_screen.py) found no
mathematical mismatch in its composition. States zero through fifteen are
known actual masks. State sixteen is pooled nonzero state. Its only
destinations are known zero and pooled nonzero; the implementation never
redistributes pooled mass into specified masks. Unsupported nonzero masks
use the same uniform-history two-state bound and forget their mask.

The fallback calls `rank_return_moment`, not the fixed-word return API.
Taking the minimum of two whole-return rank bounds is valid. The marginal
continuation may overcount return cases, but remains an upper bound on the
actual nonreturn branch. Its use does not assert that the added mass is an
actual state distribution.

The implementation takes an entrywise maximum over the four fixed phases,
not an average. Its base packet matrices already average all input
placements of each weight. Normalized block composition therefore uses

\[
 \frac{\binom{n_1}i\binom{n_2}{h-i}}{\binom{n_1+n_2}h}
 F_i^{(1)}F_{h-i}^{(2)}.
\]

No additional `binom(16,w)` factor belongs on these normalized base
matrices. Matrix order remains chronological. The regional sum then uses
`binom(a,h) tau^h`, and the final source factor is
`binom(L,a) tau^(-10a-1) y^-e`, as in (7).

The dimension `W-10a` is applied to the original feasible exact-`a` class.
Adding positive terms for patterns with fewer active rows or invalid row
weights is conservative; it does not assign that dimension to the added
patterns. Floating evaluation and optimization remain screens. This audit
establishes the conditional upper-bound interface, not the required
whole-class 40-bit margin.

## Independent audit of the realized unequal-window obstructions

Two saved K20 routes now have deterministic low-weight subspace witnesses.
This changes the immediate task: improving the seventeen-state upper bound
cannot repair these particular routes. The witness works for every choice of
the four-state linear packet coefficients, including fresh full matrices.
It is not, by itself, a lower bound on failure probability for a uniformly
random route.

For each of the 32 coordinate regions, fix one packet-aligned window `W_j`.
Let `c_i` count the coordinates of outer row `i` whose routed positions lie
inside their respective windows. Shortening GRS[32,22] row `i` to those
coordinates leaves dimension `max(c_i-10,0)`. Different rows have independent
source coordinates, so the shortened source dimension is

\[
 D_{\rm win}=\sum_i\max(c_i-10,0).
\]

Fix the route and every packet coefficient. The state at the end of each
window is a linear function of the shortened source. Requiring all 32 exit
states to vanish imposes at most `32*4=128` equations. If `D_win>128`, their
common kernel contains a nonzero source.

For such a source, the initial state is zero, and every input outside the
windows is zero. Induct over the windows: after each constrained zero exit,
the state and outputs remain zero through the following gap. Packet
alignment prevents a window boundary from cutting through an emission
block. The entire output is therefore supported on the chosen windows.
Identity current-input feedthrough makes the causal inner injective, so a
nonzero shortened source has nonzero output. This argument needs linearity,
zero preservation, finite state dimension, and causal injectivity; it does
not need invertible state refresh or a particular distribution of the
packet coefficients.

The literal fixtures are in
[`window_dimension_witnesses.py`](theory/window_dimension_witnesses.py).
The audited geometry is `L=47680`, `N=1525760`, and bad-output limit
`e=152575`. Both window families have total length 152064.

| PCG64 route seed | Full padded shortened dimension | Dimension after excluding final 18 rows | Surviving dimension after exit constraints |
|---:|---:|---:|---:|
| 1 | 136 | 136 | at least 8 |
| 2 | 133 | 132 | at least 4 |

The fixed-fixture verifier was replayed. An independent reconstruction used
inverse permutations and direct window-position comparisons, rather than
the verifier's concatenation and histogram calculation. Both methods give
the displayed dimensions, and check the aligned window endpoints and total
length. No optimization was rerun during this audit.

### Padding scope

The padded source has `22L=1048960` coordinates, whereas the requested
source has `2^20=1048576`. Under the natural row-major padding convention,
zeroing the last eighteen complete rows leaves

\[
 22(47680-18)=1048564\le2^{20}
\]

source coordinates. This also zeros twelve genuine input coordinates,
which is allowed when proving a bad subcode exists. The table recomputes
the shortened dimension on these retained rows. Both witnesses therefore
apply to the requested source size under that padding convention, not only
to the larger padded source. This does not assert invariance under an
arbitrary different linear embedding of the unpadded source.

### Bounded repair comparison

One extra redundancy coordinate and one extra state coordinate both defeat
these specific dimension inequalities. With outer redundancy eleven, the
same retained hit patterns have shortened dimensions 74 and 79. With five
state coordinates, the 32-exit budget is 160. Neither observation proves
that a changed construction has no other bad window family.

The cheapest direct candidate to screen is GRS[33,22], of rate exactly 2/3,
with the current sixteen-wide, four-state Cauchy inner. It adds one parity
evaluation and preserves the verified fast inner maps. A single new parity
transpose input needs at most seven nontrivial GF256 basis products per
four source rows. Thus an unoptimized outer upper count is `41+7=48`, and
the complete leading product estimate is

\[
 \frac{48}{88}+\frac{33}{22}
 =\frac{45}{22}\approx2.04545,
\]

about 6.51% above `169/88`. This is an arithmetic upper estimate, not a
native result. A 33-region construction changes the window budget and
route geometry, so it requires a new obstruction screen before a positive
proof or implementation effort.

State five retains the outer rate but changes the packet algebra. It does
not directly inherit the four-residue five-product circuits or all the
four-state local lemmas, and its MDS emission has smaller minimum weight.
It is therefore a less direct first repair despite raising the exit budget.
Optimizing the existing sixteen-wide kernel, or replacing its state diagonal
by a fresh four-by-four matrix, does not change the dimension witness.

Changing how regional permutations are grouped requires a new spread
condition. Sharing a permutation across all regions is particularly bad:
aligned windows can then select the same rows in every region. A grouping
change should be justified against the capacity functional above, not
assumed to eliminate it.

The current conclusion is limited but material: these two realized routes
are invalid for the four-state candidate, even at the requested unpadded
size. Two successful searches among a few pseudorandom seeds are not a
certified setup-failure rate. The next bounded task is the parity-eleven
capacity screen; no changed construction is promoted here.
