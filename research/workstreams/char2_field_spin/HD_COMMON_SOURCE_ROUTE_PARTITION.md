# A common-source route partition for the column-first bound

This note derives an exact route-event probability and a conditional mean
operator interface. It does not give a distance bound. Six tiny exact tests
pass; no production probability grid, regional bank, or native run was made.
The late-alpha application is the column transfer in
`HD_TWO_SUM_SYNDROME_RANK.md`.

## A source pair spans two emission blocks

Fix one region's source group set `G`. Its prescribed even or odd cyclic
matching contains `s` pairs with one endpoint in `G` and `f` pairs with both.
Write `N=4096` for the number of source pairs. The route permutes these pairs
uniformly and swaps the two endpoints independently and fairly.

In `HdGroupedInner.h`, `innerInput4<J>` loads source pairs `route[2*J]` and
`route[2*J+1]`. Their endpoints fill record lanes `(0,1)` and `(2,3)`.
Each lane belongs to a different RS[16,4] emission/feedback block; `J` is the
coordinate within that block. Thus a pair does **not** occupy two coordinates
of one emission block. A packet has 16 pair slots coupling blocks 0 and 1,
and 16 further slots coupling blocks 2 and 3. The existing `sparse_geometry`
formula independently records the same occupancies `(a+f,b+f,u+g,v+g)`.

An odd cyclic source pair, including the wrapping pair, has this same
destination geometry. The odd matching changes `(s,f)`, not the slot layout.

Call a regional route **good for G** when every one of its 512 emission
blocks contains at most four positions from `G`. This event refers to all
potential source positions. It is the same event for every actual outer
mask `U` contained in `G`.

## Exact positive coupon polynomial

In one coupled half, let `a,b,d` count singleton pairs oriented toward its
left block, singleton pairs oriented toward its right block, and full pairs.
Its block occupancies are `a+d,b+d`. Define the integer polynomial

\[
H_4(X,Y)=
\sum_{\substack{a,b,d\ge0;\ a+b+d\le16\\a+d\le4,\ b+d\le4}}
\binom{16}{a,b,d,16-a-b-d}X^{a+b}Y^d.
\]

There are 256 coupled halves per region. The exact probability is

\[
p_G=\Pr[E_G]
=\frac{[X^sY^f]H_4(X,Y)^{256}}
 {2^s\binom{4096}{s}\binom{4096-s}{f}}.
\]

The denominator counts placements of the fixed pair types and the fair
orientations of singleton pairs. Swaps of full or empty pairs are invisible
to this event and cancel from numerator and denominator. Without the cap,
the half polynomial is `(1+2X+Y)^16`; its regional coefficient equals that
denominator. Coefficient extraction enforces the fixed histogram, so no
independence between packet occupancies is asserted.

This also gives simple checks. If `s+f<=4`, the route is always good. One
full pair always satisfies a block cap of one: its endpoints occupy distinct
blocks. Two singleton pairs need not satisfy that cap.

For a cheap analytic upper bound on the bad probability, write
`(x)_j=x(x-1)...(x-j+1)`, with `(x)_0=1` and `(x)_j=0` when `j>x`.
Counting five-position coupons in each block gives

\[
1-p_G\le
\min\left\{1,
512\binom{16}{5}
\sum_{u=0}^{5}\binom5u
\frac{(f)_u(s)_{5-u}}{(4096)_5\,2^{5-u}}
\right\}.
\]

Any bad block contains at least one such coupon. This is a union/first-moment
bound, not an independence approximation. The helper implements the exact
positive coefficient formula; no production parameter values were evaluated.

## A conditional operator must track G and the actual mask together

Fix one actual regional mask `U` contained in `G`. Six pair types describe
their joint support, after using the independent fair swap:

| Type | Potential endpoints in G | Actual endpoints in U |
|---|---:|---:|
| 0 | 0 | 0 |
| S0 | 1 | 0 |
| S1 | 1 | 1 |
| D0 | 2 | 0 |
| D1 | 2 | 1 |
| D2 | 2 | 2 |

Let their fixed counts be `n=(n0,nS0,nS1,nD0,nD1,nD2)`. Then
`s=nS0+nS1` and `f=nD0+nD1+nD2`. Labels within each type remain fixed at
the source; uniform pair permutation makes their destination placements
exchangeable. This does not resample `G` or `U`.

For one slot, use variables `z_type` for these six counts, `pL,pR` for
potential occupancy, and `aL,aR` for actual occupancy. The positive slot
polynomial is

\[
\begin{aligned}
K={}&z_0
+\tfrac12z_{S0}(p_L+p_R)
+\tfrac12z_{S1}(p_La_L+p_Ra_R)\\
&+z_{D0}p_Lp_R
+\tfrac12z_{D1}p_Lp_R(a_L+a_R)
+z_{D2}p_Lp_Ra_La_R.
\end{aligned}
\]

Use `K^16` independently for the two coupled halves, with their own occupancy
variables and the same six type variables. For each resulting packet term,
replace its actual occupancy monomial by the existing **unpowered** packet
cap matrix `M(o0,o1,o2,o3)`. Keep its exact nonnegative coefficient. Put the
term in `Q_good` when all four potential occupancies are at most four;
otherwise put it in `Q_bad`. The actual cap can be symmetric in block
occupancies, but its matrix indices and packet order are not interchangeable.

Write `Q_all=Q_good+Q_bad`. The event-weighted regional cap means are

\[
V_{\rm good}(G,U)=
\frac{[z^n]Q_{\rm good}^{128}}{\binom{4096}{n}},\qquad
V_{\rm bad}(G,U)=
\frac{[z^n]F_{\rm bad,128}}{\binom{4096}{n}},
\]

where a positive bad-flag recurrence computes the second numerator:

\[
F_{\rm good,0}=I,\quad F_{\rm bad,0}=0,\qquad
F_{\rm good,j+1}=F_{\rm good,j}Q_{\rm good},
\]
\[
F_{\rm bad,j+1}=F_{\rm bad,j}Q_{\rm all}
                 +F_{\rm good,j}Q_{\rm bad}.
\]

All multiplications retain chronological order. This recurrence avoids
subtracting two independently rounded upper matrices. Neither the route
multinomial nor the orientation probabilities are powered.

Divide `V_good` by `p_G` and `V_bad` by `1-p_G` to obtain conditional cap
means; omit a zero-probability branch. The same `p_G` applies to every `U`.
The resulting matrices bound expectations of cap products, not probabilities
of state transitions or individual fixed refresh realizations.

## Use inside one late-alpha column

Regional routes are independent. One may condition on their vector of
good/bad flags, then compose the corresponding conditional matrices in the
original sixteen-region order. For each fixed outer comparison input,
`U_i` is its actual support in region `i`. Sum all comparison inputs belonging
to the same column coefficient **before** taking an entrywise alpha power.
Conditional Jensen then puts the common event probability outside alpha.

Alternatively, use just two column events: all sixteen regional routes are
good, or at least one is bad. A positive event-flag transfer computes the
second conditional mean without a `2^16` enumeration. For the alternating
matching, the first probability is `p_even^8*p_odd^8`. A fixed `G` can have
different even- and odd-matching histograms. Their probabilities must not
be replaced by a single matching's value.

The six-type reduction suffices for each fixed mask's conditional regional
mean. Linearity permits those means to be summed inside the complete column
coefficient. It does not make different masks independent. Directly averaging
separately powered mask means would return to the earlier, looser ordering.
Row sets and source groups remain fixed across all sixteen outer columns.

This partition cannot zero a regional 0-to-0 entry. It only ensures that a
nonempty packet starting at zero has nonzero feedback. Later packets can
return the state to zero. It also does not exclude the 0-to-1/1-to-0 boundary
paths that dominated earlier controls: low-occupancy inputs already satisfy
the good event. Therefore a small bad-event probability alone does not show
that this partition removes the current proof obstruction.

## Isolated reference and next gate

`theory/hd_shifted_common_source_partition.py` provides the exact scalar
`good_probability(s,f)` and a deliberately tiny six-type matrix reference.
The latter rejects block lengths above two. It is not a production builder.

Six tests in `hd_shifted_common_source_partition_tests.py` passed in
0.222 seconds. They exhaust tiny oriented source histograms, check the
coupled-block geometry, and compare a two-packet conditional operator with
an independent labelled-route enumeration using noncommuting rational
matrices. Identity caps recover the exact event probability. Zero-probability
branches and invalid dimensions are also checked.

Next, derive a truncated six-type production builder for one frozen class,
or a proved envelope over those types. Compare its conditional Jensen column
bound with the unconditional column mean under the same unpowered provider.
Do not start another parameter grid before measuring whether this common
event changes the troublesome state-path terms.

## A rare setup event across all 256 regions

Fix a source set `G` of size 18 with `Fmax(G)=1`, where `Fmax` is the maximum
of its even- and odd-matching full-pair counts. Write those counts as `a,b`.
Their possible pairs are `(0,1),(1,0),(1,1)`. Let `p_f` be the exact bad-region
probability for potential histogram `(18-2f,f)`, as defined above. The setup
uses independent regional pair permutations and fair swaps. Consequently,
for this fixed `G`, its total number `B_G` of bad regions has generating
function

\[
\mathbb E Y^{B_G}
=((1-p_a)+p_aY)^{128}((1-p_b)+p_bY)^{128}.
\]

This is an exact law in the ideal setup ensemble. It does not assert
independence between events for different source sets. The reproducible
native PRNG is not substituted for this ensemble in the probability claim.

There are exactly

\[
N_1=154355220469939768130077702840442839490778956275433472
\]

such source sets. With `p*=max(p0,p1)=p0`, the factorial-moment bound gives

\[
\Pr[\exists G:\ |G|=18,\ F_{\max}(G)=1,\ B_G\ge b]
\le N_1\binom{256}{b}(p_*)^b.
\]

Indeed, `1_{B_G>=b}<=C(B_G,b)`. The expectation of the latter counts
`b`-subsets of independent regional flags. A union bound then counts each
source once. No constituent-row-set or message multiplicity belongs in this
setup event.

The exact saved coupon probability is

\[
p_*=
\frac{4288483414243575547879649097273422783127}
     {69648387130532914430931262677779163653530980352}.
\]

For `b=13`, exact integer arithmetic proves

\[
N_1\binom{256}{13}(p_*)^{13}<2^{-63}.
\]

Its diagnostic FLOAT logarithm is `-63.684904805927545`. For comparison,
threshold 12 gives `-43.96208335792653`, between `-44` and `-43`. The selected
threshold 13 leaves a larger union budget and retains setups with at most
12 bad regions per source. This lemma bounds only the specified rare setup
event; it is not a distance certificate.

### Composing the remaining distance term

Do not condition on the event that *all* source sets have at most 12 bad
regions and then assume independent regional routes. Instead work source
by source. Let `S_G(Pi)` be a valid conditional first-moment upper bound on
failure for source `G`, including the union over its fixed row sets. Put
`E_G={B_G<=12}`. The class failure probability is at most

\[
\varepsilon_{13}
+\sum_G\mathbb E_\Pi[1_{E_G}\min\{1,S_G(\Pi)\}],
\qquad
\varepsilon_{13}=N_1\binom{256}{13}(p_*)^{13}.
\]

For `0<alpha<=1`, each remaining summand is at most

\[
\Pr[E_G]^{1-\alpha}
\bigl(\mathbb E[1_{E_G}S_G]\bigr)^\alpha
\le \bigl(\mathbb E[1_{E_G}S_G]\bigr)^\alpha.
\]

The unpowered expectation must therefore retain the common-source flags
through the whole word. For each actual regional support `U` contained in
`G`, use the event-weighted matrices

\[
M_g(G,U)=\mathbb E[1_{\rm good}T_{G,U}],\qquad
M_b(G,U)=\mathbb E[1_{\rm bad}T_{G,U}],
\]

where `T` is the ordered, unpowered regional cap after averaging the
independent amplitude and inner keys. A pointwise good-route cap `G_U` and
an unconditional mean cap `R_U` imply, entrywise,

\[
M_g(G,U)\le\min\{(1-p_f)G_U,R_U\},\qquad M_b(G,U)\le R_U.
\]

No subtraction of upper bounds is justified. A pointwise bad-route cap,
when available, also permits `M_b<=p_f B_U`. Preserve the exact zero-input
cross entries: an empty input cannot create or remove the zero state.

Replace each regional matrix by `M_g+Y M_b`, in its original position.
After summing the outer masks, compose all 16 column matrices and select
total `Y` degree at most 12. Resetting this degree at each column is invalid.
The matrices retain both boundary states; the initial state is zero and the
final vector is `(1,1)`.

For the homogeneous arity-eight class, fixed-row-set Holder weights may be
used pointwise at every fixed flag word. The indicator is independent of
the row sets. Thus the existing raw `W8` weights, surrogate input degree
`min(output_weight,8)`, whole-word degree floor 1296, and normalization
`255^-1152` remain applicable. The projective factor is charged once.
Only after both whole-word degree selections may one apply the displayed
scalar fractional inequality. This describes a composable interface, not
an implemented full-state distance gate.

### Marked good-packet caps and the remaining obstruction

The generic matrix in `HD_COMMON_SOURCE_ANALYTIC_GATE.md` discards all output
marks. More of the emission distance can be retained. Let `l` be a nonzero
packet's actual occupancy, with at most four inputs in each of its four
emission blocks. Then `1<=l<=16`. With `D=256^16-1` and `beta=D/255^16`, a
valid marked packet cap is

\[
A_l(z)=
\begin{pmatrix}
0&\beta z^l\\
z^{13\lceil l/4\rceil-l}/D&
\beta z^{l+13-2\min(l,4)}
\end{pmatrix}.
\]

Starting from zero emits exactly the input support. Every occupied block
has nonzero feedback, so the zero-to-zero entry is zero. Returning from a live
state requires the unique reference value `Q=-CX`, of mass `1/D`. Every
occupied block then has a nonzero state and emits at least `13-l_i`
coordinates, giving the lower-left exponent. For a live incoming state,
activating a state block changes its output lower bound from `l_i` to
`13-l_i`. Since `l_i<=4`, every additional active state block costs at least
five coordinates. Activating one block with largest `l_i` minimizes this
bound and gives the lower-right exponent. Fresh refresh contributes the
same established pointwise density factor `beta`.

Empty packets have cap `diag(1,beta*z^13)`. Ordered products of these caps
retain the weight cost of empty gaps between state birth and return.
Deleting the empty packets gives a weaker, position-free bound. Neither
bound zeros a regional return entry or forbids excursions across regional
boundaries.

The current full-state diagnostic is dominated by four excursions in a
full column, with transition counts `(8,4,4,0)` and occupancy 18 in each
active region. These transitions can occur on good routes. All 18 inputs
may lie near one end of a region while still respecting the block cap.
Therefore the rare-setup lemma and improved good `00` cap alone do not
remove that dominant term. A useful next gate must average endpoint
positions or empty gaps in the good `01/10` entries. The small marked
per-packet improvement is not evidence that such a gate will close.

### Exact implementation and reproduction

`theory/hd_shifted_common_source_bad_count.py` recomputes the exact coupon
probabilities, checks the pinned source-count partition, and verifies the
strict `2^-63` inequality by an integer comparison. Six tests pass in
0.030 seconds. They compare labelled flag enumeration with its polynomial,
check the factorial-moment identity and inequality, distinguish a global
cut from per-column cuts, and validate the production scope and margins.
No regional provider or numerical search was run.
Independent source review found no issue with the exact lemma or the three
marked good-packet entries.

From `theory`, with a fresh output path:

```powershell
C:/Python314/python.exe -m unittest hd_shifted_common_source_bad_count_tests
C:/Python314/python.exe hd_shifted_common_source_bad_count.py --output hd_shifted_common_source_bad_count_exact.json
```

| File in `theory` | SHA-256 |
|---|---|
| `hd_shifted_common_source_bad_count.py` | `7df83f4f85134c58c9f8d88d634c81f0b557530647f8952e9b72185e67db6c37` |
| `hd_shifted_common_source_bad_count_tests.py` | `80ed9e15eaf3c0dbc7db3e492891032165c9c0645e3cf9085fbfcc3fa0622a4d` |
| `hd_shifted_common_source_bad_count_exact.json` | `8c5445bf2266457799a5d56fec6cff9ff9dd03e8aec9511427bf6264cdefbbd8` |
