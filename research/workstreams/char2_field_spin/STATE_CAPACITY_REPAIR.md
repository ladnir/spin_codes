# Eight-state Cauchy repair: algebra, return rank, and cost

The four-state GRS[32,22] candidate has deterministic bad-route witnesses. Increasing the state capacity changes that obstruction without changing the outer rate or adding a routing stage. The eight-state construction below has exact scalar maps and a useful local return bound. With GRS[256,176], it has a [complete K20 certificate](FIELD256_K20_CERTIFICATE.md) and a [validated native implementation](FIELD256_NATIVE.md). The later [GRS128/87 checkpoint](FIELD128_K20_CERTIFICATE.md) has its own complete certificate and [faster K20 measurements](FIELD128_NATIVE.md).

The preferred variant uses a symmetric pair, `C=A^T`. Its return always costs eight independent diagonal rank pivots. A degree-seven asymmetric alternative has stronger one-input output weight, but weaker three-input returns. No nonzero GF256 emission scale makes every allowed-output update invertible for that asymmetric alternative.

## Fixed maps and packet interface

Use AES GF256 and its embedded subfield `H=GF16`. In AES byte notation, fix

\[
W=\{0,1,92,93,224,225,188,189\}\subset H.
\]

This is an additive three-dimensional GF2 subspace, not a GF8 subfield. Its vanishing polynomial is

\[
L(X)=X^8+X^4+X^2+X,\qquad L'(X)=1.
\]

The four packet phases remain the sixteen-point sets `H+delta`, with `delta` in `{0,2,6,4}`. Fix the eight residue nodes `V=16+W`, disjoint from all 64 packet labels. Let `Q(X)=L(X+16)`. For packet label `u` and node `v`, define

\[
A_{u,v}=\frac{Q(u)}{u+v}
       =(u+v)^7+(u+v)^3+(u+v)+1,
\qquad C=A^T.
\]

The numerator is nonzero at every packet label. Thus `A` is a row-scaled Cauchy matrix; every square minor of size at most eight is nonzero. Its image is a `[16,8,9]` MDS code. The transpose `C` is superregular as well.

The input `X` already includes the independent nonzero outer-coordinate diagonals. The packet has state `S` in GF256^8. With a fresh independent nonzero eight-coordinate diagonal `D_s` and fixed `eta=2`, apply

\[
R=D_sS,\qquad Y=X+\eta AR,\qquad S'=R+CX.
\]

All map entries are GF256 scalars acting on the record components. No general GF128 multiplication is required.

Every entry of `CA` sums a polynomial of degree at most fourteen over an affine GF16 set. These sums vanish, so `CA=0`. Consequently `CY=CX`. For fixed diagonals, a second zero-initial-state pass follows the same state trajectory and recovers the input. Identity current-input feedthrough also proves causal invertibility directly. The discarded final state does not affect either argument.

The extended adjoint is

\[
x^*=z+C^Th,\qquad s^*=D_s(h+\eta A^Tz).
\]

The diagonal acts after the backward accumulation. The forward packet is not its own extended adjoint.

## Pure-return rank bound

Fix the route, outer supports, all previous randomness, and the current common source quotient. Its current input and state coordinates are linear forms. They may be dependent. The current input and state diagonals remain fresh and independent.

Let `I` be the active input set, with `w=|I|>0`. First fix the actual nonzero support `H_0` of the incoming state, with `h=|H_0|>0`. Impose the deterministic equations `S_(H_0^c)=0`, and retain the requirement that every input form and every state form in `H_0` is nonzero on the witness.

A return `S'=0` imposes the numerical coefficient matrix

\[
 B=\begin{pmatrix}
 I_h&C_{H_0,I}\\
 0&C_{H_0^c,I}
 \end{pmatrix}
\]

on the independently scaled live forms. If `w<=8-h`, the lower block has full column rank. It forces every active input value to zero, so this branch is impossible. Otherwise the lower block has rank `8-h`, and `rank(B)=8`.

The independent-diagonal pivot lemma therefore charges `tau^-8`, for `1<=tau<=255`, uniformly over the common source quotient. Form dependencies do not invalidate this statement: at a dependent pivot, all but at most one diagonal value force a required-live form to zero. That surviving value has probability at most `1/255<=1/tau`. An independent pivot adds one equation rank.

This is an alive-rank bound, not an assumption that the state is uniformly distributed. No current diagonal is used first in an output-zero bound and then charged again here.

### Retaining the rank of the incoming zero equations

It is unnecessary to pay one unweighted branch for each possible `H_0`. Write `Z=H_0^c`. A feasible `Z` is a flat of the vector matroid defined by the state forms: any state form in their span must also vanish. A rank-`k` flat has a canonical independent `k`-element basis. Distinct flats have distinct bases, so there are at most `binom(8,k)` such flats.

Imposing the equations of a rank-`k` flat costs `tau^-k` before exposing the current diagonals. If the state-form rank is `r`, a nonzero incoming state requires `k<=r-1`. Return feasibility also requires `|Z|<=w-1`, hence `k<=w-1`. The total return rank factor is therefore at most

\[
\boxed{\tau^{-8}
 \sum_{k=0}^{\min(w-1,r-1)}\binom8k\tau^{-k}}.
\]

When `r` is not tracked, substitute `r=8`. For weights one and two, the extra factors are respectively `1` and `1+8/tau`. This improves the unweighted mask factors `1` and `9`. Rank-zero state is not an incoming-nonzero branch.

### Output-weight mark

On a return, `R=CX` is nonzero. The vector `AR` has weight at least nine. Thus

\[
\operatorname{wt}(Y)\ge \max(1,9-w).
\]

The extra lower bound one uses `CA=0`: if both `Y=0` and `S'=0`, then `R=eta CA R=0`, a contradiction. For `0<y<=1`, multiply the preceding rank factor by `y^max(1,9-w)`. This is a pure-return cover. Do not also multiply by a selected-output-zero enumerator on the same branch.

For `w=1`, the exact output weight is nine. Indeed,

\[
A_uA_v^T=
\begin{cases}
1,&u=v,\\
L(u+v)/(u+v),&u\ne v.
\end{cases}
\]

There are seven zero off-diagonal entries in the same `W` coset and eight nonzero entries in the other coset. The input coordinate contributes `1+eta!=0`. The one-input return envelope is consequently `y^9/tau^8`.

From zero incoming state, a nonempty return is impossible for `w<=8`. For `w>=9`, `C_I` has rank eight, the output is exactly `X`, and the pure-return bound is `y^w/tau^8`. Empty input preserves zero state; it cannot return from a nonzero state because the state diagonal is invertible.

These bounds can be combined with a separate selected-zero marginal bound on nonreturn branches. They do not prove a complete regional or whole-code estimate.

## The asymmetric alternative and its exact obstruction

Keeping the same `A`, instead set

\[
\widetilde C_{v,u}=(u+v)^7.
\]

This retains `tilde C A=0`, and the injection has full-rank eight-column Vandermonde structure. Its exact pairing is `A_u tilde C_v=(u+v)^7`. One-input returns have weight sixteen. Every two-by-two injection minor is nonzero: the seventh-power map permutes GF256*, since `gcd(7,255)=1`.

However, all-subset invertibility fails for every nonzero `eta` in GF256. An exact test finds 112 singular complement classes for `eta=2`, and 132 for `eta=7`, out of 32,768 classes each. The first bad masks use reordered phase-zero indices `(0,1,4,10,12)` and `(0,3,8,12)`, respectively. Their small Gram updates have ranks four and three.

This is a finite obstruction for all scalars, not just those two tests. Label scaling and Frobenius carry a bad scalar through

\[
\{h\eta^{2^j}:h\in H^*,\ 0\le j<8\}.
\]

The two displayed representatives give disjoint 120-element orbits that cover GF256 outside `H`. Every `eta` in `H*` already fails a two-label mask: choose difference `d` with `eta d^7=1`. The pairing depends only on label differences, so the same obstruction applies in every translated phase.

This invalidates the old invertible-update shortcut. It does not invalidate causal encoding or establish a bad codeword. The symmetric return proof above does not require all-subset invertibility.

The bounded three-input enumeration also favors the symmetric return law. Across all 560 placements in each phase, a torus combination of three asymmetric injection columns can vanish at four of the eight residues. The minimum diagonal-return exponent is then six, rather than eight. Phase maxima for the number of injection zeros are `(4,4,4,3)`. The reference enumerates every exceptional projective line using pairs of zero residue equations; the two-minor property makes this enumeration complete. Every enumerated exceptional return has output weight at least fourteen at `eta=2`.

## Fixed-map cost

Each phase is the union of two additive `W` cosets. Each map therefore uses two eight-point XOR-circulant blocks. The superset zeta transform converts each block into multiplication in the algebra with three commuting generators satisfying `e_i^2=0`. The constant transformed coefficient is one. The direct circuit has `3^3-2^3=19` nontrivial product sites.

Two injection blocks and two emission blocks can fill the four record lanes of each GFNI512 instruction. A direct packed count is:

| Operation per sixteen encoded records | GFNI512 products |
|---|---:|
| Both fixed maps | 19 |
| Scale eight emission-state records | 2 |
| Fresh eight-coordinate state diagonal | 2 |
| Sixteen input diagonals | 4 |
| Total | 27 |

With the existing 41-product outer per four GRS[32,22] rows, this is `257/88≈2.92045` products per source record. It is about 52% above `169/88`, before implementation effects.

Coefficient-basis sharing gives a lower scalar-product upper bound. For a transformed input, collect all required fixed coefficients from both output blocks. A GF2 basis of their span together with one needs one fewer products than its rank; the coefficient one is free. Transposing this straight-line circuit preserves its multiplication count.

For the symmetric pair, the emission circuit needs `(30,28,30,30)` scalar products in the four phases. Its transpose needs the same count. Their mean is 59 scalar products per packet, or 14.75 GFNI512 processed-byte equivalents. Adding the eight other packet products gives a full volume estimate of `223/88≈2.53409`, about 32% above the current four-state count. The worst-phase volume is `225/88≈2.55682`.

These sharing numbers are arithmetic upper bounds, not an optimized instruction schedule. Packing, extra XORs, transposes, coefficient loads, and doubled state storage can affect runtime. The explicit scalar reference implements the simpler nineteen-site blocks. No native kernel or benchmark was produced.

## Five-state and global scope

A five-state symmetric variant is algebraically possible: restrict `A` to any five residue columns, and transpose that restricted matrix for `C`. Cauchy superregularity and `CA=0` survive. Its emission distance is twelve, and the same pure-return proof gives five rank charges with the corresponding flat factor.

That restriction does not inherit the old five-product four-state circuit. It needs its own pruned fixed-map circuit and exact one-input return check. It is a separate cost/proof option, not a validated shortcut.

The earlier realized four-state windows impose at most 128 exit equations. Eight states increase that capacity to 256, so dimensions 136 and 132 no longer force a kernel for those fixtures. Other windows and routes remain unclosed. No failure probability follows from merely defeating those two witnesses.

Next, evaluate the symmetric eight-state flat-cover return operator on the complete unresolved sparse class. A successful local or small-active-row bound must not be promoted to a whole-code certificate. Native work should wait for a useful complete proof gate and a credible end-to-end cost case.

## Replay

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p state_capacity_reference_tests.py
python -B research/workstreams/char2_field_spin/theory/state_capacity_reference.py --all-subsets
python -B research/workstreams/char2_field_spin/theory/state_capacity_reference.py --small-returns
```

Eight tests pass. They check fixed maps, both adjoints and involutions, exact symmetric singleton weight, asymmetric two-minors, the complete scalar obstruction via its two orbits, return factors, and coefficient-rank counts. Independent GF4 tests enumerate actual live source witnesses with dependent input/state forms. A separate finite check verifies matroid closure and the weighted flat cover, including loops and repeated state forms.

## Rank-cap continuation and stronger return weights

The two-state cover forgets that a singleton birth has state-form rank at most one. It then bounds every later empty packet using the worst rank up to eight. The new [state_rank_cap_transfer.py](theory/state_rank_cap_transfer.py) preserves this information through a positive nine-state cover.

State zero means **actual zero witness state**. State `r` in `{1,...,8}` means actual nonzero witness state with form rank at most `r` in the current common source quotient. These are rank caps, not mutually exclusive exact-rank events. The cover is valid on nonnegative continuation vectors satisfying

\[
v_1\le v_2\le\cdots\le v_8.
\]

There is no required order between `v_0` and the other coordinates. In particular, the helper does not claim entrywise domination of an exact rank-transition matrix.

### One simultaneous upper-tail measure

Fix an actual incoming form rank `a<=r`. Select `z` output zeros at zero-input positions, before exposing any current input diagonal. If `H` is the remaining old-state form rank, the fresh state diagonal gives

\[
0\le H\le\min(a,8-z),\qquad
\mathbb E(256^H-1)\le\frac{256^a-1}{255^z}.
\]

The selected equations have rank `a-H`. Branches with `H=0` cannot retain an actual nonzero incoming witness and are discarded. They are not sent to actual-zero state or allowed to restart from the current input.

For each threshold `j>=1`, maximize

\[
\mathbb E[\tau^{H-a}\mathbf1_{H\ge j}]
\]

over the displayed finite-support moment constraints and over `1<=a<=r`. Call the result `W_j`, and set `W_9=0`. Each scalar moment optimization is a finite linear program. Its optimum is a single support point or a mixture of two support points. The helper evaluates these possibilities directly, including slack one-point optima.

The same sequence `W_1,...,W_8` simultaneously dominates every admissible tilted tail. It is nonincreasing in `j`. Therefore the differences

\[
m_h=W_h-W_{h+1}\ge0
\]

define one positive upper measure. Summation by parts proves that it dominates every admissible actual law on nonnegative nondecreasing functions of `H`. The measure need not be an attainable law of the diagonal refresh. Different ranks or different tail optimizers are not summed as separate histories.

### Packet update and composition contract

Suppose the packet has `w` active inputs. Select `i` output zeros at those active positions and `z` at zero-input positions. After the state-diagonal step above, each of the `i` active equations has its own fresh input-diagonal pivot, contributing `tau^-i`.

Those equations replace the corresponding inputs by linear forms in the residual old state. The following state is consequently contained in the span of that `H`-dimensional state image and at most `w-i` remaining input forms. This gives the deterministic cap

\[
r'\le\min(8,H+w-i),
\]

even if the input forms are dependent or the allowed-output update is singular. The pure-return branch is covered separately by the weighted-flat argument above.

For a nonreturn, add the positive coefficient

\[
\binom wi\binom{16-w}z
y^{16-i-z}(1-y)^{i+z}\tau^{-i}m_h
\]

to successor cap `min(8,h+w-i)`. Here `0<=i<=w` and `0<=z<=min(7,16-w)`. A birth from actual zero state goes to cap `min(w,8)`, not automatically cap eight. Actual-zero returns retain their separate rank-eight bound for `w>=9`.

Increasing the incoming cap increases every `W_j` and the weighted-flat return factor. Every packet matrix therefore maps the continuation cone into itself. Positive sums and products preserve this property. This justifies the same normalized regional matrix-polynomial products used by the two-state screen, with terminal continuation equal to one. No distribution on the rank caps is asserted.

### Improved return output weights

The new helper also improves the earlier coarse `max(1,9-w)` output bound. In the two `W`-coset ordering,

\[
AA^T=\begin{pmatrix}I&K\\K&I\end{pmatrix},
\qquad K_{u,v}=\frac1{u+v},\qquad K^2=I.
\]

The cross-coset matrix `K` is superregular. On a return the output is `(I+eta AA^T)X`. If all inputs occupy one half, its diagonal output has weight `w` and its opposite-half output has weight at least `9-w`. Thus the total is at least nine.

If both halves contain inputs, each half outside the active input coordinates has at least `max(0,9-w)` nonzero outputs. Together with the nonzero-return argument, this proves the uniform lower bounds

| Input weight | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 through 16 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Return output lower bound | 9 | 9 | 9 | 9 | 8 | 6 | 4 | 2 | 1 |

The exported `return_minimum_weight(w)` gives this table. It changes the output mark only; it does not improve the return rank charge or probability. The earlier scalar reference retains its explicitly coarser output envelope. The new rank-cap helper uses the improved table.

### Bounded checks and remaining gate

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p state_rank_cap_transfer_tests.py
```

Ten tests pass. They verify every generator of the continuation cone for all seventeen input weights, simultaneous tail inequalities, and all GF4 positive-rank subspaces and nonzero diagonals in a two-state analogue. The cap-eight nonreturn row total exactly matches the old scalar envelope; the gain comes from retaining successor caps. Further checks cover the symmetric block identity in all four phases and all nonzero two-input ratios for representative same-half and split-half placements.

The tests also keep exact Fraction evaluations separate from floating caches. They certify local algebra and comparison properties, not a whole-code distance bound. The next gate is an executed finite regional composition over the unresolved sparse class. No native implementation or runtime measurement is part of this checkpoint.

## Separate thirty-two-wide geometry

The sixteen-wide pure-return density screen has a large middle-weight gap; [STATE8_DENSITY_AUDIT.md](STATE8_DENSITY_AUDIT.md) records that executed baseline. A separate width-thirty-two construction retains the same eight residue nodes and maps, but combines the old phases in pairs:

\[
 (H+0)\cup(H+2),\qquad (H+6)\cup(H+4).
\]

Both sets are affine five-dimensional GF2 spaces, disjoint from the residue nodes. Their power sums vanish through degree thirty. Consequently the degree-fourteen products still sum to zero and `CA=0`. Cauchy superregularity now gives an emission code of parameters `[32,8,25]`. The symmetric eight-pivot return argument and weighted-flat cover are unchanged.

A one-input return has exact weight twenty-five. The Gram pairing vanishes at the seven other labels in the same `W` coset, is nonzero at the other twenty-four labels, and leaves nonzero self-output `1+eta`. For other nonempty input weights, the proved lower bound is `max(1,25-w)`.

The separate scalar reference leaves all sixteen-wide files unchanged:

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p state8_width32_reference_tests.py
```

Four tests pass: eighty basis maps, sixty-four exact singleton returns, 2048 Gram pairings, sixty-two power sums, thirty-two extended adjoints, and four stream involutions. The tests also verify the Cauchy identity at every entry and independently eliminate 128 representative square minors. The Cauchy determinant formula, not the sampled minors, proves superregularity.

## Joint selected-zero return cover

The pure-return weight minimum can lose too much at higher input weights. The new [state8_joint_return.py](theory/state8_joint_return.py) instead charges the return and selected output-zero equations in **one** numerical matrix. Its public `joint_return_bound(w,y,tau,width,cap=8)` accepts either checked width, sixteen or thirty-two. It returns the minimum of the complete pure-return and joint-return covers.

Fix an incoming actual state support `H_0`, with `h=|H_0|`, and an input support `I` of size `w`. All required coordinates are live forms in the common quotient. Select `i` output zeros on `I` and `z` output zeros outside `I`. Put `v=t-i-z`; this counts unselected output coordinates, not the actual output weight.

The numerical domain consists of the `h+w` independently scaled live state/input coordinates. Three lower bounds apply to the rank of the same joint equation matrix:

1. Return alone has rank eight, or its branch is impossible by Cauchy superregularity.
2. Inactive selected zeros require `z<h`; otherwise their Cauchy matrix forces the refreshed incoming state to zero. The output-zero rows have rank `i+z`, using the active input-identity columns and the `z` independent state rows.
3. The full map `(R,X) -> (R+CX,X+eta AR)` is invertible because `CA=0`. On a return, `R=CY` and `X=(I+eta AC)Y`. Its constrained kernel therefore embeds in the `v` unselected output coordinates. The joint rank is at least `h+w-v`.

Since `h>=z+1`, a uniform joint rank charge is

\[
\rho(w,i,z)=\max\{8,\ i+z,\ w-t+i+2z+1\}.
\]

No current diagonal is first used in a separate root probability and then reused in this charge. The alive diagonal-pivot lemma is applied once to the full joint matrix. It permits arbitrary dependencies among the source forms.

The incoming zero flat is fixed before those diagonals. Its rank `k` satisfies

\[
k\le\min(w-1,\ r-1,\ 7-z),
\]

where `r` is the incoming form-rank cap. The last bound follows from `h>z`: at most `7-z` state coordinates can belong to the zero flat. Counting canonical flat bases therefore gives the joint cover

\[
\sum_{\substack{0\le i\le w,\ 0\le z\le\min(7,t-w)\\
                 t-i-z\ge m_t(w)}}
\binom wi\binom{t-w}z
y^{t-i-z}(1-y)^{i+z}\tau^{-\rho(w,i,z)}
\sum_{k=0}^{\min(w-1,r-1,7-z)}\binom8k\tau^{-k}.
\]

Here `m_t(w)` is the proved minimum return weight for the chosen geometry. Selected-zero sets with fewer than `m_t(w)` unselected outputs are impossible. The positive selected-zero identity counts every actual return with its output mark; it does not require an exact output-support distribution.

Take the minimum with the earlier pure-return cover only after completing both sums. Both have the same successor: actual zero state. This minimum does not combine incompatible rank transitions or multiply two bounds that consumed the same randomness.

```sh
python -B -m unittest discover -s research/workstreams/char2_field_spin/theory -p state8_joint_return_tests.py
```

Six tests pass. They check every feasible return/selected-zero mask in an exact GF16 two-state analogue, the positive-mark identity, the inactive-zero feasibility test, the additional flat-rank cutoff, and comparison against the looser untyped selected-zero bound. An independent audit also checked the joint rank and diagonal chronology. These are local results; only a complete finite composition can establish the target distance guarantee.
