# A parity mixer followed by one crossed pair layer

This is a new outer proposal, not a transferred distance certificate.
The first mixer is `H0=I+11^T` on sixteen field coordinates. It uses
31 XORs and no field products per column. A subsequent crossed pair
layer can expand its two-row quiet inputs. Exact local counts are now
available without enumerating field assignments.

The proposal retains the current RS[16,8] rows, outer groups, paired
route, and original one-state inner map. None of the helpers in this note
runs a regional provider, distance screen, native kernel, or benchmark.
The later one-product pair variant has a different sampler and bound.

## Ensemble and original support

Fix sixteen RS[16,8] rows and their sixteen output columns. Let A denote
the set of nonzero RS rows; A remains fixed across the whole word.
Sample independent uniform nonzero D1 entries at setup, and apply H0 separately
to each column. Then pair `(r,j)` with

\[
 (r\mathbin{\mathrm{xor}}(8+\lfloor j/2\rfloor),\ j\mathbin{\mathrm{xor}}1).
\]

For the first candidate, each pair receives its own independent
`alpha` from GF(256)* and applies

\[
 (a,b)\longmapsto(a+\alpha b,\ a+2\alpha b).
\]

The smaller column is the first endpoint. An independent final uniform-nonzero
diagonal D2 follows the pair layer. All setup keys remain fixed for all
messages. The pair layer's support-only law holds for every fixed input
value pair; final D2 supplies the original inner amplitude interface.
See [the pair-layer note](HD_OUTER_PAIR_LAYERS.md) for that law and its
conditioning proof.

H0 is a symmetric involution: in characteristic two, `(11^T)^2=0`
because the length is even. Consequently, one column is zero after H0
exactly when its input is zero. The pair layer then couples two columns.
The RS floor is still nine nonzero **original** columns, not nine final
columns. A two-column transfer must mark whether zero, one, or both
original columns are nonzero.

## Two distinct counting interfaces

For one column, write `s=sum_i x_i` and `y=H0 x`. Then

\[
 y_i=x_i+s\quad(i\in A),\qquad y_i=s\quad(i\notin A).
\]

The helper distinguishes two domains:

- `allowed`: `x in F^A`, including zeros at coordinates in A;
- `forced`: every coordinate in A is nonzero, and every other coordinate is zero.

The first domain is the ambient counting measure used in an RS counting
comparison. It has mass `q^a`, where `a=|A|` and `q=256`. The second
domain has mass `(q-1)^a`. Dividing its counts by that mass gives the law
after independent uniform-nonzero D1 on a fixed exact support. These
normalizations are not interchangeable.

Let `d=q-1`. Define the numbers of nonzero b-tuples with sum zero or with
sum equal to one fixed nonzero field value:

\[
 Z_b=\frac{d^b+d(-1)^b}{q},\qquad
 N_b=\frac{d^b-(-1)^b}{q}.
\]

These formulas include `Z_0=1` and `N_0=0`. For a labelled inside output
set `B subset A`, put `b=|B|` and `O=A^c`.

### Allowed inputs

The quiet branch `s=0` has output support B and count `Z_b`.
In the nonquiet branch, divide all values by s. The inside values
`z_i=y_i/s` satisfy

\[
 \sum_{i\in A}z_i=1+(a\bmod2).
\]

Thus the output support is `O union B`, with count

\[
 (q-1)\begin{cases}Z_b,&a\text{ odd},\\N_b,&a\text{ even}.\end{cases}
\]

There is no factor `choose(a,b)` for a fixed labelled B. Such a factor
appears only when explicitly summing all B of that weight. When O is
empty, quiet and nonquiet branches may share an output support. Their
counts must be added. The full result has exact zero coefficient one.

### Forced nonzero inputs

The quiet branch has only the support A, with count `Z_a`. For `s!=0`,
the normalized inside outputs satisfy `z_i!=1`. At their nonzero
positions, they therefore lie outside `{0,1}`. For `b>0`, the numbers
summing to zero and one are equal:

\[
 U_b=\frac{(q-2)^b+(q/2-1)(-2)^b}{q}.
\]

The nonquiet support `O union B` has count `(q-1)U_b`. At `b=0`, replace
`U_0` by one for odd a and zero for even a. This separate empty-tuple
case is necessary. It gives, for example, the singleton law: a nonzero
input at i produces exactly all fifteen coordinates other than i.

For two forced nonzero inputs, the complete law is particularly small:

\[
 \Pr[\operatorname{supp}y=A]=1/255,\qquad
 \Pr[\operatorname{supp}y=\{0,\ldots,15\}]=254/255.
\]

The first event is equality of the two input values. It is a genuine
quiet branch; the cheap mixer does not remove it.

## Does the crossed layer repair two-row quiet words?

Fix `A={r,s}` and let `Delta=r xor s`. There are 56 row pairs with
`Delta<8` and 64 with `Delta>=8`. In the first class, none of the eight
crossing differences aligns with A. In the second, exactly one aligns.

Suppose first that every original column is quiet and nonzero. Its two
active values are equal, but the crossed layer's law does not assume
anything about values in different columns. For `Delta<8`, all 32
active coordinates meet zeros, so the output weight is exactly64.
For `Delta>=8`, 28 active singletons expand, while two pairs have two
active endpoints. The exact output-weight polynomial is

\[
 t^{56}\left(\frac{2}{255}t+\frac{253}{255}t^2\right)^2.
\]

Its minimum output weight is58. These are exact full-sixteen-column
examples, not assumptions that every RS word has full support.

For an all-quiet two-row RS word with L nonzero columns, both rows have
the same support. Thus `L>=9`. Relative to four outputs per active
column, one fully active aligned column pair can lose at most six
outputs. Hence the deterministic floors are

\[
 \operatorname{wt}_{out}\ge
 \begin{cases}4L\ge36,&\Delta<8,\\4L-6\ge30,&\Delta\ge8.\end{cases}
\]

The old quiet branch had only `2L` outputs. This is a concrete gain on
the low-arity obstruction, but neither floor reaches the full distance
target by itself. Nonquiet inputs and inner-state paths remain necessary.

## A small classification for scalar two-column kernels

H0 commutes with every row permutation. For one difference delta, let

\[
 f_\delta(A)=\tfrac12|A\cap(A\mathbin{\mathrm{xor}}\delta)|,
 \qquad s_\delta=|A|-2f_\delta.
\]

Permuting the eight matched row pairs, and swapping their endpoints,
maps any A to any other set with the same numbers of full and singleton
pairs. These permutations commute with the crossing. Therefore a scalar
equal-marker two-column kernel, including its original-activity marker,
depends only on `(a,f_delta)`.

There are 44 nonempty pairs `(f,s)` satisfying `f+s<=8`. Thus one common
scalar marker requires at most 44 distinct two-column kernels. The
eight differences still see one fixed A. The helper retains the ordered
vector `(f_8,...,f_15)`; it may retain only its histogram when the same
scalar pair functional and markers apply to all eight column pairs.
The resulting scalar polynomial factors commute. This does not justify
collapsing noncommuting chronological matrix operators.

Splitting A into lower and upper eight-row masks makes the profile
especially simple:

\[
 f_{8+k}(A)=|A_{low}\cap(A_{high}\mathbin{\mathrm{xor}}k)|.
\]

It is an eight-point XOR cross-correlation. Exact class multiplicities
can therefore be counted with bit operations rather than field ranks.
They multiply once after the fixed-A word calculation. For arity two,
the entire scalar class needs only representatives `(0,1)` and `(0,8)`,
with multiplicities56 and64. Arity one needs one representative with
multiplicity16.

## A one-product symmetric pair alternative

A cheaper, distinct pair sampler chooses alpha outside `{0,1}` and sets

\[
 t=\alpha(a+b),\qquad (y_0,y_1)=(a+t,b+t).
\]

Its matrix is `I+alpha*11^T` in dimension two. It has determinant one,
is an involution, and has no zero entries. The circuit uses one variable
field product and three XORs. It halves the added field products versus
the fresh-ratio pair, at the cost of one additional XOR per pair.

Its exact double-input law needs an equality flag. If `a=b!=0`, both
outputs remain nonzero. If `a!=b` and both are nonzero, each singleton
has probability `1/(q-2)`, and the double output has probability
`(q-4)/(q-2)`. Singletons always expand to two outputs.

For product markers `0<=u,v<=1`, the unequal-input expectation dominates
the equal-input expectation because their difference is

\[
 \frac{u+v-2uv}{q-2}\ge0.
\]

Consequently, the unequal-input law provides a support-only pointwise
majorant for **one** layer evaluated at such markers. It can be applied
before a common positive atomic functional; all its atoms must satisfy
the marker restriction. It is not coefficientwise domination of a
probability law, and it does not establish a multi-layer Markov law.
For arbitrary nonnegative output functions, a safe positive envelope
instead assigns mass one to the double output and `1/(q-2)` to each
singleton. Its surcharge `1+2/(q-2)` must not be discarded.

Independent final D2 still supplies the amplitude interface. However,
the earlier input-gauge pushforward is not closed in this symmetric
matrix family. The isolated GF(4) test exhibits that failure. Use full
unnormalized D1 unless a new gauge proof is supplied. This adds sixteen
input products per superblock relative to the current column-normalized
sampler, against the 128 products saved by halving one pair layer.
No native cost or distance claim follows from these static counts.

## A binary rank-five middle mixer

The parity mixer has quiet vectors of weight two. A separate candidate
retains five binary checks without the eight fixed field products of H*.
Index rows by four-bit labels, and let V have row `(1,i_0,i_1,i_2,i_3)`
at label i. Define

\[
 H_{bin}=I+VV^T=I+11^T+\sum_{k=0}^3 b_kb_k^T.
\]

Every entry of `V^T V` is zero: each relevant binary monomial has even
weight on the sixteen labels. Thus `H_bin^2=I`. V has rank five, so
`H_bin x=x` exactly when `V^T x=0`. Any three distinct affine columns
are independent, whereas four labels in an affine plane are dependent.
The minimum nonzero quiet weight is therefore four over GF(256).

The isolated circuit computes the five checks with 26 XORs, evaluates
their corrections with 15 XORs, and adds the input with 16 XORs. Its
total is 57 XORs and no field products per column, excluding diagonals
and the crossed pair layer. This saves nine XORs and eight fixed field
products relative to the H* circuit. It is a static circuit comparison,
not a native timing result.

The cheaper correction is not an MDS mixer. Its singleton output
weights are seven for seven labels, nine for eight labels, and fifteen
for one label. The new pair layer must therefore handle those labelled
supports. Neither the H* full-support singleton law nor its receipts
transfer.

For either binary mixer, two allowed input rows have a small joint
counting interface. Fix distinct rows a and b, and let C_a,C_b denote
the supports of their matrix columns. The exact branches are

| Original input weight | Output support | Count |
| --- | --- | ---: |
| 0 | empty | 1 |
| 1 | C_a | q-1 |
| 1 | C_b | q-1 |
| 2 | C_a symmetric-difference C_b | q-1 |
| 2 | C_a union C_b | (q-1)(q-2) |

The last two branches distinguish equal and distinct nonzero input
values. Binary coefficients make every shared coordinate equal to their
sum. Equal dictionary keys are added; the total is `q^2`.

`hd_binary_syndrome_outer.two_row_joint_counts(a,b,mixer=...)` exposes
these exact counts for `parity` and `binary_syndrome`. No RS factor is
included. Summing the first coordinate over all sixteen original
columns retains the aggregate input-weight floor18. This is distinct
from retaining at least nine active original columns. Both constraints
are valid and can be kept in one positive polynomial. The fixed pair
`{a,b}` remains unchanged across all columns; no H_bin class reduction
is asserted.

## Integration and verification status

`hd_parity_outer.support_counts(A, mode='allowed')` returns sparse exact
`outputmask -> integercount` data. The companion `support_branches`
keeps the parity-sum flag. `weight_coefficients` exposes the compact
per-labelled-subset formulas without expanding all masks. These APIs
are suitable for the generic pair contraction after optional densifying.

Use allowed counts for the ambient RS comparison; use forced counts
only when the exact support and its input-diagonal normalization are
explicit. The whole-group RS factor is `255^(-8a)` once, not once per
column. Input activity is marked before the crossed pair layer.
No old H* outer kernel, sampler receipt, or full-code certificate transfers.

Seven parity tests enumerate every GF(4) input for every support in four
coordinates, including both counting modes. They also check branch
overlap, involution, the two complete quiet-pair classes, and their
activity floors. Five symmetric-pair tests check its exact law, marker
majorant, invalid above-one use, and failed old-gauge transfer.
Five additional binary-mixer tests verify the circuit, involution,
quiet minimum, and all 120 row pairs for both binary mixers by exhaustive
GF(4) enumeration. The combined 17-test suite passes.

The next bounded proof test should retain all120 labelled arity-two
sets in one common positive state-aware measure, with original activity
and aggregate input weight retained. The two-representative reduction is **scalar only**; no
noncommutative orbit reduction has been proved. The first scalar control
reported arity one at -4753.45 bits and arity two at +4579.81 bits across
its prescribed slopes. These are proposal diagnostics, not certificates.

An unrestricted all-arity scalar screen is not the recommended next
step. At the current bank its positive full-arity branch alone contributes
+4661.53 bits. The 44-type classification remains mathematically useful,
but does not remove that state-collapse loss. Prefer state-aware low-
arity bounds and the portable unique-inverse high-arity comparison.
Any high-arity receipt must explicitly bind the new invertible mixer and
sampler; old receipts cannot be relabelled. None of these failed scalar
majorants proves that the construction itself fails.

## State-aware low-arity controls

The next bounded control retains every labelled singleton and two-row set
through all eight actual crossed blocks. The two-column tensor kernel
keeps both chronological inner-state chains. Its XOR contractions were
independently checked against direct state-chain sums on small instances.
All branches use the same positive five-atom measure, including its empty
entries. The selected source class is still r=18, Fmax<=1, at mark9927/10000.
These are FLOAT diagnostics, not directed or whole-code certificates.

The initial finite input-marker screen gave misleadingly poor singleton
results. Its numerator was dominated by the all-zero word: each doubling
of the marker reduced the reported bound by exactly162 bits. The scalar
comparison had instead removed low-activity coefficients exactly. Thus
the initial discrepancy was not evidence against the new state kernel.
The initial receipt is retained as `theory/hd_pair_state_H0_9927_float.json`.

The replacement uses a positive polynomial in original input weight.
For each pair block and fixed boundary states, form
`P(X)=sum_b n_b J_b^(1/r) X^(w_b)`, where n_b is the exact branch count,
w_b is its original input weight, and J_b is its r-th moment matrix entry.
Expand P(X)^r entrywise, multiply the eight polynomial matrices in order,
and retain only degrees at least9ar. Generalized Holder bounds each
coefficient separately. The aggregate floor relaxes the individual group
constraints, so it remains an upper bound. All labelled-A roots are summed
only after the whole-word tail and RS normalization.

The calculation uses positive log-domain convolution, without coefficient
subtraction or signed FFTs. It omits the separate nine-active-column
constraint; for arity one that constraint equals the input-weight floor.
For arity two, its omission is a further relaxation. The saved common
measure receives a small FLOAT log-scale lift of2^-36 to cover its recorded
rounding deficit; this does not make the result a directed enclosure.

| First mixer | Pair law | Singleton bound log2 | Two-row bound log2 |
|---|---|---:|---:|
| H0 | Exact fresh ratio | -3991.9841 | +3334.6729 |
| H0 | One-product symmetric majorant | -3991.9797 | +3334.6729 |
| Hbin | Exact fresh ratio | +926.3838 | +104.3880 |
| Hbin | One-product symmetric majorant | +926.3979 | +104.4213 |

The two controls took approximately12 seconds each, sequentially. Their
receipts are `theory/hd_pair_state_tail_H0_9927_float.json` and
`theory/hd_pair_state_tail_Hbin_9927_float.json`. A positive entry does not
prove a low-distance word exists; it means this selected-class bound fails.
H0 retains good singleton spreading but gives an inadequate two-row bound.
Hbin improves that bound considerably but gives an inadequate singleton
bound. Neither proposal is ready for a native performance campaign.

The reusable implementation is `theory/hd_pair_state_tail.py`, called by
`theory/hd_pair_state_tail_control.py`. Independent review checked the
coefficientwise Holder argument, the9ar floor, heterogeneous-group root
sum, and the controller's factors. The combined local suites pass39 tests:
31 for the state kernel, controls, and local mixers, plus8 polynomial-tail
tests. They include noncommuting state order, exact small rational
comparisons, and preservation of a nonzero coefficient at -18000 bits.

Next, separate insufficient singleton spreading from unnecessary bound
loss before adding more stages. Hbin with varied column wiring is one
unmeasured construction option; the current H* with a cheaper pair stage
is another. Do not tune the discarded finite-marker screen or count these
results as progress toward a completed35%/40-bit certificate. No native,
production, or paper files changed, and the certified fallback is intact.
