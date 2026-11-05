# Crossed pair layers inside the 16-by-16 outer block

This is an exact local-algebra proposal, not a distance certificate.
It changes the outer ensemble. An isolated native prototype now exists:
root reported basis and K20 dense-oracle checks for zero, one, and two
layers. Timing is pending and belongs to a separate native record.
The useful interface is a sixteen-mask transition on four coordinates.
Fresh random ratios make that transition valid even when H* outputs have
correlated amplitudes. No new intermediate diagonal is required.

The originally proposed fixed wiring also has an explicit limitation:
four-row quiet planes can remain closed under both layers. A second
wiring below removes this particular invariant set at the same pair count.
Neither fact establishes a whole-code distance bound.

## Setup and map

Fix one outer group. Its sixteen RS[16,8] rows form an array
`x=(x_(r,j))`, where `0<=r,j<16`. The current map applies independent
nonzero input diagonals D1, the fixed column mixer H*, and independent
nonzero output diagonals D2. All diagonal entries are sampled once at
setup from GF(256)* and are independent before any proved gauge reduction.

For an ordered coordinate pair, fix `c=2` in the AES field and define

\[
 M_\alpha(a,b)=(a+\alpha b,\ a+c\alpha b),
 \qquad \alpha\gets\mathbb F_{256}^{*}.
\]

Its determinant is `alpha(c+1)`, which is nonzero. The preferred proposal is

\[
 x\longmapsto D_2\,R_2\,R_1\,H_*\,D_1x.
\]

Each R layer applies M independently to every pair in a perfect matching.
Every pair, layer, and outer group receives a fresh independent alpha.
These keys remain fixed for all encoded messages. The final D2 remains
independent of all preceding randomness.

The preferred two-layer matching uses `(r xor1,j xor1)` first. In column
quartet `b=floor(j/4)`, the second partner is `(r xor d_b,j xor2)`, with
`d=(2,4,8,14)`. The constant second difference2 is retained only as a
control. The preferred one-layer candidate instead pairs `(r,j)` with
`(r xor(8+floor(j/2)),j xor1)`.

The native prototype places the smaller **column** at the first endpoint,
not the smaller integer coordinate `16r+j`. Reversing pair orientation
does not change the support law below. The algebra helper's local labels
may choose either orientation, but a value-level native oracle must use
the smaller-column convention. Each wiring above is a distinct
construction identity and must remain explicit in a future receipt.

## Exact support law and conditioning

Fix the message and all randomness preceding one pair. Its two input
values may be correlated and need not have any specified distribution.
Let `d=q-1`, where `q=256`. The remaining alpha gives this exact law:

| Input support | Output support probabilities |
|---|---|
| empty | empty with probability 1 |
| either singleton | both coordinates with probability 1 |
| both coordinates | each singleton with probability `1/d`; both with probability `(q-3)/d` |

For two nonzero values, the two cancellation keys are `a/b` and
`a/(cb)`. They are distinct and nonzero. Both outputs cannot vanish.
Because the law is identical for every representative of a support,
it remains valid after multiplying by any nonnegative weight from the past.

Disjoint pairs have conditionally independent support transitions.
The second layer uses fresh keys, so the same argument applies after
conditioning on the entire first-layer output. Intermediate amplitudes
need not be independent. The final D2 supplies independent uniform
nonzero amplitudes, conditional on the complete final support. Thus the
original inner forced-nonzero-amplitude interface remains available.

This reasoning is for fixed messages in a setup first moment. It does not
assert that a message chosen after seeing the keys has a fresh ratio.
Sharing alpha across pairs or layers gives a different law; a GF(4)
regression explicitly distinguishes that ensemble.

## Four-coordinate kernel

The two baseline matchings partition the 256 coordinates into 64 squares.
Label a square `(0,1,2,3)` so the first pairs are `(0,1),(2,3)` and the
second pairs are `(0,2),(1,3)`. For an input support `S subset {0,1,2,3}`,
define the exact positive polynomial

\[
 F_S(t_0,t_1,t_2,t_3)
 =\sum_T \Pr[\operatorname{supp}(R_2R_1x)=T\mid\operatorname{supp}x=S]
          \prod_{i\in T}t_i.
\]

The helper composes the two matching transitions using rational arithmetic.
It returns all sixteen labelled coefficients, not only their weight sum.
For arbitrary nonnegative markers and a fixed complete pre-layer support,
the superblock expectation is the product of its 64 square polynomials.
That product conditions on the actual support produced by H*; it does not
replace H* by an independently sampled support in each square.

The sharp minimum output weights for nonempty square inputs of weights
`1,2,3,4` are `4,2,2,1`. Some exact checks are

\[
 \Pr[\mathrm{wt}_{out}=2\mid S=\{0,1\}]=2/255,
 \qquad
 \Pr[\mathrm{wt}_{out}=2\mid S=\{0,2\}]=4/255^2,
\]
\[
 \Pr[\mathrm{wt}_{out}=1\mid S=\{0,1,2,3\}]=4/255^3.
\]

The difference between the two weight-two cases is why a weight-only
input state would lose geometry. A singleton deterministically expands
to all four outputs, while dense inputs can contract with small probability.

## A finite obstruction and a same-depth wiring change

For the baseline, every row set invariant under xor1 and xor2 remains
closed across all columns. In particular, `A={0,1,2,3}` is an H* quiet
affine plane: H* fixes the vector that is one on A and zero elsewhere.
If all sixteen columns have this support, the square profile is

\[
 48\text{ empty squares},\qquad16\text{ fully occupied squares}.
\]

No output can leave A. With markers in `[0,1]`, each fully occupied
square satisfies `F_full(t)>=product_i t_i`, since every possible output
support is a subset of the original square. Thus these layers cannot
improve this particular quiet-support term; cancellation can increase it.
This is not a proof that the whole construction fails, but it rejects a
universal claim that two fixed crossed directions remove quiet planes.

A distinct, same-depth schedule keeps the first matching unchanged. In
column quartet `b=floor(j/4)`, let the second row difference be

\[
 (d_0,d_1,d_2,d_3)=(2,4,8,14),
 \qquad (r,j)\leftrightarrow(r\mathbin{\mathrm{xor}}d_b,j\mathbin{\mathrm{xor}}2).
\]

Both matchings are involutions and commute inside each quartet. The same
four-coordinate kernel therefore applies. A row set closed for every
column must be invariant under xor1, xor2, xor4, and xor8. Only the empty
and full row sets have that property. More strongly, any three second
differences, together with1, span all four row bits. For the example A
above, the profile becomes 36 empty, 24 adjacent-double, and 4 full
squares. This proves a connectivity improvement, not a distance bound.

The preferred single layer has eight row differences `8,...,15`, one
per column pair. Any five of these differences span all four row bits.
It has a useful exact quiet-plane test. For `A={0,1,2,3}` occupied across
all sixteen columns, every pair has at most one active endpoint. The
64 nonzero inputs therefore become exactly 128 nonzero outputs for
every key and every nonzero input amplitude.

There are 140 affine four-planes A. Their two-dimensional direction
space contains either zero or two of the differences `8,...,15`.
Exactly 28 planes have zero aligned differences and the deterministic
output weight128 above. The other 112 planes have the exact weight PGF

\[
 t^{96}\left(\frac{2}{255}t+\frac{253}{255}t^2\right)^8.
\]

Their output weight is at least104. The eight double-input pairs can
each lose one coordinate; the other 48 active pairs expand to two.

These statements concern the full support `A times all16columns`.
For this floor, fix an active row set A that is an affine four-plane,
and require every column to be quiet before the pair layers. The matrix
`V_A^T` has rank3, so its kernel is the one-dimensional span of the
all-one vector on A. Every nonzero quiet column therefore has all four
coordinates active. Nonzero D1 entries preserve the original RS support;
each active RS row has the same L nonzero columns, with `L>=9`.
H* fixes those quiet columns. Relative to the
baseline eight outputs per active column, a fully active aligned
column pair loses at most twelve outputs. At most two pairs align, so

\[
 \operatorname{wt}_{out}\ge8L-24\ge48\qquad(L\ge9).
\]

The last inequality uses the original RS minimum weight, not an assumed
final-column activity floor. It is a deterministic lower bound for this
all-quiet four-plane class, not a bound for arbitrary nonquiet words.

## Which older interfaces survive

The final diagonal preserves the inner amplitude interface. The full
outer map is invertible, so a nonzero group remains nonzero. The unique-
inverse high-arity comparison also remains algebraically available:
condition on all post-D1 keys, invert the entire superblock, and use the
same RS counting mass. It needs a new construction binding and receipt.

The layers do **not** preserve independent output columns. Each square
lies inside a quartet of columns. For a fixed message, setup randomness
still separates across the four column quartets. A message-counting
argument must retain the RS constraints or explicitly relax them; fixing
the row set alone does not make RS column values independent. The former
sixteen-column kernel and its chronology cannot simply be relabelled.
Likewise, the RS minimum of nine nonempty input columns is not a claim
about nine final output columns. An immediate weaker consequence is at
least three nonempty output quartets. Prefer retaining original input
activity markers when building a stronger four-column transfer.

The helper exposes this distinction through `backward_activity_kernel`.
For caller-supplied positive input-support masses a_S and output weights
f(T), it returns coefficients of

\[
 \sum_S a_S\left(\sum_T\Pr[T\mid S]f(T)\right)Y^{\mathrm{flag}(S)}.
\]

The exponent is the original input-column OR flag. Flags from different
squares in the same quartet must be combined by OR before counting
columns. This is an exact local interface; it does not assume that the
input masses factor across H* squares. The nine-column RS floor can
therefore be carried without replacing it by a final-output floor.

The old columnwise gauge proof also needs replacement. The exact identity

\[
 M_\alpha(s_a a,s_b b)
 =s_a M_{\alpha s_b/s_a}(a,b)
\]

pushes nonzero coordinate scales through a pair. Conditional on preceding
scales, the rescaled fresh alpha is still uniform nonzero. Repeating this
triangular change of variables and absorbing the final scales into D2
provides the route to an equivalent normalized sampler. No native sampler
has been changed here.

## Construction and cost choices

The counts below are literal scalar circuit sites per 256-coordinate
superblock. They exclude RS, D1, H*, and the unchanged final D2 unless noted.
They are not performance measurements.

| Choice | Added pair circuit | Proof requirement |
|---|---|---|
| A: deterministic layers before D2 | 256 fixed products, 512 XORs | Exact joint value/fiber kernel after H*; no support-only Markov claim |
| B: fresh-ratio layers before D2 | 512 products, 512 XORs, 256 fresh alpha keys | Exact support kernel above; no new diagonals |
| Diagonal-refresh alternative | fresh NZ masks before each fixed layer and after the last | Exact same support transition, but extra diagonal/key work |
| C: replace H* by cheaper first mixing | Depends on the first mixer; B alone replaces 128 H* fixed-product sites by 512 pair-product sites | New RS-to-square kernel and no old outer receipt transfer |

For B, precomputed `(alpha,c alpha)` uses two variable products per pair.
Alternatively, compute `alpha*b` once and apply one fixed c multiplication.
Both use two products and two XORs. The final D2 still costs 256 products.
The one-layer candidate halves the added counts: 256 products, 256 XORs,
and 128 fresh alpha keys per superblock.
The first random ratio may be absorbed into independent D1 entries when
there is no intervening H*; that optimization does not apply through H*
without a separate algebraic argument.

The current native `coset<Group=16>` produces eight columns and calls
`mixStore` one column at a time. Crossed quartets fit inside one coset,
but four mixed columns require 64 ZMM values. A naive insertion would add
scratch traffic. One possible layout keeps the five low-rank H* correction
vectors for each of four columns and emits square coordinates on demand.
That suggestion requires code-generation and cost validation; it is not
an assertion of a free extra stage or a 3.4 ms implementation. The
isolated native prototype was developed separately from this helper.

## Next decisive test and verification

Before a full distance provider, use the sixteen exact square laws to
compare the fixed and varied wiring on actual quiet-plane supports under
the saved common marker measure. Keep each original A across all columns.
Reject the candidate if the added products/layout cannot fit the native
budget, or if the varied wiring does not remove the known quiet-plane
contribution. A successful local test would justify a positive four-column
kernel covering both quiet and nonquiet H* outputs; it would not itself
certify a whole code.

`theory/hd_outer_pair_layers.py` contains the exact laws and both wiring
identities. Its focused tests enumerate every GF(4) four-coordinate input
and every independent four-key choice. They also check final-diagonal
amplitude lifting, the shared-key counterexample, all 256 coordinates,
sharp local support floors, all140 affine-plane profiles, original
activity flags, and the trapped control plane. This helper ran no regional
provider, distance contraction, native benchmark, or remote job. Its
native status is limited to the separately reported correctness checks
at the start of this note.

The checked helper SHA256 is
`017c5ba20747ae8e25acf02393e9702a03b086038e1814ac2424cdfe49d66f60`;
the twelve-test source SHA256 is
`37cc4ce568eebaa3a9d35dffb5eb883e29f5f125c15d3b98803d35b11f1535aa`.
