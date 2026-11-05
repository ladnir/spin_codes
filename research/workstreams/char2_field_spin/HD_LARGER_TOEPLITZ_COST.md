# Larger Toeplitz mixers: proof transfer and circuit cost

Larger invertible Toeplitz mixers retain the established uniform-measure
proof interface, but the available circuits do not predict an encoding
win. Recursive circuits spend more products than the certified g4 outer.
Evaluation-based circuits use only 15 or 31 variable products, but their
two fixed transforms erase that saving. Park this direction unless an
explicit circuit meets the thresholds below.

This is a read-only algebra and source review. No native source changed,
and no new benchmark, numerical proof search, or certificate was run.

## Candidates and proof interface

Use GF256, with q=256. For g in {4,8,16}, each 128-input group contains
g RS[n,k] rows, where n=256/g and k=128/g. Independently for every group
and column, sample a g-by-g Toeplitz matrix uniformly conditioned on
invertibility. Keep the 256 output regions, paired route with fair swaps,
and h4 inner unchanged.

The [grouped-outer argument](HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers)
gives, for every fixed nonzero x and every y,

\[
 \Pr[T x=y]\le \rho_g:=q^{-g}(1-q^{-1})^{-1}.
\]

The shortened RS information-set argument leaves n-k independent column
hits. Thus each fixed nonzero complete group output has expected
nonzero-message count at most

\[
 \rho_g^{n-k}=q^{-128}(1-q^{-1})^{-128/g}.
\]

The corresponding density relative to uniform GF256^256 is

\[
 \beta_g=q^{128}(1-q^{-1})^{-128/g}\le\beta_4
 \qquad(g=8,16).
\]

Independent group setups give the same product domination for each fixed
source-group set and route. Consequently, proof components using only this
density and the unchanged route/inner geometry need no weaker numerical
witness. This is a transfer theorem, not permission to relabel existing
receipts. A new construction binding must identify the new sampler and
connect its implementation to this law. In particular, no current g16
two-sum receipt becomes a Toeplitz receipt merely by changing metadata.

## Recursive circuits and setup traffic

The current native `toeplitzColumn` and setup support only g=2 and g=4.
The recursion T=[[A,B],[C,A]] uses three half-size products. Its product
and binary-XOR counts satisfy

\[
 P(g)=3P(g/2),\qquad X(g)=3X(g/2)+3g/2,
 \qquad P(1)=1,\ X(1)=0.
\]

The following counts apply to four groups packed in ZMM lanes. There are
256/g column calls. The RS column excludes mixer work; each product site
means one variable byte multiplication or one fixed byte-linear GFNI map.

| g | RS products | Recursive products/column | Mixer XORs/column | RS + mixer products | Prepared mixer keys at K20 |
|---:|---:|---:|---:|---:|---:|
| 4 | 516 | 9 | 15 | 1092 | 4.5 MiB |
| 8 | 392 | 27 | 57 | 1256 | 6.75 MiB |
| 16 | 272 | 81 | 195 | 1568 | 10.125 MiB |

The keys are one byte per prepared leaf, group, and column. The counts
exclude routes, inner keys, setup rejection, and padding. The zero-specialized
RS counts are 4*129, 8*49, and 16*17. An older compiled g4 report counted
1096 sites; the static 1092 count is not a correction to that disassembly.
Final-butterfly fusion and zero-coefficient handling must be audited in
any new compiled circuit.

## The 2g-1-product construction

The short bilinear circuit is real, but interpolation belongs in setup
only on the key side. Fix 2g-1 distinct field points. Let V be their
square Vandermonde matrix, and let E be its first g columns. For a vector
t of 2g-1 Hankel diagonals, prepare lambda=V^(-T)t. Then

\[
 H(t)_{i,j}=t_{i+j},\qquad
 H(t)=E^T\operatorname{diag}(\lambda)E.
\]

A Toeplitz product is H(t)Jx, where J reverses the input coordinates.
The hot circuit evaluates E, applies 2g-1 variable products, and applies
E^T. It does not evaluate a length-(2g-1) key polynomial times a length-g
input polynomial modulo an arbitrary short transform; that shortcut can
alias coefficients into the requested middle product.

One may instead use the native novel polynomial basis. If P changes
that basis to monomials, E'=EP gives the mixer P^T H(t)P. This fixed
congruence preserves invertibility and the pointwise atom bound. It is a
new mixer distribution with the same proof interface, not the original
Toeplitz matrix in the original coordinates. It avoids runtime monomial
basis conversions.

The straightforward native additive circuit evaluates degree-below-g
polynomials at 2g points and discards the last point. Two g-point coset
FFTs use 17 fixed products for g8 and 49 for g16. Removing that final
output saves an addition, not a product in this circuit. Circuit
transposition uses the same fixed-product count. Therefore:

| g | Variable products | Two fixed maps | Total mixer products/column | RS + mixer products | Prepared variable keys at K20 |
|---:|---:|---:|---:|---:|---:|
| 8 | 15 | 34 | 49 | 1960 | 3.75 MiB |
| 16 | 31 | 98 | 129 | 2336 | 3.875 MiB |

These are explicit baseline circuit counts, not lower bounds on every
possible evaluation circuit. They omit Boolean operations, register
pressure, and key expansion. The smaller setup streams are genuine, but
they do not establish that the extra fixed maps are worthwhile.

## Why RS fusion is not free

The input evaluation E' acts across constituent rows; the RS transform
acts along columns. These fixed maps commute as a tensor product. Moving
E' before RS, however, expands g rows to 2g-1 rows. The new rows must also
be evaluated by RS.

For the explicit circuit above, moving only this input map before RS
changes the static total from 1960 to 2031 products at g8, and from 2336
to 2199 at g16. The latter saves some work, but remains well above 1092.
These counts include the extra RS rows and the smaller number of input
evaluation calls; they do not assert an optimized tensor circuit.

The output map E'^T follows independent random coefficients at each
column. It cannot commute through those diagonal keys or through the RS
transform without changing their law or the emitted coordinates. Dropping
it would emit 2g-1 coordinates instead of g. A fixed basis change after
mixing also need not preserve Hamming weight. Absorbing a map therefore
requires an exact circuit identity or an explicit new construction and
proof interface, not a change of notation.

Nor does a GFNI affine instruction combine distinct 128-bit records. It
implements a binary-linear map inside each byte. Cross-record evaluation
and transpose sums remain real operations.

## History and reconsideration threshold

The large-g negatives in [HD_GROUPED_SUBFIELD.md](HD_GROUPED_SUBFIELD.md)
used diagonal-Cauchy-diagonal mixers, not g8/g16 Toeplitz. The older
[nineteen-product convolution](GF8_FULLSTATE_COST.md#nineteen-product-convolution)
and [coefficient-sharing outer circuits](LARGE_OUTER_COST.md#outer-circuit-and-checks)
concern fixed XOR-circulant kernels. They are not implementations of this
random Toeplitz bilinear circuit. No saved g8/g16 grouped bilinear
Toeplitz implementation was found in the inspected native files or notes.

Let C be the additional fixed-product cost per column after the 15 or 31
variable products, with the original RS schedule retained. Beating the
static g4 product total requires

\[
 392+32(15+C)<1092\quad\Longleftrightarrow\quad C<6.875
\]

at g8, or

\[
 272+16(31+C)<1092\quad\Longleftrightarrow\quad C<20.25
\]

at g16. The current normalized two-sum mixers already use the same 15/31
variable products and only 6/8 fixed products per column. Matching their
arithmetic therefore requires C<=6/8 before considering XORs and layout.

The next useful gate is an exact scalar factorization with a counted,
register-bounded schedule below one of these thresholds. It should test
all bilinear basis pairs and the setup-law conversion before native work.
Without that circuit, reconsider the faster two-sum construction only
with a different correlation-preserving proof. The latest
[exact-constraint controls](HD_TWO_SUM_GF16_GROUP_SEEN.md) do not justify
another search over the same relaxed bounds. Its current
[g16 store-fused measurement](HD_TWO_SUM_GF16_NATIVE.md#matched-store-fusion-result)
is 3.397943 ms; the initial outlined-store regression is not the current
candidate. Nothing in this note predicts a new runtime or replaces that
measurement.
