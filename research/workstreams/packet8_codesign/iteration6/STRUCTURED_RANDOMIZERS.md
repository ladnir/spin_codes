# Structured 32-bit symbol randomizers

The most concrete changed-family candidate is a seven-diagonal MDS sandwich.
Its exact pointwise loss transports the existing wider24 certificate to
68.3508043660 bits. This conclusion concerns the specified random linear
family, not an unreviewed implementation or a timing prediction.

The source and tiny tests are `structured_randomizers.py` and
`test_structured_randomizers.py`. The transported endpoint receipt is
`mds_sandwich_transport_v1.json`. Frozen iteration5 sources were only read.

## Concrete family and coordinate convention

Let F be GF(256), represented as GF(16)[t]/(t^2+t+8), where GF(16) uses
x^4+x+1. The low nibble is the constant coefficient in t. The quadratic
has no root in GF(16), as the exhaustive test checks.

Use the matrix

```
H = [2 3 1 1
     1 2 3 1
     1 1 2 3
     3 1 1 2]
```

with entries in the GF(16) subfield. All 69 square minors are nonzero.
They remain nonzero over F. Thus H is invertible and every rectangular
submatrix has the largest possible rank.

At each symbol, independently sample four nonzero F-scalars for D_out
and three for D_in. Set the first coefficient of D_in to one. Define

```
R = D_out H D_in.
```

The fast transpose can apply R; the forward randomizer is then its binary
adjoint. Fresh randomness is required between symbol positions and groups.
The sampled maps are reused for every message.

The seven-coefficient family has exactly the same matrix distribution as
the family with eight independent nonzero diagonal coefficients. To see
this, divide all coefficients of D_in by its original first coefficient d,
then multiply all coefficients of D_out by d. Scalar multiplication commutes
with H, so the matrix is unchanged. The resulting seven coefficients remain
independent and uniform nonzero.

This normalization is a matrix-distribution identity. It does not omit a
whole symbol map, change the code image, or invoke message averaging.

## Fixed-input domination

Fix x!=0 with h active bytes. In the eight-diagonal description, D_in x
has independent uniform nonzero values on those h coordinates. Fix an
output support with r zero coordinates.

If r>=h, the zero-coordinate submatrix has rank h, so the specified support
has probability zero. Otherwise, choose r pivot input coordinates.
Once the other h-r nonzero coordinates are fixed, the r zero equations
determine the pivots uniquely. There are at most 255^(h-r) valid assignments.
Consequently, the probability of the specified output support is at most
255^(-r).

Conditional on the intermediate vector, D_out makes every active output
byte independently uniform nonzero. Therefore every fixed nonzero target
y has probability at most

```
255^(-r) * 255^(-(4-r)) = 255^(-4).
```

Relative to a uniform nonzero 32-bit vector, the domination factor is

```
gamma = (2^32-1)/255^4
      = 16843009/16581375.
```

The bound is tight for a one-byte input: an MDS column has four nonzero
entries, and the final diagonal makes all four outputs independently
uniform nonzero.

The binary-adjoint family has the same bound. In any fixed byte basis,
adjoint field multiplication is conjugate to ordinary multiplication by
the bytewise trace-Gram matrix. The adjoint sandwich is therefore a fixed
bytewise conjugate of a sandwich with H transposed and the diagonals
reversed. H transposed is MDS. Fixed conjugation preserves pointwise
domination by uniform nonzero vectors. No assumption identifies a field
transpose with a binary transpose.

The GF8 test enumerates all 4095 nonzero four-byte inputs for a Cauchy MDS
matrix and its transpose. Every input support satisfies the exact bound
1/7^4 after the final diagonal. A separate GF4 two-coordinate enumeration
checks the complete seven-versus-eight coefficient normalization identity
in its smaller analogue.

## Transport of the existing certificate

For each fixed outer message, every active randomized symbol contributes
at most gamma relative to its previous uniform-nonzero image law. There
are at most 16 such symbols in a group. Summing messages and applying an
independent packet permutation preserves the pointwise factor gamma^16.

The symbol maps need not be independent between messages. Linearity of
expectation is sufficient. They must be independent between positions
and groups, and independent of routing and inner setup.

For q active groups, condition on the same routing geometry used by the
iteration5 proof. The conditional first-moment bound increases by at most
gamma^(16q). Fractional clipping with its fixed witness alpha therefore
adds the factor gamma^(16q alpha). The subset union stays outside alpha.

The q1 term uses the previous exact shell argument, multiplied by gamma^16.
It does not claim that the new family retains the previous exact shell law.

The transport script authenticates the old whole receipt and corrects
each of its 256 selected dyadic endpoints. Integer exponents use positive
outward scaled arithmetic. Fractional exponents use the existing
integer-checked tangent bounds. The corrected endpoints are summed exactly.

| Quantity | Transported outward result |
|---|---:|
| Combined margin | 68.3508043660315 bits |
| q1 margin | 110.42450487393137 bits |
| Weakest occupancy | q=3 |
| Exact comparison with 2^-40 | passes |

This transport includes the old receipt's conservative 2^-200 floor on
small terms. It requires no new local census or placement calculation.

## Cost, including the fixed middle matrix

The retained tower32 multiplier costs nine GF256 products and fifteen
binary XORs per payload half before later scheduling optimizations.
The normalized sandwich uses seven arbitrary byte multiplications.
In the stated tower basis, each can be implemented by one precomputed
8-by-8 GFNI affine matrix.

This saves two variable products, but H is not free. Let
S=x0+x1+x2+x3 and define

```
u0 = 2(x0+x1)
u1 = 2(x1+x2)
u2 = 2(x2+x3)
u3 = u0+u1+u2
y_i = S+x_i+u_i.
```

Thus a GFNI implementation of H needs three, not four, fixed byte-affine
operations, plus about ten XOR/XOR3 instructions. Together with the
diagonals, this version uses ten GFNI operations per half. It may trade
some XOR work for an additional GFNI; operation counts alone do not decide
its speed.

The alternative performs the three doubled differences with integer
instructions. Multiplication by two acts independently on the nibbles:

```
t = (x >> 3) & 0x11
twice = ((x << 1) & 0xEE) xor t xor (t << 1).
```

Six ordinary SIMD operations suffice when the final three-input XOR uses
one ternary instruction. Sixteen-bit shifts are valid with these masks.
This gives roughly seven GFNI plus 28 integer instructions per half,
excluding loads and broadcasts. It reduces pressure on GFNI execution
resources but adds substantial integer work.

The hot coefficient table grows from nine bytes to 56 bytes per symbol.
At 256 groups and 16 symbols, that is 224 KiB rather than 36 KiB.
The effect of the larger table and the extra fixed work must be measured.
No timing conclusion follows from this proof.

## Two shortcuts that do not preserve the small loss

Fixing two coefficients of D_in generally fails. Choose an input supported
on those two coordinates, with their fixed ratio canceling one output row.
MDS implies the other three outputs are nonzero. The final diagonal then
gives a target atom 255^-3, a factor 255 larger than the desired envelope.
The GF8 test constructs this cancellation exactly.

The XOR-only middle H=I+J also fails the mild-domination test. For two
active equal bytes, their sum is zero and H leaves those bytes unchanged.
Independent input scales produce this event with probability 1/255.
A target with those two active output bytes consequently has mass at least
255^-3. Additional packet-support arguments might still bound this family,
but the present certificate cannot be transported with the MDS gamma.

The bounded search over all 65536 binary D in

```
D + 2J + 4uv^T,        u=v=(0,0,1,1),
```

found no matrix with every 2-by-2 minor nonzero. This rejects only that
two-product form. It is not a lower bound for all two-product circuits.

There is a simple lower bound for a more restricted implementation model.
Suppose a fixed four-byte MDS map uses only byteplane XORs and one
constant byte-linear multiplication. Restrict inputs to four binary
coefficients times a fixed byte. The multiplier's input vanishes on a
binary subspace of dimension at least three. On that subspace the complete
input/output graph is a binary code of length eight. Among seven nonzero
vectors of a three-dimensional subspace, each coordinate contributes at
most four ones. Their total weight is at most 32, so one has weight at most
four. MDS requires input-plus-output support at least five. This excludes
one-multiplier XOR-only circuits, not implementations using shifts or other
byte operations.

## Optional omission of one whole symbol map

A distinct saving is available after averaging all outer messages.
It is not a fixed-input randomizer statement.

Fix one symbol position and use the identity map there. Let A_S count
outer codewords with exact aligned-symbol support S. Global GL(8,GF16)
row transformations preserve the eight-row RS code and S. They act
transitively on any one nonzero aligned 32-bit symbol. Hence, for an
anchor position i in S and a specified nonzero value y_i, exactly

```
A_S/(2^32-1)
```

codewords have that anchor value.

Now give the other 15 symbols independent gamma-dominated maps.
For a target with i in S, summing over the compatible codewords bounds
its expected count by gamma^(|S|-1) times the fully uniform expected count.
If i is not in S, the exponent is |S|<=15. Thus the complete expected
group measure is dominated by gamma^15 times the old one.

For gamma=1, this recovers exactly the old expected group measure, even
though a fixed message no longer has a uniform image at the anchor.
No invariance of the randomizer ensemble under the row group is claimed.

Independent group setups allow products of these expected measures.
Routing conditioning is independent of outer setup, so the conditional
first moment can use this domination before fractional clipping. The
resulting factors are gamma^15 for q1 and gamma^(15q alpha) for q>=2.
This optional variant is not the family instantiated by
`mds_sandwich_transport_v1.json`, which uses all 16 maps.

The next implementation decision should compare the matched exact-map
baseline with the two MDS circuits. If one whole symbol map is omitted,
record that changed family explicitly and retain its aggregate-measure
argument. It cannot be justified by a false fixed-message uniformity claim.

The independent reviewer confirmed this omission argument and checked a
tiny GF4 MDS example exhaustively. Do not omit multiple symbols on this
argument: joint anchor labels retain row-span and proportionality invariants.

Parent-reported preliminary timings favor the exact-map scheduling variant:
roughly94microseconds versus98 for the four-fixed-GFNI H implementation and
112 for the bitwise H implementation. These are not measurements from this
workstream. The three-fixed-GFNI middle and a whole-symbol omission remain
separate bounded implementation candidates; there is no reason to broaden
the construction search before their matched tests.

The optional fifteen-map transport is reproduced by
`structured_omit_transport.py`. Its receipt uses the aggregate-measure
premise above, not a fixed-message transitivity assertion at the identity
symbol. Neither transport receipt changes the frozen baseline certificate.
