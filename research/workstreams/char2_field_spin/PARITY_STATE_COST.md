# Extra outer parity versus a larger regional state

Several inexpensive algebraic changes merit proof screening.
GRS[33,24] adds seven fixed products to the verified GRS[32,24] outer transpose.
GRS[32,23] instead removes one source dimension and adds only one fixed product.
Both increase the outer minimum distance from nine to ten.
Removing two source dimensions gives GRS[32,22], also with seven added products, and increases the outer minimum distance to eleven.
These arithmetic results do not establish a distance certificate for any complete encoder.
The earlier GRS[32,24], four-state profile remains rejected at K20 by the whole-subspace obstruction recorded in [RATE_DESIGN_REVIEW.md](RATE_DESIGN_REVIEW.md).

## The extra parity at point 32

Use AES GF256 and the systematic interpolation points 0 through 23.
The new parity coordinate is \(f(32)\), where \(\deg f<24\).
Independent Lagrange interpolation gives its coefficient row

```text
225,153,25,168,71,121,240,55,199,48,100,145,
128,251,212,207,151,144,69,220,104,24,113,204.
```

The eight constants

```text
1,225,153,25,168,71,240,100
```

form a GF2 basis of GF256.
For the single new transpose input \(x\), compute the seven nontrivial basis multiples of \(x\).
The reference then derives all 24 required multiples with 23 binary XORs.
Its complete fixed program is verified for every byte value of \(x\).
Thus the parity kernel grows from 34 to 41 GFNI512 products per four source rows.

The final store can combine the systematic value, original parity contribution, and new contribution with one three-input XOR.
The complete outer arithmetic therefore uses 41 products and 161 Boolean instructions, excluding output lane shuffles.
The original GRS[32,24] outer uses 34 products and 138 Boolean instructions under the same convention.

Adding points 33 through 35 reuses this scalar circuit with permuted output indices.
For offset \(h\in\{0,\ldots,7\}\), the coefficient at source index \(j\) for point \(32+h\) is the displayed coefficient at \(j\mathbin{\mathrm{xor}}h\).
A direct extension gives these verified upper counts:

| Outer | Rate | Outer minimum distance | Products per four source rows | Boolean instructions, including systematic addition |
|---|---:|---:|---:|---:|
| GRS[32,24] | 0.7500 | 9 | 34 | 138 |
| GRS[33,24] | 0.7273 | 10 | 41 | 161 |
| GRS[34,24] | 0.7059 | 11 | 48 | 208 |
| GRS[32,22] | 0.6875 | 11 | 41 | 170 |
| GRS[36,24] | 0.6667 | 13 | 62 | 278 |

The counts for multiple added coordinates do not claim optimal sharing between their inputs.

## Removing one source dimension costs one product

To obtain GRS[32,23], interpolate source values at points 0 through 22.
First compute the value at point 23, then apply the existing GRS[32,24] map.
The interpolation coefficients for that missing value are

\[
(\underbrace{12,\ldots,12}_{8},
  \underbrace{13,\ldots,13}_{8},
  \underbrace{1,\ldots,1}_{7}).
\]

For the transpose, let \(p\) be the 24 parity contributions produced by the existing circuit from encoded coordinates 24 through 31.
Let \(y_0,\ldots,y_{31}\) be the complete transpose input and set

\[
t=y_{23}+p_{23},\qquad u=12t,\qquad v=u+t=13t.
\]

For output indices 0 through 7, return \(y_j+p_j+u\).
For indices 8 through 15, return \(y_j+p_j+v\).
For indices 16 through 22, return \(y_j+p_j+t\).
This complete outer transpose uses 35 products and 139 Boolean instructions per four source rows.
The 23-symbol row stride requires a different output-store schedule; no claim about its shuffle or cache-line cost is made here.

## Removing two source dimensions costs seven products

For GRS[32,22], interpolate at points 0 through 21 and reconstruct the values at points 22 and 23 before applying the existing GRS[32,24] map.
All constants in this subsection denote bytes in AES GF256; addition is XOR.
Define

\[
c_j=\begin{cases}12&0\le j<8,\\13&8\le j<16,\\1&16\le j<22.\end{cases}
\]

The two missing-value rows are

\[
f(22)=\sum_{j=0}^{21}c_j(j+23)f(j),\qquad
f(23)=\sum_{j=0}^{21}c_j(j+22)f(j).
\]

These follow from the two interpolation parity checks
\(\sum_{j=0}^{23}f(j)/Q'(j)=0\) and
\(\sum_{j=0}^{23}jf(j)/Q'(j)=0\), where
\(Q(X)=\prod_{j=0}^{23}(X+j)\) and \(\deg f<22\).
The two missing points have equal derivative values, and
\(Q'(23)/Q'(j)=c_j\) on the remaining points.
The reference independently checks both rows by direct Lagrange interpolation.

Let \(p\) again be the existing circuit's 24 parity contributions, and put

\[
t_{22}=y_{22}+p_{22},\quad t_{23}=y_{23}+p_{23},\quad
s=t_{22}+t_{23},\quad w=t_{22}+22s.
\]

The missing-value transpose contributes \(c_j(w+js)\) to output \(j\).
Compute the six products \(2s,4s,8s,16s,32s,64s\), then one product \(12w\).
All other required multiples follow by XOR, including \(13w=12w+w\).
The 22 corrections form three affine cubes, truncated to six entries in the last group:

| Output group | Base | Three directions | Number of outputs |
|---|---|---|---:|
| 0 through 7 | \(12w\) | \(12s,24s,48s\) | 8 |
| 8 through 15 | \(13w+104s\) | \(13s,26s,52s\) | 8 |
| 16 through 21 | \(w+16s\) | \(s,2s,4s\) | 6 |

Here \(12s=4s+8s\), \(24s=8s+16s\), \(48s=16s+32s\);
add \(s,2s,4s\) respectively to obtain \(13s,26s,52s\).
Finally \(104s=8s+32s+64s\), and \(22s=2s+4s+16s\).
Generating each cube by adding one direction to an earlier entry takes one binary XOR per non-base output.
The complete correction uses seven products and 32 Boolean instructions, counting a three-input XOR as one instruction.

The final outputs are \(y_j+p_j+c_j(w+js)\), for \(0\le j<22\).
Including the existing parity circuit, the two missing-input XORs, and 22 final three-input XORs gives **41 products and 170 Boolean instructions per four source rows**.
The correction's coefficient pairs, together with the two free input forms, have GF2 dimension nine.
An XOR/fixed-GF256-multiplication circuit starts with dimension two, and each product introduces at most one new signal dimension.
Thus seven products are necessary for this isolated correction in that circuit model; this is not an optimality claim for the complete outer transpose or for other bitwise instruction models.

Compared with GRS[34,24], the outer has the same minimum distance eleven, fewer products per four rows, but fewer source symbols per row.
The full state-four arithmetic volumes are nearly equal: \(177/88\) for GRS[32,22] versus \(385/192\) for GRS[34,24], a 0.307% increase for the shortened outer.
Their rates are respectively \(11/16\) and \(12/17\).
The 22-symbol output stride is 352 bytes for 128-bit records, rather than 384 bytes for 24-symbol rows; packing and store costs remain unmeasured.

## Five and six state symbols

For state size \(s\in\{4,5,6\}\), use the actual embedded GF16 labels \(u\) and

\[
v_s(u)=(1,u,\ldots,u^{s-1})^T,
\quad C_s=(v_s(u))_{u\in\mathrm{GF}(16)},
\quad A_s=C_s^T.
\]

Every entry of \(C_sA_s\) is a power sum of degree at most ten, so \(C_sA_s=0\).
The nonzero emission code has minimum weight \(17-s\): increasing the state strengthens some rank constraints but reduces this weight guarantee.

The translation matrix has entries

\[
(T_{\beta,s})_{i,j}
=\binom{i}{j}\bmod2\;\beta^{i-j}\quad(j\le i),
\]

with zero entries above the diagonal.
Its Boolean-cube shear factorization costs \(\sum_{i=0}^{s-1}\operatorname{popcount}(i)\) fixed products and the same number of XORs.
This gives four, five, and seven products for states four, five, and six.
Reversing and transposing the shear factors gives the adjoint with the same cost.
The phase offsets and once-per-region orientations remain as in [REGION64_CANDIDATE.md](REGION64_CANDIDATE.md).

The existing additive basis

\[
1,X,X^2+X,X^3+X^2,S_2(X),XS_2(X),
\quad S_2(X)=X^4+7X^2+6X
\]

uses abstract GF16 constants, embedded before arithmetic in AES GF256.
The reference verifies the monomial-to-additive conversion for each state size.
Truncating the six-state factorization in `SubfieldWide.h` gives the following per-step GFNI counts:

| State | Both emission/moment maps | Input diagonal | Internal phase change | Nonzero emission minimum |
|---:|---|---|---|---:|
| 4 | 12 ZMM | 4 ZMM | 4 XMM | 13 |
| 5 | 14 ZMM + 2 XMM | 4 ZMM | 5 XMM | 12 |
| 6 | 16 ZMM + 4 XMM | 4 ZMM | 7 XMM | 11 |

Here one step consumes sixteen 128-bit records.
The XMM products in the maps convert between monomial and additive bases.
Phase changes also need four, five, or seven XOR128 instructions.
No random state matrix is multiplied at an internal phase boundary.

## Full normalization by source size

Let the outer have length \(b\), dimension \(k\), and let \(K=kL\) be the total source size.
The encoded length is \(N=bL\), so there are \(bK/(16k)\) inner steps.
For each profile, include both inner maps, input diagonal, internal phase changes, and the complete outer transpose.
The leading GFNI counts per source symbol are:

| Outer, state | GFNI512/source | GFNI128/source | Processed-byte equivalent | Minimum record traffic, bytes/source |
|---|---:|---:|---:|---:|
| 32/24, 4; rejected reference | 1.68750 | 0.33333 | 1.77083 | 80.000 |
| 33/24, 4 | 1.80208 | 0.34375 | 1.88802 | 82.000 |
| 32/23, 4 | 1.77174 | 0.34783 | 1.85870 | 82.783 |
| 34/24, 4 | 1.91667 | 0.35417 | 2.00521 | 84.000 |
| 32/22, 4 | 1.92045 | 0.36364 | 2.01136 | 85.818 |
| 36/24, 4 | 2.14583 | 0.37500 | 2.23958 | 88.000 |
| 32/24, 5 | 1.85417 | 0.58333 | 2.00000 | 80.000 |
| 32/24, 6 | 2.02083 | 0.91667 | 2.25000 | 80.000 |

The processed-byte equivalent is GFNI512 count plus one quarter of the GFNI128 count.
It measures arithmetic volume; it is not an instruction-throughput model.
For a given platform, the two vector-width counts should remain separate.
The table does not hide the extra inner work caused by a longer outer code.

The exact correction for phase products is to subtract \(b\) phase changes: each region's first step uses its prepared entry matrix.
For a zero initial state, the \(b-1\) remaining dense regional boundary maps cost:

- state four: four GFNI512 products per boundary;
- state five: five GFNI512 and five GFNI128 products per boundary;
- state six: six GFNI512 and six GFNI256 products per boundary.

Add these boundary counts and divide by \(K\) when evaluating a concrete size.
The final state is discarded, so there is no final matrix conversion.
Regional orientation coefficients are prepared outside encoding.

Minimum record traffic is \(16(3N+K)\) bytes: inner input read, plane write, plane read, and output write.
The regional temporary-buffer schedule adds \(32N\) bytes.
Routing indices and input diagonal bytes add their own traffic; neither changes the comparison between state sizes at a fixed outer.
For GRS[32,22], minimum traffic is \(944/11\) bytes per source symbol; the regional temporary-buffer schedule raises it to \(1456/11\).
Cache write allocation, lane shuffles, register spills, and the shortened output strides remain unmeasured.
Concrete power-of-two payloads must also account for padding to a multiple of \(64k\).

An extra parity coordinate is therefore a smaller arithmetic change than moving to state five.
GRS[32,23] is smaller still in arithmetic volume, but has slightly more record traffic than GRS[33,24].
Which change repairs the whole-code obstruction is a separate proof question.

## First budget for a fully superregular packet map

A separate proposal replaces the packet map with a Cauchy matrix on the input and state coordinates, with a fresh state refresh per packet.
For width eight and two state symbols, a straightforward packed implementation costs:

- 12 GFNI512 products for the eight-by-eight XOR-circulant input block, using three four-by-four block products;
- four GFNI512 products for each of the two cross maps;
- two GFNI512 products for the input diagonal;
- two GFNI256 products for the state block and two more for the random GL2 refresh.

This totals 22 GFNI512 and four GFNI256 products per eight outputs.
With GRS[33,24], the processed-byte equivalent including the outer is approximately 4.55208 per source symbol.
This is an upper count for that elementary schedule, not a lower bound on the construction.
A faster joint Cauchy circuit would need to improve that budget before native work is justified by performance.

## Exact checks and next gate

The reference verifies all 256 scalar inputs for each of the eight extra-parity row permutations.
It also checks combined GRS[33,24], GRS[34,24], and GRS[36,24] transpose inputs, and the complete shortened GRS[32,23] transpose.
For GRS[32,22], it checks all 65,536 pairs of missing-value transpose inputs against independent interpolation, 64 complete combined-input transposes, the seven-product correction count, the complete 41-product/170-Boolean count, and the GF2 dimension-nine certificate.
For state sizes four, five, and six, it checks all field translation offsets against direct binomial matrices, their adjoints, \(T_\beta^2=I\), composition, \(C_sA_s=0\), and the additive-basis conversion.
All checks pass.

```sh
python -B research/workstreams/char2_field_spin/theory/parity_state_cost_reference.py --check
```

The next gate is a positive whole-code proof for one of these changed outer/state choices.
No new native kernel, remote run, benchmark, or complete distance certificate accompanies this study.
