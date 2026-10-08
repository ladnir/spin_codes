# Width-eight GF256 packet maps

The preferred local candidate has eight input records and four state records.
It repairs the existing GF4-circulant construction without adding multiplication sites.
Both coordinate maps are fully superregular, and the exact return minima are

\[
(8,6,5,4,3,2,2,1)
\]

for input weights one through eight.
The packed arithmetic schedule uses nine GFNI512 instructions per eight 128-bit records.
This note establishes local algebra and an instruction schedule, not a whole-code distance certificate or measured speedup.
The separate [whole-code checkpoint](GF8_PACKET_CERTIFICATE.md) composes these maps with GRS[128,83] at K20 and proves a 65-bit setup-failure margin at 10% distance.

The reproducible references are
[`gf8_packet_maps.py`](theory/gf8_packet_maps.py) and
[`gf8_packet_maps_tests.py`](theory/gf8_packet_maps_tests.py).
The profile identifier is
`gf256-repaired-cubic-w8-s4-eta2-c3-regional-nonzero-diagonals`.
All existing construction, proof, and native files remain unchanged.

## Packet interface and exact constants

Arithmetic uses the AES field \(F=\operatorname{GF}(256)\), with polynomial `0x11b`.
Integer constants below denote field elements, and addition is XOR.
Each record contains sixteen parallel field elements; every lane uses the same scalar map.

For input \(u\in F^8\) and incoming state \(s\in F^4\), define

\[
x=D_xu,\qquad r=D_ss,\qquad
z=x+2Ar,\qquad t=r+Cx.
\]

Here \(z\) is the output and \(t\) is the successor state.
The proof ensemble samples every diagonal entry of \(D_x,D_s\) independently from \(F^*\).
The fixed matrices have shapes \(A\in F^{8\times4}\) and \(C\in F^{4\times8}\).

Use the actual embedded subfield

\[
H_4=(0,1,188,189),\qquad V=16+H_4=(16,17,172,173).
\]

Let \(R=(0,92,176,236)\) and \(\Delta=(0,2,6,4)\).
For phase \(p=2a+b\), where \(0\le a<4\) and \(b\in\{0,1\}\), the ordered labels are

\[
U_p=(\Delta_a+R_{2b}+H_4)\;\Vert\;(\Delta_a+R_{2b+1}+H_4).
\]

These are the eight halves of the old four 16-label phases.
They partition the same 64 labels and avoid \(V\).
The local matrices are

\[
A_{u,j}=(u+v_j)^3+1,\qquad
C_{j,u}=(u+v_j)^3+3.
\]

Thus \(C=A^T+2J\), where \(J\) is the all-one matrix of the required shape.
The construction is asymmetric; an implementation must not substitute the symmetric packet adjoint.
The exact adjoint is

\[
u^*=D_x(z^*+C^Tt^*),\qquad
s^*=D_s(t^*+2A^Tz^*).
\]

The eight-phase schedule keeps the old label geometry.
A four-phase restriction is not needed for the proposed implementation.

## Why the old injection needs repair

The old circuit used \(C_{j,u}=(u+v_j)^3\).
It has four zero square minors in each of phases 0, 1, 4, and 5 of the split schedule.
For example, phase-zero labels at positions 0 and 6, and residue positions 1 and 2, give

\[
\begin{pmatrix}161&242\\243&252\end{pmatrix}
\]

with determinant zero in the AES field.
The old all-subsets update property does not imply full coordinate superregularity of \(C\).
Changing the additive constant to 3 removes these defects for all eight phases.

For each phase, the tests exhaust every square minor of both maps:

\[
2\sum_{j=1}^4\binom8j\binom4j=988.
\]

All 7,904 determinants are nonzero.
This is a finite exhaustive certificate for these explicit constants, not a sampled rank test.

## Feedback and return weights

The emission matrix also has a direct Cauchy proof.
Define

\[
Q(X)=\prod_{v\in V}(X+v)=X^4+X+16^4+16.
\]

For every label \(u\),

\[
A_{u,j}=\frac{Q(u)}{u+v_j}.
\]

The labels and residues are distinct and disjoint, so the Cauchy determinant formula makes every square minor of \(A\) nonzero.

Each \(U_p\) is an affine binary space of size eight.
Summing a polynomial of degree at most six over \(U_p\) gives zero.
One proof uses Lagrange interpolation: the derivative of the additive vanishing polynomial is constant, and the coefficient of degree seven is zero.
Every entry of \(CA\) is the sum of a polynomial of degree at most six.
Consequently,

\[
CA=0.
\]

The extended packet matrix \(\begin{psmallmatrix}I&2A\\C&I\end{psmallmatrix}\) is therefore invertible.
The return matrix \(M:=I+2AC\) satisfies \(M^2=I\).
Direct field algebra gives

\[
(AC)_{u,w}=(u+w)^3+3.
\]

Thus all eight phases have the same ordered return matrix.
Every entry of \(M\) is nonzero, proving singleton return weight eight.

A return event means \(r\ne0\) and \(t=0\).
It forces \(r=Cx\), and the packet output is \(Mx\).
For this width-eight map, the displayed minima hold even without the restriction \(Cx\ne0\):

| Exact input weight | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Minimum weight of \(Mx\) | 8 | 6 | 5 | 4 | 3 | 2 | 2 | 1 |

The tests prove the lower bounds with 11,329 exact shortened-matrix rank checks.
For every input support \(I\) of size \(w\) and output support \(J\) of size one below the claimed minimum,
the matrix \(M_{J^c,I}\) has column rank \(w\).
No nonzero input supported on \(I\) can produce an output supported on \(J\).

A second enumeration establishes attainability at every weight.
For a matrix \(B\) with \(w\) columns, its kernel contains a vector with all coordinates nonzero exactly when

\[
\operatorname{rank}B<w,
\qquad
\operatorname{rank}(B\text{ with column }i\text{ removed})=\operatorname{rank}B
\quad\text{for every }i.
\]

The reverse implication uses \(|F|=256>w\): the union of \(w\) proper coordinate hyperplanes cannot cover the kernel.
The enumeration applies this criterion to every candidate input/output support pair.

## Five-product packed schedule

Each four-label block of either map is a symmetric XOR-circulant.
Its four coefficients \(g\) have XOR sum one, including after either constant shift.
For a four-vector \(x\), write

\[
Zx=(x_0+x_1+x_2+x_3,\ x_1+x_3,\ x_2+x_3,\ x_3).
\]

The transform satisfies \(Z^2=I\).
If \(Zg=(1,g_1,g_2,g_3)\) and \(Zx=(z_0,z_1,z_2,z_3)\), XOR convolution is \(Zw\), where

\[
\begin{aligned}
w_0&=z_0,\\
w_1&=z_1+g_1z_0,\\
w_2&=z_2+g_2z_0,\\
w_3&=z_3+g_1z_2+g_2z_1+g_3z_0.
\end{aligned}
\]

There are five fixed-product sites.
Pack four independent circuits into the four 128-bit lanes of each ZMM:

\[
(C_0x_0,\ C_1x_1,\ A_0(2r),\ A_1(2r)).
\]

Each site becomes one GFNI512 instruction.
The five sites reuse three fixed coefficient ZMMs per phase: the \(g_1\), \(g_2\), and \(g_3\) tables.
Changing the injection constant only changes those tables.
The adjoint swaps the emission and injection block roles and scales the four-record injection result by 2.

The resulting arithmetic count per packet is:

| Component | GFNI512 |
| --- | ---: |
| Four paired circulant maps | 5 |
| Eight input diagonal records | 2 |
| Four state diagonal records | 1 |
| Eta scaling | 1 |
| Total per eight output records | 9 |

The two zeta transforms and their central combination use twelve XOR/XOR3 instructions.
In the forward map, two output additions and one state XOR3 raise the packet total to fifteen.
The adjoint first XORs its two injections, scales that sum by eta, and XORs the result with the incoming state.
Its packet total is sixteen XOR/XOR3 instructions, or thirty-two per sixteen outputs; the forward total is thirty per sixteen outputs.
Two standard 4-by-4 record transposes require at most sixteen 128-bit-lane shuffle instructions.
The initial transpose uses the two input ZMMs and two state copies: eta-scaled in the forward map, unscaled in the adjoint.
The final transpose separates the two injections and two emissions.
Eight phases require 24 fixed coefficient ZMM tables, or 1,536 bytes.

These counts exclude loads, stores, diagonal-matrix expansion, spills, and compiler scheduling.
They are not a native-kernel performance result.
The Python reference checks 49,152 basis/scalar packed-map cases against independent dense multiplication, and 256 full-packet adjoint identities.

## General scaled-Cauchy fallback

The same file constructs four-phase symmetric maps for \((b,s)=(8,4),(16,8),(32,8),(32,16)\).
Let

\[
U_p=pb+\{0,\ldots,b-1\},\quad
V=128+\{0,\ldots,s-1\},\quad
Q(X)=\prod_{v\in V}(X+v),\quad D=Q'.
\]

Here addition is field addition, not integer addition; the displayed integer intervals are binary subspaces.
Set

\[
A_{u,j}=\frac{Q(u)}{D(u+v_j)},\qquad C=A^T.
\]

Every square minor is nonzero by the scaled-Cauchy formula.
The column polynomials have degree \(s-1\); their pairwise products have degree at most \(2s-2\le b-2\).
The same affine-space power-sum argument gives \(CA=0\).

For \(b=2s\), a valid return bound at input weight \(w\le s\) is

\[
\max\{1,\min\{s+1,2(s+1-w)\}\}.
\]

Each half is an involutory circulant, and the cross-half Gram block is scaled Cauchy.
An input in one half gives weight at least \(s+1\); inputs in both halves give at least \(2(s+1-w)\).
For larger weights, invertibility gives the bound one.
For \((b,s)=(32,8)\), singleton return weight is 25, and the generic conditional return bound is \(\max(1,25-w)\).
That latter bound uses \(Cx\ne0\); it is not an unconditional minimum of \(Mx\) when \(w>8\).

| Geometry | Dense GFNI512 per packet | Dense GFNI512 per record |
| --- | ---: | ---: |
| 8/4 | 19 | 19/8 |
| 16/8 | 70 | 35/8 |
| 32/8 | 138 | 69/16 |
| 32/16 | 268 | 67/8 |

These dense counts include input/state diagonals and fold eta into the fixed map constants.
The symmetric 8/4 fallback also has the five-product paired-circulant schedule, reducing its count to nine.
Its return table is only \((5,5,4,2,1,1,1,1)\), so the repaired cubic map is preferable at the same arithmetic count.

Run `python research/workstreams/char2_field_spin/theory/gf8_packet_maps_tests.py` for all fourteen exact tests.
The next gate is a whole-code sparse/middle/dense certificate for the repaired eight-phase map and one regional route.
A native implementation and serial benchmark should follow only if that gate closes.
