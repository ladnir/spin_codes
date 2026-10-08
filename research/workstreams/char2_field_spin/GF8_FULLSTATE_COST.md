# Full-state8 packet: exact normalization and proposed cost

The width8/state8 candidate replaces the cubic packet's four state records with eight.
It has an MDS return map, while retaining two structured Cauchy maps.
This note gives an exact scalar implementation and a proposed four-packet adjoint schedule.
It does not certify whole-code distance or report native performance.

## Packet and local guarantees

All bytes denote elements of GF(256) in the AES polynomial basis, with reduction polynomial `0x11b`.
Addition is XOR.
For `a` outside the binary subspace `U={0,...,7}`, define the symmetric XOR-circulant

\[
(G_a)_{i,j}:=(a+i+j)^{-1},\qquad 0\le i,j<8.
\]

Set \(A:=G_{128}\), \(C:=G_{192}\), and \(\eta:=2\).
Fix nonzero input and state diagonals \(D_x,D_s\).
For eight input records \(v\) and eight persistent state records \(s\), the forward packet computes

\[
x:=D_xv,\quad r:=D_ss,\qquad y:=x+2Ar,\quad t:=r+Cx.
\]

It emits \(y\) and retains \(t\).
Each record contains sixteen independent field coordinates; the scalar equations apply to every coordinate.

Both \(A\) and \(C\) are superregular by the Cauchy determinant formula.
For nonzero input of weight \(w\), a zero incoming state gives outgoing-state weight at least \(9-w\).
The emitted birth weight remains exactly \(w\); this construction does not add direct input diffusion.
On a return to zero state, \(r=Cx\), so the emitted vector is

\[
y=Mx,\qquad M:=I+2AC.
\]

The exact reference checks all 12,869 square minors of \(M\); every minor is nonzero.
Consequently its exact minimum output weights, for input weights one through eight, are

\[
(8,7,6,5,4,3,2,1).
\]

The full 16-by-16 packet matrix is invertible.
Unlike the cubic packet, this construction has \(CA\ne0\).
Any sparse or middle proof must therefore justify its transfer bound for this packet, not import the old zero-product hypothesis.

## Exact scalar normalization

Let \(S_a:=\sum_{u\in U}(a+u)^{-1}\).
Direct field arithmetic gives \(S_{128}=15\) and \(S_{192}=74\).
Define

\[
\widehat A:=A/15,\quad \widehat C:=C/74,\quad
\widehat\eta:=2\cdot15\cdot74=13,\quad \widehat s:=s/74.
\]

The equivalent packet is

\[
x:=D_xv,\quad \widehat r:=D_s\widehat s,\qquad
y:=x+13\widehat A\widehat r,\quad
\widehat t:=\widehat r+\widehat Cx.
\]

Its emitted vector is identical and \(t=74\widehat t\).
Scalar state rescaling commutes with every diagonal realization.
Thus a whole traversal from zero state emits exactly the same vector, not merely an equal-distance code.

For the adjoint, let \(p\) be the dual emitted vector and \(h\) the dual next state in normalized coordinates.
The reverse step is

\[
u:=D_x(p+\widehat C^Th),\qquad
g:=D_s(h+13\widehat A^Tp).
\]

Here both normalized matrices are symmetric.
The corresponding original terminal dual state is \(h/74\), and the normalized previous dual state is 74 times the original one.
This dual rescaling matters when checking a packet with nonzero boundary state.

The normalized convolution kernels are

```text
Ahat: 40 172 107 122 26 56 235 93
Chat: 201 147 239 156 153 97 131 83
```

Each kernel sums to one.
Partial fractions also give the independently checked identity

\[
G_aG_b=(S_a+S_b)G_{a+b},\qquad a+U\ne b+U.
\]

In particular, \(AC=(15+74)G_{64}\).

## Nineteen-product convolution

For an eight-element vector \(v\), define its superset zeta transform by

\[
(Zv)_m:=\sum_{j\supseteq m}v_j,
\]

where the indices are three-bit masks.
The transform uses twelve XORs and satisfies \(Z^2=I\).
For a normalized kernel \(k\), set \(g:=Zk\); then \(g_0=1\).
After setting \(z:=Zv\), compute the following assignments in descending index order:

```text
z7 ^= g1*z6 + g2*z5 + g3*z4 + g4*z3 + g5*z2 + g6*z1 + g7*z0
z6 ^= g2*z4 + g4*z2 + g6*z0
z5 ^= g1*z4 + g4*z1 + g5*z0
z4 ^= g4*z0
z3 ^= g1*z2 + g2*z1 + g3*z0
z2 ^= g2*z0
z1 ^= g1*z0
```

Applying \(Z\) again gives the XOR-circulant product.
Descending order preserves every input needed by later assignments.
There are exactly nineteen multiplication sites and thirteen XOR/XOR3 sites between the transforms.
The complete circuit uses nineteen fixed multiplications and thirty-seven XOR/XOR3 instructions.
The transformed kernels are

```text
Z(Ahat): 1 179 167 39 148 101 182 93
Z(Chat): 1 61 163 207 40 50 208 83
```

All fourteen nonidentity constants are nonzero and distinct from one.
The counts do not rely on accidental zero coefficients.

## Four-packet reverse schedule

List four consecutive packets in reverse traversal order, newest first, as \(p_0,p_1,p_2,p_3\).
Let \(h_0\) be the terminal dual state.
The schedule computes

\[
m_j:=\widehat A^Tp_j,\qquad
h_{j+1}:=D_{s,j}(h_j+13m_j),\qquad
u_j:=D_{x,j}(p_j+\widehat C^Th_j).
\]

Only the middle recurrence is sequential.
All four input moments can precede it; all four state emissions can follow it.
In particular, packet \(j\) emits from \(h_j\), not \(h_{j+1}\).

Each packet occupies two ZMM registers, with four 128-bit records per register.
Two 4-by-4 record transposes place coordinate \(i\) of all four packets in one ZMM.
The eight resulting ZMMs carry four independent copies of the nineteen-product circuit.
Each transpose uses eight `VSHUFI32X4` instructions.
Packing and unpacking one batch of maps therefore uses 32 record shuffles.

The implementation plan is explicit:

1. Retain the eight ZMMs containing the original packet inputs.
2. Pack a working copy, evaluate all four input moments, and unpack it.
3. Advance the four state steps in order, replacing consumed moment slots with their incoming states.
4. Pack those four saved states, evaluate their four emissions, and unpack them.
5. Add the retained inputs, apply their diagonals, and store the outputs.

No extra route, heap allocation, or caller scratch is required by this register schedule.
The four-packet batch size divides the existing regional sizes.
A future implementation must still handle its API's smaller tails explicitly.

| Operation | Per four packets / 32 records | Per 16 records |
| --- | ---: | ---: |
| Two fixed map batches, GFNI512 | 38 | 19 |
| Input diagonal, GFNI512 | 8 | 4 |
| State diagonal, GFNI512 | 8 | 4 |
| Eta scaling, GFNI512 | 8 | 4 |
| Total GFNI512 | 62 | 31 |
| XOR/XOR3 | 90 | 45 |
| 128-bit-record shuffles | 64 | 32 |

For comparison, the current paired cubic8 adjoint uses 18 GFNI512, 32 XOR/XOR3, and 32 record shuffles per sixteen outputs.
The full-state candidate trades thirteen extra GFNI512 instructions and thirteen extra Boolean instructions for stronger state births and MDS returns.
Its state critical path is one addition followed by a diagonal multiplication per packet.
Input-moment multiplication does not depend on the preceding state.

A destructive hand schedule can keep at most eight original-input registers, eight working registers, and two carry-state registers live.
Seven map constants and two product temporaries raise this to 27 ZMMs; one diagonal-expansion temporary gives 28.
The constants for the two maps need not coexist.
The 4-by-4 transposes can reuse their sources with two temporaries.
This is a proposed register allocation, not compiler-generated evidence of spill-free code.

The arithmetic table excludes payload loads/stores, routing, diagonal expansion, spills, and the outer transform.
It is not a prediction of wall-clock speed.

## Exact receipt and next gate

Run

```text
python research/workstreams/char2_field_spin/theory/gf8_fullstate_cost.py
```

The checked scalar reference covers:

- all 12,869 square return minors and rank16 of the extended packet;
- 4,096 basis/scalar circuit cases against independent dense8 multiplication;
- 8,256 forward/adjoint basis/scalar or random cases for the state normalization;
- the same 8,256 cases against an independent dense16 matrix and its literal transpose;
- 128 four-packet batches against sequential reverse traversal, including nonzero terminal state and independent nonzero diagonals.

These checks passed locally. No native kernel was built or benchmarked.
The next gate is a complete sparse/middle/dense certificate at the requested higher outer rate.
Only a positive proof result should trigger native implementation and serial measurements against paired cubic8 and the existing width16 controls.
