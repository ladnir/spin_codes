# Fixed seven-output full-Cauchy packet

This candidate keeps the factored width8/state8 packet but retains only seven of its eight emitted coordinates.
The stronger outer code is GRS[128,84], giving design rate \(84/112=3/4\).
There is one regional routing stage, no second permutation sampled for puncturing, and no multiplication in a larger field.
The local transfer and scalar implementation model pass the independent audit below.
This document does not itself establish a whole-code failure bound or native performance.

## Field, source, and setup

Use \(\mathbb F=\operatorname{GF}(256)\) in the AES polynomial basis, with reduction polynomial `0x11b`.
Integers `0..255` name field elements in that basis; field addition is XOR.
An application record contains sixteen field elements in sixteen byte positions.
The same scalar code and the same setup are applied independently to those positions.

For requested source length \(K\), define

\[
L:=64\left\lceil\frac{K}{64\cdot84}\right\rceil,\qquad
K_{\rm pad}:=84L,\qquad N_{\rm in}:=128L,\qquad
T:=N_{\rm in}/8,\qquad N_{\rm out}:=7T=112L.
\]

Pad the source with zero records to length \(K_{\rm pad}\).
For each row \(r\in\{0,\ldots,L-1\}\), interpret its 84 source records as coefficients of a degree-below-84 polynomial.
Evaluate at the field nodes `0..127` to obtain one GRS[128,84] word.
Every nonzero scalar outer row has at least 45 nonzero coordinates.
Fixed nonzero coordinate multipliers would preserve the proof, but this specification uses multiplier one.

The coefficient basis may be the normalized additive basis used by the existing full128 transform.
To fix its meaning, put

\[
S_i(X):=\prod_{u=0}^{2^i-1}(X+u),\quad
\beta_i(X):=S_i(X)/S_i(2^i),\quad
b_j(X):=\prod_{i:\,j_i=1}\beta_i(X).
\]

Here \(j_i\) are the bits of \(j\), and the product for \(j=0\) is one.
Each \(b_j\) has degree \(j\), so the first 84 span exactly the required polynomial space.
In four-row source tiles, coefficient \(j\) of row \(r\) has source index
\(336\lfloor r/4\rfloor+4j+(r\bmod4)\).

The setup samples the following independent objects:

- one permutation \(\pi_c\gets S_L\) for each outer coordinate \(c=0,\ldots,127\);
- one input scalar \(d_{\ell,i}\gets\mathbb F^*\) for each packet \(\ell<T\) and input coordinate \(i<8\);
- one state scalar \(e_{\ell,i}\gets\mathbb F^*\) for each packet \(\ell<T\) and state coordinate \(i<8\).

The routed stream is \(v_{cL+t}:=\operatorname{outer}_{\pi_c(t),c}\).
Process consecutive blocks of eight records as packets, in increasing stream order.
The state is carried across region boundaries; it is not reset once per region.
All sixteen byte positions of a record use the same permutation and scalar coefficients.
No factor of sixteen is needed in the setup-failure union bound: any nonzero record-valued source has a nonzero scalar lane, whose support is contained in the record support.

## Explicit packet maps

Form the 16-by-16 Cauchy matrix

\[
Q_{i,j}:=(i+(128+j))^{-1},\qquad 0\le i,j<16,
\qquad Q=\begin{pmatrix}B&A_0\\C_0&E\end{pmatrix},
\]

with eight-by-eight blocks.
Let \(C_1:=E^{-1}C_0\).
The row sums of \(B,A_0,C_1\) are 15, 161, and 46.
Define scalar field elements \(l:=161^{-1}=124\) and \(r:=46^{-1}=162\), and set

\[
F:=lrB,\qquad A:=lA_0/2,\qquad C:=rC_1,\qquad \eta:=2.
\]

These are exactly `gf8_mds_packet.fullstate_matrices(normalized=True)`.
Write \(A_e:=2A\).
All three matrices are symmetric XOR-circulants; their exact factorization is

\[
F=242A_eC.
\]

The kernels of \(A_e\) and \(C\), respectively, are

```text
172  59 114 102 161 248  77 150
 52 172 158  17 209  32 250  29
```

Both sum to one, permitting the checked nineteen-product circulant circuit.

Initialize the eight-record persistent state to \(s_0=0\).
For packet \(\ell\), define \(D_{x,\ell}:=\operatorname{diag}(d_{\ell,0},\ldots,d_{\ell,7})\) and similarly \(D_{s,\ell}\) from the state scalars.
Compute

\[
x:=D_{x,\ell}v_\ell,\qquad h:=D_{s,\ell}s_\ell,\qquad m:=Cx,
\]
\[
y_\ell:=A_e(h+242m)=Fx+A_eh,\qquad
s_{\ell+1}:=h+m.
\]

Emit only \((y_{\ell,0},\ldots,y_{\ell,6})\).
Coordinate seven is always punctured; the puncture is fixed, not random.
Discard the final state \(s_T\).

For an explicit packed serialization, group four chronological packets and store

\[
\operatorname{out}_{28g+4j+a}:=y_{4g+a,j},
\qquad 0\le j<7,\quad0\le a<4.
\]

This is only a fixed reordering of the retained coordinates.
It preserves weight and makes each group occupy seven aligned ZMM-sized record groups.
Packets within each group are stored in forward order; the adjoint processes those four packet slots in reverse order.

## Local proof audit

The normalized full packet matrix is obtained from \(Q\) by invertible operations on all eight bottom rows and nonzero diagonal row/column scalings.
Its bottom-right block is identity, so it is not itself fully superregular.
The properties used below concern selected top rows and, for returns, all bottom rows.

Fix input support \(I\), incoming-state support \(H\), and selected retained-output zero rows \(Z\subseteq\{0,\ldots,6\}\).
Let their sizes be \(w,h,z\).
The output-only coefficient matrix has rank \(\min(z,w+h)\).
Appending every state-return equation gives rank \(\min(8+z,w+h)\).
The original Cauchy matrix proves these ranks for every support choice; no symmetry or randomness of the puncture is required.

Eliminating \(\min(w,z)\) input pivots leaves a superregular state-only matrix of rank \(\max(0,z-w)\).
This is the state-root profile used in the continuation entry of `gf8_punctured_screen.py`.
The root helper is applied to fresh state diagonals after fixing the incoming source forms.
The input-pivot charge uses the independent current input diagonals.

The seven-by-eight birth and return maps have exact minimum output weight

\[
b(w)=\max(0,8-w),\qquad 1\le w\le8.
\]

In particular, \(w=8\) can emit zero.
The screen retains this case in the birth term and in the live-to-zero return term with \(z=7\).
For a canonical incoming zero-coordinate flat of rank \(k\), the joint return condition requires \(k\le w-z-1\).
Its return sum uses exactly this cutoff and one joint charge \(\tau^{-(8+z)}\).
It does not reuse those diagonals for another independent probability charge.

Two other zero-output cases matter:

- With zero input and live state, all seven retained outputs can vanish while the state remains live. The continuation term for seven roots is present.
- With zero input and zero state, output and next state remain zero. Its transfer entry is exactly one.

A nonzero zero-to-zero input packet is impossible here because \(C\) is invertible.
The generic screen includes such a branch only when \(w>s+z\), which never occurs at width/state eight.
The regional product starts at zero state and sums both terminal states, matching the discarded terminal state in the construction.

The input support coefficients remain \(\binom8w\).
The positive output expansion uses \(\binom7z y^{7-z}(1-y)^z\), not \(\binom8z\).
The outer source bound changes to dimension 84, while the packet count and inverse-input weight use \(N_{\rm in}\).
The output threshold and its Chernoff charge use \(N_{\rm out}\).
These distinctions are implemented consistently in the reviewed screen.

The independent audit found no mathematical defect in these local substitutions.
The supplied three tests pass, including exact agreement with the old transfer when no row is punctured.
Floating search scores remain proposals until directed arithmetic and complete class coverage are checked.

## Global row rank and the dense preimage argument

Fix the state diagonals and regard the inner as a map from the post-input-diagonal stream to retained outputs.
Call its matrix \(H\in\mathbb F^{N_{\rm out}\times N_{\rm in}}\).
Every diagonal packet block is the first seven rows of \(F\), with rank seven.
Choosing the first seven input columns in every packet produces a block lower-triangular square matrix with invertible diagonal blocks.
Thus \(H\) has full row rank for every realization of the state diagonals.

The complete inner before routing is \(HD_x\), where \(D_x\) is the diagonal containing all input scalars.
The matrix \(H\) is independent of the input scalars and the regional permutations.
For each fixed output \(y\), its preimage set under \(H\) contains exactly

\[
256^{N_{\rm in}-N_{\rm out}}
\]

vectors.
Conditioning on \(H\), the existing regional membership theorem applies separately to every fixed preimage.
A union bound therefore adds \(8(N_{\rm in}-N_{\rm out})\) bits to the failure exponent and uses an output ball of length \(N_{\rm out}\).
For a dense-class bound, count only preimages in the specified inverse-input weight range; bounding their number by the full preimage count is safe.
Their dependence on the conditioned state diagonals does not consume the still-independent input diagonals or routing randomness.

Full row rank of \(H\) does not imply that its restriction to the outer code is injective.
The distance proof must include nonzero sources mapping to output zero.
A successful distance certificate then gives injectivity, dimension \(K_{\rm pad}\), and design rate three quarters.

The new preimage factor requires a new dense replay and a new join.
The former 0.28 dense cutoff cannot be imported.
An independent floating replay of the same endpoint formula begins passing near raw-input density 0.44890 for this seven-output candidate.
That number is a diagnostic, not a certified endpoint.

## Geometry and implementation accounting

At requested \(K=2^{20}\), the concrete geometry is

| Quantity | Value |
| --- | ---: |
| Outer rows \(L\) | 12,544 |
| Padded source dimension | 1,053,696 |
| Raw inner input records | 1,605,632 |
| Retained output records | 1,404,928 |
| Forbidden output weights | 0 through 140,492 |
| Requested finite rate | approximately 0.746356 |
| Padded design rate | exactly 3/4 |

For comparison, changing only the unpunctured outer to [128,84] gives rate 84/128, not three quarters.
Puncturing reallocates redundancy between the outer and inner while preserving the target overall design rate.
It also raises raw inner and routing work by about one seventh at fixed source length.
The exact K20 raw-work increase relative to the [128,96] control is about 14.62%, including padding.

Four packets consume 32 raw records and emit 28 records.
The factored arithmetic remains 62 GFNI512 and 90 XOR/XOR3 instructions per four packets before any puncture-specific simplification.
This is about 35.43 GFNI512 per sixteen retained records, not 31.
The eight state coordinates and their updates cannot be punctured.

The fixed packed serialization can remove one pair of 4-by-4 record transposes:

1. Forward: load and pack the 32 routed inputs, then evaluate their four injection maps.
2. Unpack the moments and advance the four state steps in chronological order.
3. Pack the four modified histories and evaluate the four emission maps.
4. Store the seven retained coordinate registers; omit coordinate seven without unpacking outputs.

This uses six 4-by-4 transposes, or 48 record-shuffle instructions, per four-packet batch in the proposed layout.
In the adjoint, load seven packed registers and set the missing eighth register to zero.
Traverse packets in reverse chronology, using the same factored equations with \(p=(p_0,\ldots,p_6,0)\):

\[
m:=A_e^Tp,\qquad
u:=D_xC^T(h+242m),\qquad
h_{\rm previous}:=D_s(h+m).
\]

Initialize the terminal dual state to zero because the forward terminal state is discarded.
Discard the final initial-state dual because the forward initial state is fixed to zero.
The inner adjoint produces all 32 raw records for the existing inverse route and outer adjoint.
No dynamic allocation or additional global route is required by this layout.
Caller-owned scratch still needs space for the larger raw stream; exact aliasing is an implementation obligation, not promised by this scalar model.

The cost counts exclude diagonal expansion, routing, loads/stores, spills, and the outer transform.
The existing full128 coefficient transform can represent the degree-below-84 outer by zero-extending its source and selecting the first 84 adjoint coefficients.
Its circuit work per row need not fall merely because fewer source coefficients are retained.
The larger row count must be included in any whole-call comparison.
This algebraic audit did not produce a native encoder or benchmark. The
subsequent [native checkpoint](GF8_PUNCTURED_NATIVE.md) implements this
layout, validates aliasing, and reports complete timed calls.

## Independent scalar receipt

Run

```text
python research/workstreams/char2_field_spin/theory/gf8_punctured_packet_audit.py
```

All 64 global-map cases passed: full row rank, an explicit invertible pivot minor, factorization independent of input diagonals, and literal-matrix adjoint agreement.
They cover four output widths, one/two/four/five packets, and independent nonzero diagonal samples.
All 256 four-packet cases passed the fixed packing, puncture zero-fill, and reverse-state chronology checks, including nonzero terminal dual states.
The universal rank claim follows from the packet pivot argument above, not from sampling these examples.

The subsequent [whole-code certificate](GF8_PUNCTURED_WHOLE_CERTIFICATE.md)
closes the joined sparse/middle/dense gate with an independent coverage
audit. The native checkpoint above then passes the implementation gate.
