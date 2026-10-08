# Full-Cauchy width8/state4 packet

The cubic packet leaves a zero-state birth's emitted weight unchanged.
The full-state8 alternative strengthens returns, but also leaves that birth weight unchanged.
This candidate instead makes both the emitted birth map and the return map MDS, using four state records.
Its local algebra and packed scalar reference are checked below.
Whole-code distance remains a separate proof obligation.

## Concrete packet

Work over GF(256) in the AES polynomial basis, with reduction polynomial `0x11b`.
Byte addition means XOR.
Define the 12-by-12 Cauchy matrix

\[
T_{i,j}:=\frac1{i+(128+j)},\qquad 0\le i,j<12,
\qquad
T=\begin{pmatrix}B&A_0\\C_0&E\end{pmatrix},
\]

where the first block row and column have size eight.
The twelve row nodes and twelve column nodes are disjoint.
The Cauchy determinant formula therefore makes every square minor of \(T\) nonzero.

Set \(C:=E^{-1}C_0\) and \(A:=A_0/2\).
For an eight-record input \(v\), four-record state \(s\), and nonzero diagonal matrices \(D_x,D_s\), compute

\[
x:=D_xv,\quad r:=D_ss,\qquad
y:=Bx+2Ar,\quad t:=r+Cx.
\]

The packet emits \(y\) and retains \(t\).
Its extended matrix before the diagonals is

\[
P=\begin{pmatrix}B&A_0\\C&I\end{pmatrix}
 =\operatorname{diag}(I_8,E^{-1})T.
\]

This packet uses one regional route, as before.
It adds a local input mixer, not another route.
No identity \(CA=0\) is assumed or needed for the local claims here.

The four state coordinates form the binary additive space \(\{0,1,2,3\}\).
In this geometry, \(E_{i,j}=(128+i+j)^{-1}\), its row sum is the field element 2, and \(E^{-1}=E/4\).
The exact reference verifies that inverse directly.

## Birth, return, and state-map properties

A zero incoming state gives \(y=Bx\).
Since \(B\) is superregular, input weight \(w\) implies emitted weight at least \(9-w\).
This is the exact minimum for each \(1\le w\le8\).

On return to zero state, \(r=Cx\), so

\[
y=Sx,\qquad S:=B+A_0E^{-1}C_0.
\]

For equally sized row and column sets \(I,J\subseteq\{0,\ldots,7\}\), the Schur determinant identity gives

\[
\det S_{I,J}
=\frac{\det T_{I\cup\{8,\ldots,11\},\,J\cup\{8,\ldots,11\}}}{\det E}.
\]

Both determinants on the right are nonzero.
Thus every square minor of \(S\) is nonzero, and its exact return minima are also

\[
(8,7,6,5,4,3,2,1).
\]

The emission matrix \(A\) is a scaled Cauchy submatrix.
The injection matrix \(C\) is also superregular, although the row operation defining it does not preserve arbitrary minors in general.
To prove this specific claim, observe that \([C,I]=E^{-1}[C_0,E]\) generates the same MDS code as the last four rows of \(T\).
Complete any square minor of \(C\) with complementary identity columns; its determinant is a nonzero four-column minor of \([C,I]\).
In particular, an input of weight \(w\le4\) creates state weight at least \(5-w\).

The exact script exhausts all 12,869 square minors of each of \(B,S\), and all 494 square minors of each of \(A,C\).
All are nonzero.

## Joint ranks available to a probability proof

Fix selected zero-output rows \(Z\), active input columns \(I\), and active incoming-state columns \(H\).
Write their cardinalities as \(z,w,h\).
The selected output equations have coefficient matrix

\[
[B_{Z,I},(A_0)_{Z,H}].
\]

This is a Cauchy submatrix, so its rank is \(\min(z,w+h)\).
On return, append all four outgoing-state equations.
The stacked matrix then has rank

\[
\operatorname{rank}
\begin{bmatrix}
B_{Z,I}&(A_0)_{Z,H}\\ C_{:,I}&I_{:,H}
\end{bmatrix}
=\min(z+4,w+h).
\]

Indeed, multiplying only the four bottom rows by \(E\) recovers the corresponding submatrix of \(T\).
This argument permits arbitrary \(I,H,Z\); it does not require contiguous supports.
Nonzero diagonal column scalings preserve both ranks.
Any nonzero solution of the joint equations therefore requires \(w+h>z+4\).

For \(z>w\), eliminate all \(w\) input columns using any \(w\) selected output rows.
The residual \((z-w)\)-by-4 state matrix is superregular by the same Schur-minor identity.
Its rank is \(\min(z-w,4)\).
Thus four residual state equations force the incoming state to zero.
These statements supply local rank facts; they do not by themselves justify a state-history moment operator.

## Normalized packed circuit

The implementation-oriented reference scales the first and second four-record output cosets by
\(L=(100,8)\), and the corresponding input cosets by \(R=(200,16)\).
Write \(L,R\) also for the resulting diagonal eight-by-eight matrices.
The returned API uses

\[
F:=LBR,\quad \widehat A:=LA_0/2,\quad \widehat C:=CR,\quad \eta=2.
\]

The extended matrix is \(\operatorname{diag}(L,I)P\operatorname{diag}(R,I)\).
Every preceding support and rank claim is preserved.
Input scaling can be coupled through the nonzero input diagonals, while output scaling preserves Hamming weight.
This is a support-equivalent normalization, not a claim that unadjusted seeded encoders emit identical bytes.

Every four-by-four block is an XOR-circulant on \(\{0,1,2,3\}\).
The two kernels of the effective emission map \(2\widehat A\) and two kernels of \(\widehat C\) are

```text
effective A, upper: 107 132  37 203
effective A, lower: 249  16 191  87
C, left:           245 232 196 216
C, right:          175 180 150 140
```

Each kernel sums to one, so four independent circulants fit the existing five-product circuit's four 128-bit lanes.
Eta is already folded into the effective emission constants; it adds no separate multiplication.

The four kernels of \(F\), in row-major block order, are

```text
155  26 175  95
 72  34 245  61
 72  34 245  61
  1 127 255 154
```

Their sums are 113, 162, 162, and 27, so the unit-sum shortcut does not apply.
For a general four-point kernel, let \(g:=Zk\) and \(z:=Zv\), where \(Z\) is the four-point superset zeta transform.
Compute

\[
u_0=g_0z_0,\quad
u_1=g_0z_1+g_1z_0,\quad
u_2=g_0z_2+g_2z_0,\quad
u_3=g_0z_3+g_1z_2+g_2z_1+g_3z_0,
\]

and return \(Zu\).
Four parallel copies require nine GFNI512 instructions, twelve XOR/XOR3 instructions, and two four-record transposes.
The two equal off-diagonal kernels share constants, but this count does not assume fewer multiplication sites.

The complete four-state packet uses the following fixed arithmetic in either direction:

| Operation | Per eight output records | Per sixteen output records |
| --- | ---: | ---: |
| Paired state maps, GFNI512 | 5 | 10 |
| Input mixer, GFNI512 | 9 | 18 |
| Input/state diagonals, GFNI512 | 3 | 6 |
| Total GFNI512 | 17 | 34 |
| XOR/XOR3 | 27 | 54 |
| Record shuffles | 32 | 64 |

The scalar adjoint reference computes the input mixer from the emitted dual input alone.
Its state recurrence can therefore be scheduled independently of that mixer.
The table counts two explicit packed circuits, without assuming shared transposes or a future batching optimization.
It excludes payload loads/stores, routing, diagonal expansion, spills, and the outer transform.
No native instruction schedule or speed claim is established.

For context, paired cubic8 uses 18 GFNI512 per sixteen records, and the proposed full-state8 schedule uses 31.
The older dense width16/state8 accounting is 70.
The new packet costs more arithmetic than cubic8 but strengthens both emitted births and returns while retaining four state records.

## Receipt and next gate

`theory/gf8_mds_packet.py` exposes `matrices(normalized=True)`, returning \((F,\widehat A,\widehat C,2)\).
It also exposes the extended and return matrices, selected-support rank, and the packed scalar forward/adjoint reference.
Run

```text
python research/workstreams/char2_field_spin/theory/gf8_mds_packet_tests.py
```

The width8/state4 tests pass.
Besides the exhaustive minor checks, the receipt includes 1,834 selected-rank cases covering every cardinality triple.
It verifies the exact Cauchy factorization used for the universal rank argument.
All 6,208 packed basis/scalar or random forward/adjoint cases agree with an independently assembled dense12 matrix and its literal transpose.

The next step is the new sparse/middle/dense proof screen at the requested higher rate.
The old identity-feedthrough transfer cannot be reused unchanged because \(F\ne I\).
A native implementation should wait for that probability bound; none was built or benchmarked here.

## Full-Cauchy width8/state8 extension

Use the same construction with sixteen row nodes `0..15` and sixteen column nodes `128..143`.
Partition both axes after eight coordinates.
Then \(B=E=G_{128}\) and \(A_0=C_0=G_{136}\), where \(G\) now denotes an eight-point XOR-circulant.
The state row normalization remains \(C=E^{-1}C_0\).
Here the row sums of \(B,A_0,C\) are 15, 161, and 46.
Since \(E^2=15^2I\), its inverse is \(E/15^2\).

The emission and injection maps are now square MDS matrices.
Both emitted birth and return minima remain \(9-w\).
The outgoing state on a zero-state birth also has weight at least \(9-w\), rather than the four-state bound \(\max(0,5-w)\).
The selected-output rank remains \(\min(z,w+h)\).
Appending all eight outgoing-state rows gives rank \(\min(z+8,w+h)\).
The output-only residual state matrix has rank \(\min(\max(0,z-w),8)\).
These follow from the same Cauchy and Schur arguments, with eight bottom rows.

For the packed reference, use scalar output and input scales \(L=124\) and \(R=162\).
The effective emission and injection kernels then sum to one; the feedthrough kernel sums to 242.
Thus the two state maps each use nineteen-product eight-point circuits.
The input mixer uses a nineteen-product unit-sum circuit followed by eight multiplications by 242.
All maps can batch four independent packets in four 128-bit lanes.

The following three-map schedule is a valid baseline; the next section improves it algebraically.
For reverse traversal, first compute all four emission-input moments \((2\widehat A)^Tp_j\).
Next overwrite the retained raw inputs with \(F^Tp_j\).
Advance the four state updates, replacing each consumed input moment with its incoming state.
Finally compute \(\widehat C^Th_j\) in that working bank and add the retained feedthrough results.
This order avoids keeping three eight-ZMM banks simultaneously.
It retains the pointwise adjoint chronology: output \(j\) uses \(h_j\), and the update produces \(h_{j+1}\).

| Operation | Per four packets / 32 records | Per 16 records |
| --- | ---: | ---: |
| Input mixer, GFNI512 | 27 | 13.5 |
| Both state maps, GFNI512 | 38 | 19 |
| Input/state diagonals, GFNI512 | 16 | 8 |
| Total GFNI512 | 81 | 40.5 |
| XOR/XOR3 | 127 | 63.5 |
| Record shuffles | 96 | 48 |

The fractional per-sixteen figures are amortized over the complete 32-record batch.
A hand schedule has sixteen data registers and two carry-state registers, plus at most eight constants and two temporaries.
Allowing one additional diagonal-expansion temporary gives a proposed peak of 29 ZMMs.
Compiler allocation and speed are untested.

Compared with the four-state full-Cauchy schedule, the additional arithmetic is 6.5 GFNI512 instructions per sixteen records.
The batch uses fewer record shuffles, while doubling persistent state and strengthening state births and return constraints.
This is a credible proof-headroom tradeoff, not evidence that the eight-state encoder is faster.

The scalar API is `fullstate_matrices(normalized=True)` with the same \((F,A,C,2)\) interface.
`fullstate_packet` and `fullstate_four_packet_adjoint` model the circuits and chronology.
All 12,869 minors of each of the four square maps are nonzero.
The additional test checks 8,256 basis/scalar or random forward/adjoint cases against a dense16 matrix and its literal transpose.
It also checks 128 four-packet batches against sequential traversal.
The following factored schedule is also checked by those tests.

## Removing the eight-state feedthrough circuit

Write \(A_e:=2\widehat A\) for the effective eight-state emission map.
The concrete normalized matrices satisfy the exact identity

\[
F=242A_e\widehat C.
\]

This identity follows from the eight-point Cauchy convolution product: the kernels of \(A_e\) and \(\widehat C\) are scaled \(G_{136}\) and \(G_8\).
Their product is scaled \(G_{128}\), as is \(F\).
Its row sum is one, while \(F\) has row sum 242.
The scalar reference also checks every matrix entry of this identity.

Consequently the exact forward packet can compute

\[
m:=\widehat Cx,\qquad
y:=A_e(r+242m),\qquad t:=r+m.
\]

The adjoint can compute

\[
m:=A_e^Tp,\qquad
u:=D_x\widehat C^T(h+242m),\qquad
g:=D_s(h+m).
\]

These are the same packet maps, not a new ensemble or a changed basis.
Four reverse packets require only two batches of nineteen-product circulants.
After the input moment is consumed, its register slots can hold \(h_j+242m_j\).
The state update still uses \(h_j+m_j\), so those two combinations must not be interchanged.

| Factored width8/state8 operation | Per 32 records | Per 16 records |
| --- | ---: | ---: |
| Both circulant batches, GFNI512 | 38 | 19 |
| History scaling, GFNI512 | 8 | 4 |
| Input/state diagonals, GFNI512 | 16 | 8 |
| Total GFNI512 | 62 | 31 |
| XOR/XOR3 | 90 | 45 |
| Record shuffles | 64 | 32 |

There is no need to retain the original input after computing \(m\).
A proposed allocation uses eight working registers, two carry-state registers, seven constants, two temporaries, and one diagonal-expansion temporary: twenty ZMMs.
The exact Python checks compare both factored directions against the original dense matrix and the unfactored reference.
They also compare the two four-packet schedules.
This supersedes the 40.5-GFNI arithmetic baseline above, without establishing native speed.

## Full-Cauchy width8/state16 extension

Take the 24-by-24 Cauchy matrix with row nodes `0..23` and column nodes `128..151`, again splitting after eight.
The sixteen state labels `8..23` form the affine additive space
\(8+\operatorname{span}_{\mathbb F_2}(1,2,4,24)\).
Thus \(E\) is an XOR-circulant on sixteen indices, despite the nonstandard highest basis vector.
Its row sum is 16 and \(E^{-1}=E/16^2\).
The exact script checks both inverse products and the complete state-row factorization.

The injection matrix \(C=E^{-1}C_0\) is 16-by-8 and superregular by the same systematic-MDS argument.
A zero-state birth has outgoing-state weight at least \(17-w\).
In particular, no nonzero packet of width eight can start and finish at zero state.
An arbitrary return requires at least \(17-w\) live incoming-state coordinates.
Emitted birth and return minima remain \(9-w\).
The joint selected-output/return rank is now \(\min(z+16,w+h)\).

The two eight-coordinate state blocks have original emission row sums 161 and 79, and injection row sums 184 and 24.
Rescale the state blocks by \(\beta=(184,24)\), and rescale outputs by 234.
The input coordinates are unchanged.
The normalized injection blocks \(C_0',C_1'\) both sum to one.
The effective emission blocks \(A_0',A_1'\) sum to one and 110, respectively.
The normalized feedthrough map sums to 177 and satisfies

\[
F=177A_0'C_0'.
\]

All four blocks are eight-point XOR-circulants.
Each state rescaling is scalar within its block, so it commutes with every state diagonal.
As before, the output scaling preserves support, but an unadjusted seeded byte-for-byte encoder comparison is not implied.

In reverse traversal, compute

\[
m_0:=(A_0')^Tp,\quad m_1:=(A_1')^Tp,
\]
\[
u:=D_x\bigl((C_0')^T(h_0+177m_0)+(C_1')^Th_1\bigr),\qquad
g_b:=D_{s,b}(h_b+m_b)\quad(b=0,1).
\]

The second emission block uses a unit-sum nineteen-product circuit followed by scaling by 110.
The first state history is modified by 177; the second is unmodified.
For a four-packet batch, compute the first moments in a working bank and overwrite the raw-input bank with the second moments.
Both state recurrences can then replace their consumed moments with the required histories.
No third input or feedthrough bank is needed.

| Factored width8/state16 operation | Per 32 records | Per 16 records |
| --- | ---: | ---: |
| Four unit-sum circulant batches, GFNI512 | 76 | 38 |
| Nonunit emission scaling, GFNI512 | 8 | 4 |
| Modified-history scaling, GFNI512 | 8 | 4 |
| Input/state diagonals, GFNI512 | 24 | 12 |
| Total GFNI512 | 116 | 58 |
| XOR/XOR3 | 180 | 90 |
| Record shuffles | 128 | 64 |

A proposed allocation uses sixteen working registers, four carry-state registers, seven constants, two temporaries, and one diagonal-expansion temporary: thirty ZMMs.
This is close enough to the architectural limit that native register allocation must be inspected before making performance claims.
The arithmetic remains below the older dense width16/state8 count of seventy GFNI512 per sixteen records.
The additional shuffles, longer state, routing, memory traffic, and compiler behavior remain unmeasured.

`largestate_matrices()` returns the normalized \((F,A,C,2)\) interface, with the effective emission blocks stored as \(2A\).
The exact receipt checks 12,869 birth minors, 12,869 return minors, full rank24, and the entire Cauchy factorization.
It checks 2,048 random selected-support ranks; the universal rank statements follow from that factorization and the Cauchy formula.
It does not exhaust the much larger sets of rectangular state-map minors.
There are 12,352 basis/scalar or random forward/adjoint comparisons against a literal dense24 matrix and its transpose, plus 128 four-packet comparisons.
All four tests in `gf8_mds_packet_tests.py` pass.

## Transfer-interface audit

The initial `gf8_feedthrough_mds_screen.py` uses the four-state construction.
Its continuation branch eliminates \(\min(z,w)\) input pivots before applying the remaining state-only equations.
Its state-root rank and charge, \(\max(0,z-w)\), agree with the Schur property above when below four.
Four or more roots exclude a live incoming state.

For an incoming coordinate-zero flat of rank \(k\), at most \(4-k\) incoming coordinates remain live.
The joint-return condition then requires \(z\le w-k-1\).
The screen uses this cutoff and charges the four return equations together with all selected output equations.
It does not multiply separate probability charges using the same current input or state diagonals.
The zero-to-zero branch also imposes the stronger necessary condition \(w>4+z\).

The local map assumptions match the screen, and its three exact tests pass on independent replay.
This audit relies on the separately established state-rank upper-tail and live-kernel lemmas for the probability step.
It does not upgrade floating-point search output to a complete directed certificate.
The corresponding eight-state transfer must replace the state dimensions consistently, including the root-tail helper and coordinate-flat counts.

The subsequent partial-live return helper removes the coordinate-flat enumeration from returns.
An independent audit of `PARTIAL_LIVE_RANK.md` and `partial_live_rank.py` found no soundness issue.
The dimensions \(d,j\) must refer to the joint coordinate image, not the larger source space: the hidden source kernel cancels from the rank gain.
For each fixed input-live vector, support-restricted Cauchy rank gives the required probability even when some state coordinates vanish.
All current diagonal randomness is included in that one event.
The pointwise lower bound on the count of input-live kernel vectors produces a valid one-mean linear program; feasible vertices and two-point mixtures suffice.
Maximizing over every permitted image dimension avoids an unproved monotonicity shortcut.

All four supplied partial-live tests passed on replay.
An additional independent GF(8) check used a full4 Cauchy packet with two input and two state coordinates.
It checked 198 moments across 22 joint image spaces, all \(7^4\) nonzero diagonals, zero/one/two selected outputs, and tilts 1, 3, and 7.
The image spaces included coordinate zeros and correlated low-rank input/state forms; every comparison passed.
This supports the local helper, not a directed whole-code certificate.
