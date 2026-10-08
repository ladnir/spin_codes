# A cube circuit for the shared eight-point map

The shared FFT circuit can be replaced by an exact circuit with fewer field
products. For shift 64, the proposed schedule uses 59 packed GFNI arithmetic
calls per 64-record packet, including the genuine-field state refresh. The
corresponding FFT schedule uses 82. The reduction is 28%, not a measured speedup.
The later [native comparison](HD_BLOCKWISE_SHARED_NATIVE.md) measures7.65ms
for this circuit versus6.54ms for the same-map FFT. A63-product whole-vector
rewrite takes7.00ms. The arithmetic reduction does not translate into a win.

The useful reduction comes from sharing Boolean cube sums, not merely from
factoring a rank-four matrix. A rank-four tile residual does **not** determine
the existing four feedback coordinates. The combined map has rank six.

## Map and field conventions

All fixed coefficients below lie in AES GF256. Each 128-bit record contains
eight GF65536 symbols, using the existing quadratic extension of GF256. A fixed
GF256 coefficient scales both bytes of every symbol. Four independently scaled
records therefore fit in one 512-bit GFNI multiplication.

Let `b_j(x)` be the global normalized novel basis from
[the shared-map note](HD_BLOCKWISE_SHARED_MAPS.md), for `0 <= j < 8`.
Write `E_alpha[u,j]=b_j(alpha XOR u)`. A packet contains eight input tiles
`X_alpha`, with `alpha=0,8,...,56`, and incoming state `Q` of length four.
For a fixed shift `beta` outside `0..7`, the reference circuit is

\[
 z_\alpha=E_\alpha^{-1}X_\alpha,\qquad
 Y_\alpha=E_{\alpha+\beta}(z_\alpha+(Q,0)),\qquad
 h=\sum_\alpha z_\alpha[4:8].
\]

Here field addition on byte labels is XOR. After emitting all tiles, the next
state is `M Q+h`, with the same genuine-GF65536 rank-one mixer as before.
Every tile uses the incoming state, not a partially updated state.

The working choice is `beta=64`. Its output nodes `64..127` are disjoint from
the packet's input nodes `0..63`. This removes the particular adjacent-tile
overlap of shift 8. It is not a whole-code distance certificate.

The feedthrough matrix of every tile is the same:

\[
 F_\beta=E_\beta E_0^{-1}.
\]

The circuit below preserves this map, the shifted state-emission map, and the
original feedback exactly. No input conjugacy or support-changing relabeling
is assumed.

## What the rank-four factorization does and does not give

Set `z=x^2+beta*x`. Every polynomial `p` of degree below eight has a unique form

\[
 p(x)=a(z)+x b(z),\qquad \deg a,\deg b<4.
\]

The eight polynomials `1,z,z^2,z^3,x,xz,xz^2,xz^3` have distinct leading degrees.
They therefore form a basis. Translation by `beta` fixes `z`, so

\[
 p(x+\beta)+p(x)=\beta b(z).
\]

Let `V` recover the four coefficients of `b` from the input evaluations.
Let `U` evaluate `beta*b(z)` at `x=0..7`. Then

\[
 F_\beta=I+UV,\qquad VU=0,\qquad
 \operatorname{rank}U=\operatorname{rank}V=4.
\]

The map `x -> x^2+beta*x` is injective on `0..7`: its only possible nonzero
collision difference is `beta`, which is outside that additive subspace.
Consequently the fixed vectors form an MDS `[8,4,5]` evaluation code.

There is also a four-coordinate pair factor. Write

\[
 F_\beta-I=\begin{bmatrix}E&D\\D&E\end{bmatrix},\qquad B=D^{-1}E.
\]

The matrices are XOR-circulant, and direct algebra gives `B^2=I`. Thus one can use

\[
 U=\begin{bmatrix}D\\BD\end{bmatrix},\qquad V=\begin{bmatrix}B&I\end{bmatrix}.
\]

For input halves `(x_0,x_1)`, compute `r=B*x_0+x_1`, then `w=D*r`, and emit
`(x_0+w,x_1+B*w)`. This is an exact low-rank implementation interface.
It is not automatically cheaper: it contains two four-point transforms by `B`
and one by `D`, before accounting for feedback and state emission.

The existing tile feedback is `C_alpha=E_alpha^{-1}[4:8,:]`.
Its kernel consists of evaluations of polynomials of degree below four.
The intersection with the fixed-polynomial space is `span{1,z}`, of dimension
two. It follows that

\[
 \operatorname{rank}\begin{bmatrix}V\\C_\alpha\end{bmatrix}=6.
\]

Thus two additional independent coordinates are needed beyond `V X` to recover
feedback. The emission space also intersects the fixed space in dimension two.
Neither all feedback nor all state emission can be folded into the same
four-coordinate residual without changing the map.

## Feedthrough from cube differences

Index a tile by three bits. Let `P_i` permute its coordinates by XOR with `2^i`,
and let `epsilon_i=I+P_i`. These operators commute and satisfy `epsilon_i^2=0`.
Since `F_beta` is XOR-circulant, it has the expansion

\[
 F_\beta=I+\sum_{\varnothing\ne S\subseteq\{0,1,2\}}
 a_S\prod_{i\in S}\epsilon_i.
\]

If `k_j=F_beta[0,j]`, the exact coefficients are

\[
 a_S=\sum_{j:\,j\supseteq S}k_j.
\]

For `beta=64`, the coefficients in mask order `0..7` are
`[1,142,4,62,188,137,192,80]`. The tile's first row is
`[16,105,170,110,165,217,144,80]`.

The vector `epsilon_S X` repeats each distinct value `2^|S|` times. Compute only
the values whose coordinate has zero bits on `S`. There are

\[
 3\cdot4+3\cdot2+1=19
\]

such cube sums. Multiply each by its corresponding `a_S`, then expand and XOR
the products into the original eight inputs. Three groups of four products use
three GFNI calls. The remaining seven products use two more calls. This gives
five GFNI calls per tile, or 40 per packet. Packing is local to each tile;
the count does not depend on buffering several packets.

All 19 cube sums are retained long enough to supply feedback. They contain
more information than the rank-four residual `V X`.

## Feedback from the same sums

For a six-bit mask `T` of size at most two, define a packet moment

\[
 m_T=\sum_{x:\,x\supseteq T}X_x.
\]

There are 22 moments. Each is an XOR of cube sums already computed for
feedthrough. Split `T` into its low three bits `T_l` and high three bits `T_h`.
For every tile index containing `T_h`, use the cube sum with difference mask
`7 XOR T_l` and representative coordinate `T_l`. The difference mask is nonzero
because `|T_l|<=2`.

The first four novel basis polynomials are

\[
 b_0(x)=1,\quad b_1(x)=x,\quad
 b_2(x)=(x^2+x)/6,\quad b_3(x)=x b_2(x).
\]

As functions of the six coordinate bits, each has Boolean degree at most two.
Their multilinear coefficients convert the 22 moments into

\[
 \mu_j=\sum_x b_j(x)X_x,\qquad j=0,1,2,3.
\]

The feedback is then the XOR-only conversion

\[
 h=D_0\mu,\qquad
 D_0=\begin{bmatrix}
 1&1&1&1\\1&0&1&0\\1&1&0&0\\1&0&0&0
 \end{bmatrix}.
\]

For every tile, the exact identity is `C_alpha=D_0 E_alpha[:,0:4]^T`.
The helper checks it at all eight offsets.

Ignoring zero and unit coefficients, the three nonconstant moment conversions
need respectively 5, 4, and 19 record products. Their total is 28, packed into
seven GFNI calls. This packing can mix products destined for different moments;
the subsequent XOR reductions must keep those destinations separate.

The local feedthrough still uses each tile's cube sums separately. Replacing
those sums by the four global feedback values would change the map to a
different construction and could restore the quiet-subcode defect.

## State emission and combined count

Expand the four functions `b_j(beta XOR x)` in their multilinear Boolean basis.
For each coefficient mask, multiply its constants by `Q_j` and XOR the results.
Reuse identical `(j,constant)` products across masks. Evaluate the resulting
degree-two Boolean polynomial at all 64 coordinates and add it to feedthrough.

For `beta=64`, the numbers of distinct nontrivial products involving
`Q_1,Q_2,Q_3` are 6, 5, and 21. The total is 32 record products, or eight packed
GFNI calls. `Q_0` contributes only an XOR. A recursive Boolean evaluation needs
at most 120 record XORs after the coefficient values are available; it needs no
further field products.

| Per 64-record packet | Packed GFNI arithmetic calls |
| --- | ---: |
| Eight local feedthroughs | 40 |
| Global feedback conversion | 7 |
| State-emission coefficient products | 8 |
| Genuine-GF65536 rank-one refresh | 4 |
| **Total, shift 64** | **59** |
| Reference shared FFT circuit, shift 64 | 82 |

Shift 8 gives 58 calls, but has the documented adjacent-tile overlap. Shift 128
gives 59, like shift 64. The working choice remains shift 64.

This count includes the entire packet map, not only its feedthrough. It excludes
loads, coefficient expansion, packing, shuffles, XORs, and spills. In particular,
22 moment records occupy 352 bytes. Keeping them live with emission products,
tile values, and constants may cause register pressure. The cube expansion also
needs more cross-record shuffles than the FFT circuit. A native improvement
cannot be inferred from the arithmetic count alone.

## Validation and next decision

The scalar prototype is
[`theory/hd_blockwise_shared_factor.py`](theory/hd_blockwise_shared_factor.py).
Its seven tests pass:

- both rank-four factorizations and the rank-six feedback interface;
- all 12,869 square minors of the shift-64 tile;
- every tile input basis vector;
- exact cube moments and feedback identities;
- all 68 input/state basis vectors of the complete packet;
- nonzero random states and both components of the GF65536 extension;
- the product counts and Boolean-degree restrictions.

The implementation-independent matrix checks preserve branch number nine.
They do not establish a full-code distance bound or failure probability.

The scalar checkpoint preceded native implementation. The later timing note
records both register-level schedules and their regression. Preserve this
factorization as an exact-map reference; do not continue tuning it as the
selected path to5ms. A structural change to the state and feedback maps may
permit more sharing, but requires a new proof interface and obstruction checks.
