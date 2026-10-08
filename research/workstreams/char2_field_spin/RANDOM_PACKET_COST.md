# Fresh Toeplitz packet maps: exact law and arithmetic cost

The inexpensive candidate uses a two-symbol state and fresh independent rectangular Toeplitz maps per packet.
It needs no random state refresh, input diagonal, or dense random feedthrough.
At width sixteen, its two maps cost twelve GFNI512 products per packet.
With the random Toeplitz[48,32] outer, the leading full arithmetic estimate is 153/64 GFNI512-equivalents per source symbol.
This is a construction and cost reference, not a whole-code distance certificate or a measured implementation.

## Packet interface and exact transition

Work over \(\mathbb F_q\), with \(q=16\) for the proposed implementation.
A packet receives \(X\in\mathbb F_q^t\) and state \(R\in\mathbb F_q^s\), initially zero.
For each packet, sample independent Toeplitz matrices
\(A\in\mathbb F_q^{t\times s}\) and \(C\in\mathbb F_q^{s\times t}\), with independent field elements on their diagonals.
All matrices are fresh across packets and independent of the fixed input word and preceding history.
The primary map is

\[
Y=X+AR,\qquad Q=R+CX.
\]

Emit \(Y\) and retain \(Q\) as the next state.
The state remains nonzero across an empty packet whenever its incoming value is nonzero.
Unlike a Schur update using \(CY\), this update does not allow an empty packet to kill the state.

For every nonzero fixed vector \(v\), a fresh rectangular Toeplitz matrix maps \(v\) uniformly onto its full output space.
To see this, take the first nonzero coordinate \(j\) of \(v\).
Select the consecutive diagonal variables containing the \(j\)-th column entries.
Their coefficient matrix is triangular with nonzero diagonal \(v_j\).
Conditioning on the remaining diagonal variables therefore leaves a uniform output.

Consequently, no uniform-incoming-state assumption is needed.
Condition on any preceding history and on the resulting fixed incoming state.
If \(R\ne0\), then \(AR\) is uniform; if \(X\ne0\), then \(CX\) is uniform.
The two quantities are independent because \(A\) and \(C\) are independent.

Let \(z\) mark output Hamming weight and put

\[
Z(z)=\left(\frac{1+(q-1)z}{q}\right)^t.
\]

Rows and columns of the following transfer matrices index zero and nonzero state, in that order.
For an empty packet, the exact transfer is

\[
T_0(z)=\begin{pmatrix}1&0\\0&Z(z)\end{pmatrix}.
\]

For a fixed nonempty packet of input weight \(w\), it is

\[
T_w(z)=
\begin{pmatrix}
q^{-s}z^w&(1-q^{-s})z^w\\
q^{-s}Z(z)&(1-q^{-s})Z(z)
\end{pmatrix}.
\]

The first row uses the actual input weight because the zero-state output is exactly \(X\).
Replacing \(z^w\) by the nonzero-uniform output enumerator without adding a suitable feedthrough would change the law.

The encoder is causally invertible even if the full packet matrix is singular.
Given the incoming state and \(Y\), recover \(X=Y+AR\), then update \(Q=R+CX\).
For the extended adjoint, with output dual \(Z_i\) and next-state dual \(H_{i+1}\), use

\[
X_i^*=Z_i+C_i^T H_{i+1},\qquad
H_i=H_{i+1}+A_i^T Z_i.
\]

Thus the transpose has the same two-map arithmetic cost.
Zero terminal dual gives the transpose of the encoder that discards its final state.

## Field-linearity for the consuming symbol type

All packet coefficients lie in GF16.
The maps therefore extend linearly to any compatible symbol field containing GF16, including the 8-, 32-, 64-, and 128-bit symbol representations.
The abstract GF16 coefficients must first be embedded into the byte basis used by GFNI.
No multiplication of arbitrary large-field data records is required.

If the entire encoder, including its outer code, is defined over GF16, its minimum symbol distance is unchanged by such a scalar extension.
Expand a nonzero extension-field message in a GF16 basis.
At least one nonzero component message produces a GF16 codeword whose support is contained in the extension-field output support.
Conversely, the GF16 codewords remain available after extension.
This observation justifies a GF16 word-counting proof for a GF16-defined encoder; it does not justify substituting a GF256 outer into that proof.

## Two-state circuit

For a two-by-two Toeplitz block, write

\[
M=\begin{pmatrix}a&b\\c&a\end{pmatrix}.
\]

With coefficient pre-additions prepared outside encoding, compute

\[
p=a(x_0+x_1),\quad q=(a+b)x_1,\quad r=(a+c)x_0,
\qquad Mx=(p+q,p+r).
\]

This uses three scalar GF16 products and three data XORs.
For even packet width \(t\) and state size two, split the emission map into \(t/2\) such row blocks.
Split the injection map into \(t/2\) column blocks and sum their two-coordinate outputs.
The blocks retain the overlaps prescribed by the original Toeplitz taps.
Each map uses \(3t/2\) scalar products, so both maps use \(3t\).

For \(t=16\), this is 48 scalar products, packed as twelve GFNI512 products on 128-bit records.
The inner arithmetic volume is therefore \(3/4\) GFNI512-equivalents per encoded record.
The output addition \(Y=X+AR\) and state addition \(Q=R+CX\) use XORs, not extra field products.
Lane interleaving, reduction of the injection contributions, and coefficient expansion remain real implementation costs.

The two maps contain seventeen independent GF16 diagonal values each.
Their eight prepared block triples contain 24 coefficient values per map.
Byte-sized prepared coefficients require 48 bytes per packet; packing raw nibbles could reduce storage but adds decoding work.
Expanding every scalar coefficient to sixteen bytes would instead consume 768 bytes per packet.
A native design should load compact coefficients and form GFNI operand vectors in registers, rather than silently assume that expanded traffic is free.

The state recurrence contains no multiplication:

\[
R_i=R_0+\sum_{j<i}C_jX_j.
\]

Its independent increments admit a short batched XOR prefix scan.
This also permits packing across consecutive packets without requiring separate independent streams.
The reverse adjoint has the same property.
These are scheduling opportunities, not a benchmarked speed claim.

## Optional three-product randomized birth map

If the proof needs uniform nonzero output at a zero-state birth, the two-symbol feedthrough can supply it cheaply.
For \(t=2\), sample

\[
B=\begin{pmatrix}a&b\\c&a\end{pmatrix}
\quad\text{subject to}\quad a^2+bc\ne0,\quad(b,c)\ne(0,0),
\]

uniformly from the accepted triples, independently of \(A,C\).
Use \(Y=BX+AR\), leaving \(Q=R+CX\) unchanged.
This family consists of the nonscalar invertible trace-zero matrices.
It has \((q-1)(q^2-1)\) elements and is invariant under conjugation by GL2.
For \(X=(1,0)^T\), each nonzero output \((a,c)^T\) has exactly \(q-1\) acceptable choices of \(b\).
Conjugation invariance gives the same count for every nonzero \(X\).
Thus \(BX\) is exactly uniform on the nonzero two-vectors.

At \(q=16\), the family has 3825 matrices and uniform rejection accepts 3825/4096 of all triples.
The three-product Toeplitz circuit applies without an additional fixed multiplier.
This is preferable here to assuming that a quadratic-extension multiplier has an XOR-only three-product realization over GF16.
The irreducible-quadratic reduction in the elementary extension-field circuit introduces another fixed multiplication.

With \(t=s=2\), all three maps together require nine scalar products per packet.
They fit two GFNI512 operations and one GFNI128 operation, excluding lane formation and recombination.
Alternatively, four consecutive packets use nine full GFNI512 operations after the state prefix scan.
The arithmetic volume is \(9/8\) GFNI512-equivalents per encoded record.
The birth row of the transfer matrix then replaces \(z^w\) by
\(\big((1+(q-1)z)^2-1\big)/(q^2-1)\).
Nothing else in the transition law changes.

## Full source-normalized comparison

Let an outer code have length \(b\), dimension \(k\), and \(K=kL\) source symbols.
The inner processes \(N=bL\) encoded symbols.
Counts below include both inner maps and the outer transpose, with no input diagonal or state refresh.
One GFNI512-equivalent is four scalar products on 128-bit records.

| Outer, GF16 coefficients | Inner | Inner products/source | Outer products/source | Total/source | Minimum record bytes/source |
|---|---|---:|---:|---:|---:|
| Random Toeplitz[48,32] | identity feedthrough, width16/state2 | 1.12500 | 1.265625 | 2.390625 | 88.000 |
| Random Toeplitz[32,24] | identity feedthrough, width16/state2 | 1.00000 | 0.843750 | 1.843750 | 80.000 |
| Fixed GRS[16,12] | identity feedthrough, width16/state2 | 1.00000 | 0.145833 | 1.145833 | 80.000 |
| Fixed GRS[16,12] | randomized birth, width2/state2 | 1.50000 | 0.145833 | 1.645833 | 80.000 |

These rows are cost alternatives, not statements that their outer distances or full-code proof bounds are adequate.
The identity-feedthrough inner has the same leading product volume for every even packet width when the state has two symbols.
Changing packet width still changes the transition enumerator, packet-control costs, state reductions, and SIMD utilization.

The Toeplitz[48,32] parity transpose contains two sixteen-by-sixteen Toeplitz blocks.
The existing binary-recursive circuit uses \(2\cdot3^4=162\) GFNI512 products per four source rows.
For Toeplitz[32,24], three eight-by-eight blocks use \(3\cdot3^3=81\).
The verified fixed GRS[16,12] parity circuit uses seven products per four source rows.
Hence the first table total is exactly \(153/64\); the remaining totals are \(59/32\), \(55/48\), and \(79/48\).

For context, the earlier GF256 regional GRS[34,24]/state4 estimate was \(385/192\), approximately 2.00521.
The Toeplitz[48,32] alternative has about 19.2% more field-product volume, despite its cheaper inner.
That comparison does not transfer the GF16 proof law to the GF256 construction.
The outer circuit now dominates the Toeplitz[48,32] arithmetic budget.

Minimum logical record traffic is \(16(3N+K)\) bytes for inner input, plane write/read, and final output.
A regional temporary-buffer schedule adds \(32N\) bytes.
Routing indices, matrix coefficients, cache write allocation, and lane rearrangement add traffic or instructions beyond that minimum.
No measured elapsed-time or bandwidth estimate follows from the product counts alone.

## Cost of the separate Toeplitz[64,40], five-state backup

For width16/state5, the elementary four-plus-one split costs 52 scalar products per rectangular inner map.
The elementary outer parity-transpose decomposition costs 351 GFNI512 products per four source rows.
These give \(351/160+2.6=4.79375\) GFNI512-equivalents per source symbol at rate5/8.
Two explicit XOR-only reorganizations reduce that upper count, but not enough to make it competitive with the regional estimate.

For the forty-by-twenty-four outer transpose, split into eight-symbol blocks and reverse the three input blocks.
The resulting five-by-three block matrix is Hankel.
Its first three block rows form a symmetric three-by-three matrix and need six block products.
The last two rows need a symmetric two-by-two product plus two products for the remaining column: five block products.
All eleven products are eight-by-eight Toeplitz products, including after coefficient pre-additions.
The outer therefore uses \(11\cdot27=297\) GFNI512 products per four source rows.
Transposing this circuit gives the same count for the twenty-four-by-forty direction.

For either inner map, split the dimension sixteen as five plus five plus five plus one.
Reversing a five-symbol input block converts its square Toeplitz block to a symmetric Hankel block.
A symmetric five-by-five matrix-vector product uses fifteen scalar products: one for each diagonal correction and one shared product for each off-diagonal pair.
The three square blocks and remaining row or column therefore use \(3\cdot15+5=50\) scalar products per map.
Both maps use 100 scalar products per packet, or 25 GFNI512-equivalents.
Realizing exactly 25 full-vector instructions requires packing the two maps' partial tails together; otherwise the instruction count can remain26.

The improved leading arithmetic volume is

\[
\frac{297}{160}+\frac{25}{16}\frac{64}{40}
=4.35625
\]

per source symbol, approximately 2.166 times the regional \(177/88\) estimate.
It excludes any additional input diagonal, coefficient expansion, shuffle, and routing costs.
A read-only scalar check compared both reorganizations with dense GF16 Toeplitz multiplication on 64 independent random examples per shape and verified their product counts.
These are modest circuit improvements, not a credible current reason to start a native backup implementation.

## Exact checks and next gate

The scalar reference exhausts all 255 nonzero two-vectors and all 4096 Toeplitz triples.
It checks the three-product circuit, uniform images of the full Toeplitz family, and uniform nonzero images of the 3825-element invertible subfamily.
It checks the width16 emission tap-map rank for all nonzero two-vectors, plus direct rectangular-map comparisons at widths 2, 4, 8, 16, and 32.
Variable-packet stream checks cover causal inversion and the extended adjoint with nonzero initial and terminal states.
The exact source-normalized cost identities are checked too.

```sh
python -B research/workstreams/char2_field_spin/theory/random_packet_cost_reference.py --check
```

All checks pass.
The next gate is a whole-code proof using the stated fresh-matrix law and an outer whose arithmetic cost is acceptable.
No native kernel, remote execution, or benchmark was added.
