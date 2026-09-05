# Construction

Fix a binary linear outer code

\[
C:\mathbb F_2^{128}\longrightarrow\mathbb F_2^{256}
\]

with weight spectrum \(A_0,\ldots,A_{256}\). The current experiments use a
real-valued BCH-like spectrum model with minimum distance 38. They do not
identify an explicit code with that spectrum.

For every outer block \(j\), setup samples an independent coordinate
permutation \(\pi_j\in S_{256}\). The encoder computes

\[
\widetilde c_j:=\pi_j(C(m_j)).
\]

It places the words \(\widetilde c_j\) in a matrix and transposes the matrix.
The 256 coordinate positions become 256 regions. Setup independently
permutes the bits inside every region. The encoder concatenates the regions
and partitions the result into \(t\)-bit inner inputs.

## Inner interface

The inner encoder has \(s\) state bits. Define \(d:=t+s\). At step \(i\),
the encoder combines the input and state into

\[
v_i:=(x_i,q_i)\in\mathbb F_2^d.
\]

It applies a setup-sampled linear map \(M_i\), then splits the result:

\[
(y_i,q_{i+1}):=M_i(v_i).
\]

The initial state is zero. The encoder discards the terminal state.

### FieldMulStepConv

Identify \(\mathbb F_2^d\) with \(\mathrm{GF}(2^d)\) in a polynomial basis.
For \(a\in\mathrm{GF}(2^d)\), let \(L_a\) denote multiplication by \(a\).
Setup samples independent coefficients \(a_i\) from the complete field and
defines

\[
M_i:=L_{a_i}^{T}.
\]

This orientation makes the deployed transposed step equal to ordinary field
multiplication \(L_{a_i}\).

Fix nonzero \(v\). The map \(a\mapsto L_a^T v\) is injective. Otherwise,
some nonzero \(a\) would satisfy \(L_a^T v=0\). Since \(L_a\) is surjective,
this equality would make \(v\) orthogonal to every vector. Hence \(v=0\), a
contradiction. The coefficient map is therefore bijective. Thus \(M_i(v)\)
is uniform in \(\mathbb F_2^d\).

The coefficient description contains \(d\) bits per step.

### ToeplitzStepConv

Setup samples an independent uniform \(d\)-by-\(d\) Toeplitz matrix \(T_i\)
at every step and defines \(M_i:=T_i\). A Toeplitz matrix contains \(2d-1\)
independent diagonal bits.

For every fixed nonzero \(v\), the map from the diagonal bits to \(T_i v\)
has rank \(d\). Therefore, \(T_i v\) is uniform in \(\mathbb F_2^d\). The
transpose of a Toeplitz matrix is Toeplitz, so the transposed evaluator keeps
the same convolution structure.

The coefficient description contains \(2d-1\) bits per step.

## Transposed evaluator

The PCG evaluates the transpose of the complete binary encoder on 128-bit
elements. It processes the inner steps in reverse order. At step \(i\), it
computes

\[
\begin{pmatrix}
\bar x_i\\
\bar q_i
\end{pmatrix}
=M_i^T
\begin{pmatrix}
\bar y_i\\
\bar q_{i+1}
\end{pmatrix}.
\]

The evaluator has two implementation strategies. The bitsliced strategy
physically transposes a padded 128-by-128 bit tile, applies the scalar map to
128 parallel lanes with vector carryless multiplication, and transposes the
tile back. The direct strategy leaves the 128-bit elements in place and
applies the map with block XORs. Three direct Toeplitz evaluators are
implemented: four-bit combination tables, contiguous AVX2 diagonal XORs,
and compile-time Karatsuba convolution. All three keep every value in
128-bit block form. Both layout strategies compute the same transposed
linear map. The physical tile transpose is an implementation choice, not a
step in the construction.

## Distance interface

For both inner variants,

\[
v=0\Longrightarrow M_i(v)=0,
\qquad
v\ne0\Longrightarrow M_i(v)\text{ is uniform in }\mathbb F_2^d.
\]

The random dense inner has the same one-vector transition law. Consequently,
the existing RandomStepConv weight enumerator applies to both structured
variants. The equality concerns each fixed outer word. A first-moment union
bound does not require independence between different outer words.

The maps and permutations are sampled once during setup. Online encoding is
deterministic.

## Claim scope

The current calculation covers one active outer block and uses a modeled
spectrum. It does not establish the distance of the complete construction.
The next proof step must sum all multi-block weight configurations induced by
the spectrum and the independent coordinate permutations.
