# Construction

Fix a binary linear outer code

\[
 C:\mathbb F_2^{128}\longrightarrow\mathbb F_2^{256}
\]

with weight spectrum \(A_0,\ldots,A_{256}\). The current calculations use a
modeled complement-symmetric even spectrum with minimum distance 38. They do
not identify an explicit code with that spectrum.

The message contains \(2^{13}\) outer blocks. For each outer block \(j\),
setup samples an independent coordinate permutation \(\pi_j\in S_{256}\).
The encoder applies \(\pi_j\) to the corresponding outer-code output. It
then transposes the resulting \(2^{13}\)-by-256 bit matrix. The transpose
produces 256 regions of length \(2^{13}\). Setup samples one independent
permutation of the bit positions in every region. The encoder concatenates
the permuted regions.

## Checkpoint accumulator

The inner state is \(q\in\mathbb F_2^{64}\), initially zero. Number the
concatenated input bits by \(p\), starting at zero, and set
\(\ell(p):=p\bmod 64\). For input bit \(x_p\), the encoder computes

\[
 q_{\ell(p)}\gets q_{\ell(p)}+x_p,
 \qquad
 y_p\gets q_{\ell(p)}.
\]

The implementation batches 32 consecutive input bits. An epoch contains 32
batches, hence 1024 input bits and 16 visits to every state lane.

Identify the state space with \(\mathrm{GF}(2^{64})\) in a fixed polynomial
basis. For every boundary between consecutive epochs, setup samples an
independent

\[
 a_e\gets \mathrm{GF}(2^{64})^*.
\]

Let \(L_{a_e}\) denote multiplication by \(a_e\). At the boundary, the
forward encoder replaces the state by

\[
 q\gets L_{a_e}^{T}q.
\]

The terminal state is discarded. The random permutations and field
multipliers are sampled once during setup. Online encoding is deterministic.

## Transposed evaluator

The deployed evaluator acts on 128-bit blocks and traverses the inner chain
in reverse. At a checkpoint, it applies ordinary field multiplication
\(L_{a_e}\) to the 64 block-valued state lanes. Within an epoch, the
transpose of one accumulator update uses one block XOR:

\[
 z_p:=\bar y_p+\bar q_{\ell(p)},
 \qquad
 \bar x_p:=z_p,
 \qquad
 \bar q_{\ell(p)}:=z_p.
\]

The evaluator keeps values in 128-bit block form. It does not bitslice the
complete vector. A checkpoint implementation may pre-expand the linear map
or use a fixed-width carryless-multiplication circuit; that choice does not
change the construction.

## Geometry at the selected parameters

Every transposed region has length 8192 and contains eight complete epochs.
The full inner has 2048 epochs. An outer word of weight \(w\) selects a
uniform \(w\)-subset of the 256 regions. With one active outer block, every
selected region contains one impulse at a uniform position.

Changing the state width, checkpoint interval, accumulator rule, or
checkpoint map defines a different candidate.
