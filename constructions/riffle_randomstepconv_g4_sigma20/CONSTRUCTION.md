# Riffle RandomStepConv g=4 sigma=20

This candidate keeps the current packetized outer construction and replaces
the four parallel accumulators with one time-varying random convolution.  The
inner state couples all four bits of a packet.

## Parameters

The message contains (2^{20}) bits.  Write

\[
g:=4,
\qquad
\sigma:=20.
\]

The outer encoder produces (N=524352) packets in
\(\mathbb F_2^g\).  The inner encoder emits (N) packets, so its binary
output length is (gN=2097408).

## Inherited outer construction

The candidate inherits these mechanisms from **Riffle
ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4**:

1. the systematic \([16386,16384,3]\) field code over
   \(\mathbb F_{2^{64}}\), with the ShiftAlpha64 double-parity schedule;
2. the extended binary BCH \([128,64,22]\) encoder for every field symbol;
3. an independent uniform permutation of the 128 bits inside every BCH
   block;
4. division of each permuted BCH block into 32 four-bit packets; and
5. one uniform global permutation of all (N) packets.

The current candidate uses two field parity symbols.  The P4, P8, and P16
outer variants remain separate hypothetical constructions.

## RandomStepConv inner

Let (x_1,\ldots,x_N\in\mathbb F_2^g) be the globally permuted packets.
The encoder maintains a distinct state

\[
s_i\in\mathbb F_2^\sigma.
\]

Setup independently samples

\[
M_i\gets
\mathbb F_2^{(g+\sigma)\times(g+\sigma)}
\qquad (1\le i\le N).
\]

The sampled matrices remain fixed for every encoding operation.  Initialize
\(s_1:=0^\sigma\).  At step (i), compute

\[
(y_i,s_{i+1}):=M_i(x_i,s_i),
\]

where (y_i\in\mathbb F_2^g).  The inner output is

\[
y:=y_1\|\cdots\|y_N.
\]

The baseline discards (s_{N+1}).  Appending the terminal state changes the
termination rule and therefore defines a different candidate.

For every fixed nonzero (v\in\mathbb F_2^{g+\sigma}), the value (M_i v)
is uniform in \(\mathbb F_2^{g+\sigma}\) over the setup randomness of (M_i).
Consequently, a live step emits a uniform output packet and an independent
uniform next state.  A step is live exactly when \((x_i,s_i)\ne0\).

The construction is linear after setup.  The setup ensemble includes
rank-deficient matrices.  A minimum-distance proof must therefore include
zero output among the bad events; such a proof also establishes injectivity
on the outer code.

## Boundary of this candidate

The following changes require distinct names:

- storing the last \(\sigma\) output bits instead of a distinct state;
- sampling invertible rather than unrestricted matrices;
- reusing one time-invariant matrix;
- retaining the terminal state;
- changing (g) or \(\sigma\); or
- changing the outer parity count.

