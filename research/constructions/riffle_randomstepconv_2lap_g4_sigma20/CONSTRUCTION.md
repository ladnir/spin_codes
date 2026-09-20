# Riffle RandomStepConv-2Lap g=4 sigma=20

This candidate removes the free starting boundary of the one-lap
RandomStepConv inner. The encoder traverses the same packet sequence twice.
The first traversal initializes the state, and the second traversal produces
the retained output.

## Parameters

The message contains (2^{20}) bits. Set

\[
g:=4,
\qquad
\sigma:=20,
\qquad
N:=524352.
\]

The outer encoder produces (N) packets in \(\mathbb F_2^g\). The inner
encoder retains one packet per position. Its binary output length remains

\[
gN=2097408.
\]

The first traversal increases encoding work but does not increase the codeword
length.

## Inherited outer construction

The candidate inherits the complete outer construction of **Riffle
RandomStepConv g=4 sigma=20**:

1. the systematic \([16386,16384,3]\) code over \(\mathbb F_{2^{64}}\),
   with the ShiftAlpha64 double-parity schedule;
2. the extended binary BCH \([128,64,22]\) encoder for every field symbol;
3. an independent uniform permutation of the 128 bits inside every BCH
   block;
4. division of each permuted BCH block into 32 four-bit packets; and
5. one uniform global permutation of all \(N\) packets.

Both traversals use the same globally permuted packet sequence. There is no
second packet permutation.

## Two-lap RandomStepConv inner

Let

\[
x=(x_1,\ldots,x_N)\in(\mathbb F_2^g)^N
\]

be the globally permuted packet sequence. Setup independently samples two
matrix families

\[
A_i,B_i\gets
\mathbb F_2^{(g+\sigma)\times(g+\sigma)}
\qquad (1\le i\le N).
\]

All matrices remain fixed for every encoding operation.

Initialize \(s_1:=0^\sigma\). The burn-in lap computes

\[
(u_i,s_{i+1}):=A_i(x_i,s_i)
\qquad (1\le i\le N).
\]

The encoder discards every \(u_i\) and carries \(s_{N+1}\) across the wrap.
Set \(t_1:=s_{N+1}\). The retained lap computes

\[
(y_i,t_{i+1}):=B_i(x_i,t_i)
\qquad (1\le i\le N).
\]

The codeword contributed by the inner encoder is

\[
y:=y_1\|\cdots\|y_N.
\]

The encoder discards \(t_{N+1}\). The two families \((A_i)_i\) and
\((B_i)_i\) are independent. In particular, the retained lap does not reuse
the burn-in maps.

An implementation need not compute the discarded packets \(u_i\). Sampling
only the \(\sigma\) state rows of each \(A_i\) gives the same distribution for
the retained output.

## Live and off steps

A step is **live** when its input-state pair is nonzero. At a live step, the
output packet and next state are independent and uniform:

\[
(y_i,t_{i+1})\gets
\mathbb F_2^g\times\mathbb F_2^\sigma.
\]

An active input packet therefore starts a new live episode, regardless of the
incoming state. A live episode turns off with probability

\[
q:=2^{-\sigma}
\]

after each live step. Once the state is zero, zero input packets produce zero
output until the next active input packet.

The burn-in lap changes the geometry of an off prefix. Before every retained
position, the encoder has traversed one complete copy of the packet sequence.
Thus a retained off interval belongs to a gap between active packets. It is
not an uncharged prefix before the first active packet.

## Probability space and claim scope

The code ensemble samples the local BCH permutations, the global packet
permutation, and both matrix families during setup. A minimum-distance
argument bounds the probability, over this setup, that any nonzero outer word
has low retained weight.

The construction is linear after setup. The matrix ensemble includes
rank-deficient maps. A distance proof must include zero output among the bad
events.

## Boundary of this candidate

Each of the following changes defines a distinct candidate:

- retaining outputs from both laps;
- reusing the burn-in matrices during the retained lap;
- resetting or randomizing the state at the wrap;
- retaining either terminal state;
- sampling invertible rather than unrestricted matrices;
- replacing the distinct state with prior output bits;
- changing \(g\), \(\sigma\), or the outer parity count; or
- applying another packet permutation between the laps.
