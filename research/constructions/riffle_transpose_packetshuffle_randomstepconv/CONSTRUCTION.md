# Riffle TransposePacketShuffle-RandomStepConv g

Let the message contain \(n\) bits. Fix an even outer block length \(B\), a
packet width \(g\), and an inner state width \(\sigma\). Define

\[
L:=\frac{2n}{B},
\qquad
P:=\frac{L}{g}.
\]

Setup independently samples \(L\) uniform rate-half linear injections

\[
C_j:\mathbb F_2^{B/2}\longrightarrow\mathbb F_2^B.
\]

Place their outputs in an \(L\)-by-\(B\) bit matrix. Partition the \(L\)
outer-block indices into \(P\) fixed groups of size \(g\). After transposing
the matrix, each bit coordinate is a row of \(L\) bits. In every row, the
bits from one fixed block group form one \(g\)-bit packet. Setup independently
samples a uniform permutation of the \(P\) packets in each of the \(B\) rows.
The encoder emits the permuted rows consecutively.

The one-lap inner encoder has state \(s_i\in\mathbb F_2^\sigma\). Setup
independently samples an unrestricted binary linear map

\[
M_i:\mathbb F_2^{g+\sigma}\longrightarrow
     \mathbb F_2^{g+\sigma}
\]

at every packet position. Starting from \(s_1=0\), the encoder computes

\[
(y_i,s_{i+1}):=M_i(x_i,s_i).
\]

The codeword is the concatenation of the packets \(y_i\); the final state is
discarded. The outer injections, row permutations, fixed block partition,
and inner maps remain fixed for all encodings.

The random outer injections are a constituent-code model, not an
implementation recommendation. The packet permutation work is \(B\)
permutations of \(P=L/g\) items.

## Current proof model

Condition on \(a\) active outer blocks. Replacing each uniform nonzero outer
word by independent uniform bits costs the exact factor

\[
(1-2^{-B})^{-a}.
\]

A fixed group containing \(r\) active blocks produces a zero packet with
probability \(2^{-r}\). Otherwise its averaged RandomStepConv transition is
the common nonzero-input transition.

The current evaluator packs the \(a\) active blocks into
\(\lfloor a/g\rfloor\) full groups and, when needed, one partial group. This
is conservative if merging two partially occupied groups cannot decrease the
probability of a low-weight inner output. Under the independent-bit
relaxation, a merged packet is active exactly when either source packet is
active. A coupling can place that merged activity at one of the two original
packet positions. Its time support is therefore contained in the split
support. The averaged RandomStepConv process is monotone under adding input
support: an added impulse can start a live episode but cannot remove output
from an existing episode.

`proof/PACKING_DOMINATION.md` proves the coupling, including the uniform
without-replacement packet order. The packed calculation is therefore a
valid analytic upper model. Complete outward-rounded occupation sums are
certified for the \(g=4\) cells in `TRADEOFF_G4.md` and the \(g=8\) cells in
`TRADEOFF_G8.md`.
