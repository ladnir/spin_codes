# Riffle S-Stripe RandomStepConv g=8

The global packet permutation permits a dense outer word to move almost
entirely into a late suffix. This family replaces that permutation with a
region-balanced permutation. The packet width remains eight bits.

## Parameters

Let the message length be \(n\). Setup independently samples each outer block
uniformly from the set of linear injections

\[
C_i:\mathbb F_2^{B/2}\longrightarrow\mathbb F_2^B.
\]

Define

\[
L:=\frac{2n}{B},
\qquad
q:=\frac{B}{g},
\qquad
g:=8.
\]

Thus there are \(L\) outer blocks and \(q\) packets in each encoded block.
The construction parameter \(S\) must divide \(q\).

## S-stripe permutation

Write the encoded packets as an \(L\)-by-\(q\) matrix. Rows correspond to
outer blocks, and columns correspond to packet coordinates within a block.

Partition the \(q\) columns into \(S\) classes of equal size. For each class
\(r\), setup samples an independent uniform permutation of the \(Lq/S\)
packets in that class. The encoder emits the \(S\) permuted classes as
consecutive regions.

Consequently, every outer block supplies exactly \(q/S\) candidate packets
to every region. The permutation does not guarantee that these packets are
nonzero.

The following settings name distinct instances:

- \(S=1\): the global-permutation control;
- \(S=2\): the half-stripe instance;
- \(S=4\): the four-stripe instance; and
- \(S=q\): the full-stripe instance, with one packet coordinate per region.

Every region uses fresh independent permutation randomness. Reusing one
permutation across regions defines a different construction.

## RandomStepConv inner

Let \(x_1,\ldots,x_N\in\mathbb F_2^g\) denote the emitted packet sequence,
where \(N=Lq=2n/g\). The inner encoder maintains
\(s_i\in\mathbb F_2^\sigma\) and initializes \(s_1:=0^\sigma\).

Setup independently samples

\[
M_i\gets\mathbb F_2^{(g+\sigma)\times(g+\sigma)}
\qquad (1\le i\le N).
\]

At position \(i\), the encoder computes

\[
(y_i,s_{i+1}):=M_i(x_i,s_i).
\]

The codeword is \(y_1\|\cdots\|y_N\). The encoder discards the final state.

## Ensemble and scope

Setup samples the outer injections, the \(S\) region permutations, and all
inner maps. These values remain fixed for all encodings. The present random
outer family is a proof model, not an implementation recommendation.

For a deterministic outer constituent, ordinary weight spectra do not by
themselves control weight in every region. Such a replacement requires either
a within-block coordinate permutation or a regional weight profile.

The current objective is a first-moment bound at relative distance 9% for a
\(2^{20}\)-bit message. A numerical diagnostic is not a minimum-distance
proof until it covers every outer occupation and validates floating-point
rounding.
