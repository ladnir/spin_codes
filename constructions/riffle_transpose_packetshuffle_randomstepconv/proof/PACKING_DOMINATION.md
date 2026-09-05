# Packing domination

This note proves the packing step used by the packet-version evaluator. The
claim applies after replacing each active outer word by independent uniform
bits. The existing conditioning factor accounts for that replacement.

Fix one transposed row. Suppose a fixed packet group contains \(r\) active
outer blocks. Its packet is nonzero with probability

\[
p_r:=1-2^{-r}.
\]

Packets from different groups use disjoint outer bits. Their nonzero
indicators are therefore independent.

## One packing move

Consider two group occupancies \(r\) and \(s\), where

\[
g>r\geq s>0.
\]

A packing move changes \((r,s)\) to \((r+1,s-1)\). Let \(K\) and \(K'\)
denote the numbers of nonzero packets before and after the move. Both random
variables take values in \(\{0,1,2\}\).

Set \(u:=2^{-r}\) and \(v:=2^{-s}\). The move replaces \((u,v)\) by
\((u/2,2v)\). Hence

\[
\Pr[K=0]=uv=\Pr[K'=0].
\]

The probability of two nonzero packets decreases because

\[
\begin{aligned}
\Pr[K=2]-\Pr[K'=2]
&=(1-u)(1-v)-(1-u/2)(1-2v)\\
&=v-u/2\\
&>0.
\end{aligned}
\]

Thus \(K'\) is stochastically dominated by \(K\). Adding the independent
nonzero-packet count from all other groups preserves this domination.

Repeated packing moves transform every occupancy vector of total weight
\(a\) into

\[
(\underbrace{g,\ldots,g}_{\lfloor a/g\rfloor},
  a\bmod g,0,\ldots,0),
\]

where the remainder entry is omitted when it is zero. Therefore the packed
occupancy minimizes the number of nonzero packets in stochastic order.

## From counts to input support

The row permutation is independent and uniform. Conditioned on exactly
\(h\) nonzero packets, their positions form a uniform \(h\)-subset of the
row's \(P=L/g\) packet positions.

Let \(H'\leq H\) be a coupling of the packed and original nonzero-packet
counts. Sample a uniform \(H\)-subset \(S\). Then sample a uniform
\(H'\)-subset \(S'\) of \(S\). The marginal distribution of \(S'\) is uniform
among all \(H'\)-subsets, and \(S'\subseteq S\). Applying this coupling
independently to all transposed rows couples the complete packed input support
as a subset of the original support.

## Monotonicity of RandomStepConv

It remains to compare two averaged inner executions whose input supports
satisfy \(S'\subseteq S\). Record only whether the inner state is zero. At a
position where both the input and state are zero, the output and next state
are zero. Otherwise, the random step map makes the output uniform in
\(\mathbb F_2^g\) and the next state uniform in \(\mathbb F_2^\sigma\). These
two values are independent.

Couple the two executions inductively. When both executions are active, give
them the same output and next-state zero indicator. When only the larger
execution is active, sample its output and next state independently. The live
state of the packed execution is always bounded by the live state of the
original execution. Its output weight is also bounded at every position.
Consequently, its total output weight is no larger under this coupling.

For every integer threshold \(d\),

\[
\Pr[W_{\mathrm{original}}\leq d]
\leq
\Pr[W_{\mathrm{packed}}\leq d].
\]

The probability is over the relaxed outer bits, row permutations, and inner
maps. The comparison concerns one fixed nonzero message. An arbitrary
coupling suffices because the first-moment argument uses only each message's
marginal bad-output probability.

Therefore the packed occupancy gives a valid upper bound for every active
outer-block set of size \(a\).

