# Riffle ParityShear-12 SplitState s=19

This variant targets 11 percent relative distance with a small change to the
outer code and the smallest state size that closes the current regular-class
calculation.

For each 256-bit outer block, setup samples a target coordinate (j) and a
12-element subset (S) of the other coordinates. The outer codeword (c)
is replaced by (c'), where

\[
c'_j=c_j+\sum_{u\in S}c_u,
\qquad
c'_v=c_v\quad(v\ne j).
\]

This linear shear is invertible because (j\notin S). It changes weight by
at most one. Because 12 is even, it fixes the all-one codeword. In the
transposed evaluator, the shear costs exactly 12 block XORs per outer block.

The inner is SplitState-PreAddMul-RM2Sub with step size 128 and state size
19. The selected image code has minimum distance 48. The corresponding
kernel has minimum distance six and contains no word of weight four.

At output length (2^{21}) and relative distance 11 percent, the modeled
regular-class margins are:

- 59.4859 bits for occupation one;
- 108.8061 bits for occupations 2 through 100; and
- 46.3294 bits for occupations 101 through 8192.

Thus the complete regular occupation range closes above 40 bits. State size
18 fails a medium occupation by about 7,568 bits. A ten-input shear fails the
all-active occupation by about 109 bits. Therefore (s=19) and shear size 12
are the minimum tested parameters that satisfy the current regular bound.

This is not yet an end-to-end theorem. It uses the modeled outer spectrum,
nearest binary64 arithmetic, and an expected spectrum over independently
sampled shears. The unchanged all-one class has large margins at the sampled
pure checkpoints. General mixtures with multiple all-one blocks remain a
separate proof gate.

See `PROOF_STATUS.md` and `PERFORMANCE_ESTIMATE.md` for the exact scope.
