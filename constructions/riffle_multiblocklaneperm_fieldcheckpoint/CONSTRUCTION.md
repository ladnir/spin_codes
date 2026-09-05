# Riffle MultiBlockLanePerm-32 FieldCheckpoint

This candidate replaces each cyclic group shift in MultiBlockStripeFresh-32
with a uniform permutation of the 32 block lanes.  Setup samples an independent
lane permutation for every 32-block group and every outer coordinate.

The route retains 256 groups per transposed region.  After the lane
permutations, each lane receives an independent permutation of its 256 group
positions.  The two large permutation dimensions are therefore both 256.

If one group contains (k) active blocks, their lanes form a uniform
(k)-subset.  This law equals (k) iid lane choices conditioned on distinct
lanes.  The conditioning event has probability

\[
 p_k=\frac{(32)_k}{32^k}.
\]

The lane permutations are absorbed into the precomputed route.  They add no
encoder operation relative to the cyclic-shift candidate.

