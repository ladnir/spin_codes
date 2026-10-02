# Follow-up construction screens

The completed result is the GL15 construction in
[`PAIRED_S15_CONSTRUCTION.md`](PAIRED_S15_CONSTRUCTION.md), with a
62.04117326862128-bit whole-code setup-failure bound. The experiments below
are alternatives, not premises of that certificate. No further full replay
was started after the GL15 implementation met the timing target.

## Fourteen state bits with paired quadratic rows

Simply deleting two rows from the certified sixteen-bit map did not close
the current bound. A redesigned fourteen-bit map instead retains the seven
affine functions 1,z0,...,z5 and uses seven pairs of quadratic monomials.
Each quadratic monomial occurs in at most one pair. This retains copied
expansion coefficients and requires seven XORs to combine feedback channels.

For every candidate below, output coordinate z is the six-bit integer in
least-significant-bit order. Define A by evaluating the fourteen functions,
set C=A^T, and use a fresh independent uniform GL(14,2) update per physical
64-bit step. The outer code and routing are unchanged from the completed
paired-map construction.

`search_seven_pairs.py` ran a bounded search over 30,000 pairings. Several
had no rank-two quadratic combination. The fixed candidate in
`seven_pair_maps.py` has a fresh q=1 bound of 61.4418304126 bits. Its refined
q=3,...,128 floating estimates sum to only a 37.962-bit bound, so this
candidate did not justify a full replay.

A second bounded search rebuilt the actual local profiles for sixteen
pairings with that same spectrum. Candidate12 in
`seven-pair-profile-search.json` gave the strongest tested critical score.
Its quadratic functions, after 1,z0,...,z5, are

\[
\begin{aligned}
g_7&=z_0z_5+z_1z_2,&g_8&=z_0z_3+z_1z_5,\\
g_9&=z_1z_3+z_2z_4,&g_{10}&=z_1z_4+z_2z_5,\\
g_{11}&=z_0z_1+z_3z_4,&g_{12}&=z_0z_4+z_3z_5,\\
g_{13}&=z_0z_2+z_4z_5.
\end{aligned}
\]

Fresh enumeration of all 16,384 states gives weights
0,24,28,32,36,40,64 with respective counts
1,1,072,3,840,6,558,3,840,1,072,1. The map has full rank, CA=0, and rank-four
restriction on every consecutive four-coordinate packet. Its map identity is
`ac71c1cd3ab7b19f9a4eba67680ff1ae380eb92ccb1e4e254a4dc1d62bfec408`.

The search recomputed each candidate's local operators with outward
arithmetic, then evaluated the following global objectives in floating point:

| Occupancy q | 36 | 38 | 40 | 42 | 44 |
|---|---:|---:|---:|---:|---:|
| Estimated margin, bits | 47.3744 | 47.3280 | 48.2788 | 50.5505 | 53.9369 |
| Rational tilt | .1325 | .14 | .1475 | .155 | .1625 |

These five objectives do not give a whole-code bound. In particular, the
61.4418-bit q=1 result belongs to the earlier pairing, not candidate12.
Equal expansion spectra do not imply equal local profiles or certificates.

The next optional construction step is to freeze candidate12 in a separate
map module and screen its q=1, q=2, and all middle occupancies. A fresh full
replay is warranted only if those screens leave room for the union bound
and its implementation offers a measured gain over certified GL15.

## Less frequent or weaker state refresh

The skipped-refresh screens keep expansion and feedback at every physical
step but replace selected refresh matrices by the identity. The exact
identity specialization has lazy weight one: it never creates a fresh
uniform nonzero state. The screen's empty-input operator preserves exact
expansion-weight classes. Its other operators use fixed-state cancellation
and fiber bounds. They do not reuse a uniform-refresh return probability.

Alternating GL16 and identity updates gave a fresh q=1 bound of 62.3354 bits,
but the current middle-occupancy envelope was too weak. Replacing only every
fourth refresh by identity also failed this envelope. These screens provide
no full certificate and do not establish that either code has low distance.

The transvection screen uses a different update distribution. Sample
u uniformly from F2^16 minus zero, and v uniformly from u-perp, including
v=0. Set R=I+u*v^T. Since v^T*u=0, R is invertible. For every fixed nonzero
a, the exact law of R*a is one half the atom at a plus one half the uniform
distribution on nonzero vectors. A composition of r independent updates
has lazy weight 2^-r. Unconstrained rank-one updates do not have this law.

For r=1 the fresh q=1 margin was only 20.8381 bits. For r=2 it reached
40.3191 bits, but the middle-occupancy estimates were inadequate: q=3 gave
31.5678 bits and q=8 gave -138.8014 bits on the tested tilt grid. The current
proof therefore does not justify either update family as a replacement.
Further work would need a stronger state envelope, not a numerical transfer
from the complete GL15 or GL16 receipts.

Reproducer modules are `screen_refresh_stride.py`,
`screen_refresh_period.py`, and `screen_transvection.py`. Their outputs are
explicitly marked `whole_code_certificate=false`.
