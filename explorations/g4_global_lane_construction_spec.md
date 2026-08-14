# Global-lane construction for the g=4 rerun

## Purpose

The prior construction chose a punctured data lane independently in each
selected band-zero tile.  That rule does not provide one conditioned row
shared by the three sloped tile partitions.  The global-lane rule supplies the
required matching.  It also introduces correlation between puncture lanes, so
the prior g=4 result does not transfer without a complete outward rerun.

## Block and tile indices

Let

\[
I:=\mathbb Z_{256}\times\{0,\ldots,63\}
\]

index the 16,384 data-message blocks.  Write a block as \((t,\ell)\).
For bands \(b\in\{0,1,2\}\), fix slopes

\[
(s_0,s_1,s_2):=(0,9,20)
\]

and define the incident tile

\[
\tau_b(t,\ell):=t+s_b\ell\pmod {256}.
\]

For a fixed lane \(L\), define

\[
M_L:=\{(t,L):t\in\mathbb Z_{256}\}.
\]

For every band \(b\), the map \(t\mapsto t+s_bL\) is a translation of
\(\mathbb Z_{256}\).  Therefore, \(M_L\) meets every tile in band \(b\)
exactly once.  The 64 sets \((M_L)_{L=0}^{63}\) are disjoint and partition
\(I\).

The script `scripts/certify_global_lane_puncture_layout.py` verifies these
statements by exact enumeration.  It also certifies the coordinate marginal,
the graph-image distribution, and the public sampling interface.  The
committed script digest is
`2e0f70cde90f94974f4bd86931905d975eaea8a8c1bc727ab5ba0fa1f810aaf2`.

## Sampling interface

The construction samples its puncture locations as follows.

1. Sample one lane \(L^*\) uniformly from \(\{0,\ldots,63\}\).
2. Sample a uniform 128-subset \(H\subseteq\mathbb Z_{256}\) of band-zero tiles.
3. For each \(u\in H\), select block \((u,L^*)\).
4. In each selected block, sample one of its 42 band-zero coordinates uniformly and puncture that coordinate.
5. Map the 128 graph coordinates bijectively to the 128 punctured positions using the construction's graph-placement randomness.

The rule in step 3 uses \(s_0=0\).  Thus block \((u,L^*)\) is the unique
member of matching \(M_{L^*}\) incident to band-zero tile \(u\).

For a fixed block \((t,\ell)\), the probability that step 3 selects the block is

\[
\Pr[L^*=\ell]\Pr[t\in H]
=\frac1{64}\frac{128}{256}
=\frac1{128}.
\]

This equality preserves the prior per-block selection marginal.  It does not
make selections independent.  Conditioned on \(L^*\), every selected block
has lane \(L^*\).

## Proof interface

The one-conditioned-row outer may condition the matching \(M_{L^*}\).  Each
band tile contains exactly one conditioned block.  Every punctured block also
belongs to that matching.  This alignment permits the exact puncture and graph
calculation after an outward implementation binds the global-lane rule.

The matching certificate proves only incidence and marginal statements.  A
proof step that uses independent lane choices across tiles requires a new
argument.  In particular, the old independent-lane puncture adjustments do
not transfer from their one-block marginals.

The selected total-weight witness uses a separate global-lane argument.
`explorations/g4_global_lane_exact_graph_puncture_total_weight_lemma.md`
conditions on the lane and selected tiles.  It then applies exact block
factorization, the exponential bound, and Jensen's inequality.  That lemma
does not assume independent lane choices.

The global-lane rule does not by itself certify a two-row argument.  A future
second-row construction must state its own sampling rule and prove its
interaction with all punctured blocks.

## Rerun status

`G4_GLOBAL_LANE_RERUN_MANIFEST.json` supersedes the old savepoint as the
construction binding for the next run.  It reuses the prior ledgers and BSP
batches only as candidate proof inputs.  The old final report remains an
integrity reference and is not a certificate for this construction.

Before a full rerun, the reproduction gate must verify the construction spec,
matching certificate, proof inputs, and proof-code digests.  The full
certifier must then replay every used witness under the global-lane semantics.
