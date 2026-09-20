# Performance result

The integrated state-size-19 ParityShear-12 implementation passes the dense
inner and independent staged end-to-end correctness oracles. On Peach, a
21-trial run pinned to CPU 15 has an end-to-end median of 10.747357 ms. The
materialized reverse inner has a median of 3.032319 ms.

The prior state-size-16 integrated baseline has a 10.821225 ms median. The
new median is 0.073868 ms lower. Because the measurements came from
different sessions, this comparison establishes no observed end-to-end
regression; it is not a speedup claim.

The complete timing samples, construction law, correctness oracles, compiler
flags, and source hashes are recorded in
`benchmark/parityshear12_s19_integrated_peach_7950x.json`.

## Earlier estimate

The existing integrated (s=16) implementation has a 10.821225 ms median
on Peach. Its reverse inner takes 2.596211 ms.

The preliminary (s=19) operation ledger has the following interior-epoch
costs in 128-bit block XOR equivalents.

| Component | s=16 | s=19 estimate | Increase |
|---|---:|---:|---:|
| Grouped (A) | 499 | 611 | 112 |
| Pruned RM stage for (B) | 346 | 346 | 0 |
| Quadratic projection for (B) | 42 | 55 | 13 |
| Field multiplication | 92 | 124 | 32 |
| Syndrome injection | 16 | 19 | 3 |
| Total inner epoch | 995 | 1155 | 160 |

The grouped-(A) value comes from a 30,000-candidate partition search. It is
not yet an optimality claim. The exact generated Paar circuit uses 525 XORs
for (A), but the analogous (s=16) circuit was slower because it spilled
registers. Grouped tables remain the preferred implementation starting
point.

ParityShear-12 adds

\[
8192\cdot12=98304
\]

block XORs to a full call. This cost is small compared with the inner and
outer kernels and should fuse with the outer permutation pass.

Before implementation, scaling the inner operation count suggested that the (s=19) state would raise
the measured reverse-inner time by about 0.42 ms. This projects an integrated
runtime near 11.2--11.3 ms from the 10.821 ms baseline. The estimate is not a
benchmark. The measured integrated result above is materially better.
