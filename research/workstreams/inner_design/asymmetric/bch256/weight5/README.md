# IMT with weight-five feedback: full BCH-256 certificate at K=2^20

The [IMT inner](../../../IMT.md) now has a full finite distance certificate for the
fixed BCH [256,128,d>=38] outer at K=2^20, t=128, and s=19. The bound certifies
relative distance greater than 0.10 with setup failure below 2^-40.
Its diagnostic margin is **50.0620882264 bits**.

This result uses the same balanced expansion A as the quarter-rate candidate,
but replaces weight-three feedback with the fixed `weight5_seed0` map.
It does not change the supported encoder. The [isolated optimized transpose](implementation/README.md)
is now bound to these exact maps and the tested outer. It measures 10.155 ms
on Peach, versus 11.169 ms for the supported build in the same series.
The asymptotic theorem and complete optimized forward encoder remain open.

## Construction and probability space

Each 128-bit input block X emits X+Aq, then updates q to Mq+BX. State starts
at zero, persists across regions, and has no final flush. For each epoch,
the shared setup independently samples u uniformly from the nonzero 19-bit
vectors and v uniformly from u's orthogonal complement, including zero.
The mixer is M=I+uv^T. The row and region permutations retain the independent
uniform setup distributions of the BCH-256 model in the parent directory.

The maps A and B are fixed, not sampled in the distance claim. `candidate.py`
reconstructs B from the retained search recipe; each receipt also contains
its exact ordered columns and their hash. B has 128 distinct weight-five
columns and rank 19. The image of B^T has minimum weight 25. The kernel of B
has no words of weights one through three and has 74 words of weight four.
The expansion image retains weights 48, 56, 64, 72, and 80.

The bad setup event is existence of a nonzero message whose output weight
is at most floor(N/10), where N=2K=2^21. The stated probability is over
setup, not a fresh random map for each message.

## Complete occupancy coverage

Let Q count nonzero outer rows. There are 8192 outer rows.

| Part of the union | Covered Q | Margin of that part (bits) |
|---|---|---:|
| Single nonzero row | 1 | 50.062272096 |
| Sparse range | 2 through 511 | 62.999999909 |
| Dense cover, 260 rectangles | 512 through 8192 | 78.909866698 |
| Full union | 1 through 8192 | 50.062088226 |

The producer uses 256-bit outward arithmetic. Every retained bound passed
512-bit replay; Q=1 replay uses linear epoch iteration as an alternate
evaluation order. The sparse range uses the existing directed positive-fold
routine and candidate-specific fixed-weight transfers. The dense cover uses
candidate-specific complete moment bounds within the existing outer-only
two-tilt reduction. It never combines entries from different transfer bounds.

The earlier weight-three tree supplies only subdivision geometry and witness
proposals. All bounds were reevaluated for B. After 120 refinements, all 260
leaves passed. The verifier reconstructs the partition and sums the retained
dyadic bounds exactly. The full union is dominated by Q=1.

The outer bounds remain rigorous spectrum envelopes, not a guessed exact
BCH-256 spectrum. Verification authenticates the existing outer receipts
and reconstructs the shell caps as in the previous migration verifier.
It does not rerun the original LP solvers or provide a separate formalization
of their proofs.

## Smaller sizes and design implications

The same weight-five map does not yet have full certificates at K=2^16 or
K=2^18. Replayed sparse coverage reaches Q=92 and Q=255, respectively.
Their Q=1 margins are 43.592833502 and 50.189181114 bits.

At normalized density coordinate v=21/128, the retained dense witness gives
about -1559 bits for K=2^16,Q=218 and -1246 bits for K=2^18,Q=870.
Balanced feedback B=A^T also fails these point checks, at about -1549 and
-1207 bits. These are failed upper bounds, not counterexamples to the codes.
Increasing feedback density alone is therefore not an established solution
for the shorter instances. Tighter finite-length bounds remain worth testing.

This is a rate-specific outcome: weight-three feedback already has the
quarter-rate certificate, while weight-five feedback now closes the half-rate
K=2^20 instance. Neither result calls for replacing every configuration at once.

## Files and checks

- `candidate.py`: reconstructs both spectra, kernel, fiber caps, and low-input
  cancellation histograms; provides an isolated engine and dense checker.
- `screen.py`: bounded Q=1, sparse, and difficult-point checks and replay.
- `ranges.py`: contiguous sparse coverage and replay.
- `cover.py`: bounded dense partition search and replay.
- `verify_full.py`: receipt authentication, exact coverage, and exact union.
- `test_candidate.py`, `test_union.py`: nine passing checks, including exact
  small-weight fibers, independent cancellation histograms, alternate Q=1,
  and rejection of gaps, weak leaves, and altered unions.

The completed local records are `Q1_v2.json`, `RANGE_M20.json`, and
`COVER_M20.json`, with their `_replay.json` files. The aggregate record is
`FULL_M20_VERIFIED.json`. `Q1.json` is superseded: it predates the normalized
JSON key representation and is not an input to the full verifier.
The `RANGE_M16`, `RANGE_M18`, `DENSE_POINTS`, `SPARSE_POINTS`, and
`BALANCED_CONTROL` records and replays document the bounded checks.

Data files remain local and uncommitted, following the repository owner's
instruction. This code-and-documentation checkpoint is not a self-contained
certificate bundle. The scripts need the local map and outer-bound inputs.

From the repository root, with those inputs available:

```text
python workstreams/inner_design/asymmetric/bch256/weight5/verify_full.py
python -m unittest discover -s workstreams/inner_design/asymmetric/bch256/weight5 -p "test_*.py" -v
```

To produce new records, select unused output paths. `screen.py --mode q1`
produces the Q=1 bounds; `ranges.py --m 20 --last 511` produces the sparse
range. `cover.py --m 20 --minimum 512 --nodes 200 --seconds 180` accepts
`--seed workstreams/inner_design/asymmetric/bch256/DENSE_ADAPTIVE_M20.json`
for the old geometry. Every producer requires `--output`; replay uses the
same path with `--verify`. Replay paths are also write-once. The full verifier
uses the fixed retained filenames listed above. Do not run Python with `-O`.

## Recommended next step

The isolated optimized transpose, implementation binding, and serial comparison
are complete. Next integrate the certified rate-specific profiles behind an
explicit opt-in. Keep the supported default while the remaining migration
gates are open. The smaller sizes need tighter finite-length bounds; the
present evidence does not establish that denser feedback would close them.
