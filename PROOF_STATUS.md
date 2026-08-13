# Proof status

Updated: 2026-08-13

This file is the current proof ledger. It separates theorem-facing results,
complete computer-assisted certificates, and exploratory diagnostics. The
chronological research log through 2026-08-13 is preserved in
`PROOF_HISTORY_2026-08-13.md`.

## Status labels

- **Theorem-facing:** the manuscript states the result and identifies its
  mathematical assumptions.
- **Certified:** an outward-rounded verifier accepts a complete finite ledger.
- **Diagnostic:** binary64 computation or an incomplete regional ledger.
- **Candidate:** an unvalidated inequality, implementation, or proof route.

## Executive summary

| Line | Status | Current conclusion |
| --- | --- | --- |
| Random sliding dense construction | Theorem-facing | The asymptotic `0.109` checkpoint passes its sampled-grid verifier. |
| Structured RM/EBCH construction | Theorem-facing and checkable | A binary `[2^21,2^20,d_min>=188744]` code exists with a theorem-safe first moment at most `2^-37.27`. |
| Structured EBCH/XOR-parity construction | Theorem-facing and checkable | A binary `[2^21,1048512,d_min>=188744]` code exists with a theorem-safe first moment at most `2^-41.29`. |
| Riffle with packet width `g=2` | Certified | The optimized end-to-end first-moment bound is at most `2^-62.7194713852`. |
| Riffle with packet width `g=4` | Incomplete | No end-to-end theorem is claimed. The current leading diagnostic profile remains about `1409.7` bits short after convergence corrections and the best implemented exact-graph outer family. |
| Riffle with packet width `g=8` | Not started end to end | This is the intended deployment width. The proof methodology must first close at `g=4`. |

## Frozen Riffle target

The active construction uses the binary systematic extended BCH code
`[128,64,22]` for both local-code layers. The packet width `g` changes only
the granularity of the permutation between those layers. It does not shrink
either BCH code.

For `g=4`, the fixed parameters are

```text
N = 2^21
K = 2^20
d = floor(0.09 N) = 188743
M = N/4 = 524288 packets
profile a = (a0,a1,a2,a3,a4)
sum_j a_j = M
sum_j j a_j >= 21
```

The complete feasible profile set has

```text
binom(M+4,4) - 717
```

members. Its base-two logarithm is approximately `71.4151`.

Let `Z_d` count nonzero messages whose encoded word has Hamming weight at
most `d`. The finish line is a complete outward certificate for

```text
E[Z_d] <= 2^-40.
```

A uniform allocation would require each profile contribution to be at most
`-111.41506501642425` bits. The final proof may instead use a sharper
cell-local union ledger.

The probability space includes the setup randomness declared by the
construction: independent lane bijections, the packet permutation, and the
recursive state permutations. The packet-profile transitivity lemma and its
implementation invariant are recorded in
`explorations/riffle_group_chain_proof.md`.

## Complete `g=2` certificate

The `g=2` proof covers every integer packet profile. It combines convex
triangles with exhaustively enumerated terminal profiles. The verifier checks
the exact mesh, support eligibility, exact rational mixtures, terminal
ownership, source hashes, and outward witness evaluations.

The original global-max accounting gives the outward interval

```text
[-40.09421546136248628735, -40.09413878241724607147].
```

The optimized cell-local accounting gives

```text
[-62.71954806419136023328, -62.71947138524163556114].
```

Thus the current `g=2` certificate has at least `62.7194713852` bits of
first-moment security. This is `22.7194713852` bits beyond the requested
40-bit threshold.

Canonical artifacts:

| Artifact | SHA-256 |
| --- | --- |
| `out/g2_triangle_hybrid_complete.json` | `d45cf8ffacdabf7d74669db60eea958b01a5ee30fc72b290b4cab8ca986e2333` |
| `out/g2_triangle_hybrid_outward_parallel.json` | `c03d188708ee5dbc65c586cd7d35eb5cb554ed12af47b89f77d4d032f1580875` |
| `out/g2_triangle_hybrid_outward_cell_sum.json` | `e7245efa0c7ff6bf6cf25189f69501bbc26805ee8e7a27e5e7178d18ebebfd1d` |

The aggregation ledger inside the final report has SHA-256
`43ecbb47c2b61e9cb20890053e879a6adf863e14923a16f4a188a9d88e2e4b40`.

## Sound `g=4` proof architecture

The following parts are established proof machinery. They do not constitute
the missing end-to-end certificate by themselves.

1. **Packet-profile orbit.** Conditional on `a`, the permuted binary word is
   uniform on an orbit of size

   ```text
   Q_g(a) = M! / product_j a_j! * product_j binom(g,j)^a_j.
   ```

   This statement uses the independent within-packet lane randomization. A
   packet permutation alone is insufficient.

2. **Exact-support partition.** The profile domain is partitioned by its
   positive support. Sparse witnesses with zero fugacities cannot cross a
   support boundary.

3. **Exact source geometry.** Each support stratum has an exact integer hull
   and an exact rational simplicial mesh. The verifier checks determinants,
   facet incidence, boundary facets, and total volume.

4. **Per-cell BSP.** An exact rational binary space partition may refine one
   source cell without conforming to adjacent source cells. Each BSP tree must
   cover its parent cell exactly. Closed shared boundaries may be overcounted.

5. **Fixed witnesses and mixtures.** A leaf uses one fixed witness or one
   fixed rational convex mixture at every vertex. Convexity then covers the
   complete leaf. Different witnesses at different vertices do not suffice.

6. **Finite-atlas pricing.** The minimax LP uses dual column generation. It
   scans every eligible atlas column before declaring convergence. The final
   verifier does not trust the discovery LP; it independently evaluates the
   selected rational mixture.

7. **Cell-local union bound.** Each certified cell bound is multiplied by a
   rigorous integer-profile count upper bound. The verifier performs an
   outward log-sum-exp over cells and supports.

Relevant mathematical notes are
`explorations/g2_triangular_cover_math.md`,
`explorations/g4_delaunay_cover_math.md`, and
`explorations/riffle_group_chain_proof.md`.

## Current `g=4` frontier

The best available source-cell discovery ledger processes the leading 1024
source cells. It selects the best of four exact BSP trees per cell, preserves
their partitions, and reattaches every leaf over an atlas of 5587 frozen
witnesses.

```text
processed source cells                 1024
largest unprocessed source-cell term   +12078.982162066526
largest updated source-cell term       +16514.50512918122
leading cell                           s1fr033931
leading rounded profile                [429359,24531,38395,28039,3964]
```

This is a diagnostic regional ledger. It is neither a complete global cover
nor an outward certificate.

Canonical local artifacts:

| Artifact | SHA-256 |
| --- | --- |
| `out/g4_cell_bsp_top1024_best_collatz_merged_round4_reattached.json` | `69ff23973c9d38a1a93ac0378b4b38e834b83b0b913004340bb8870ad027fd17` |
| `out/g4_cell_bsp_top1024_best_collatz_merged_round4_reattached_failures.json` | `de996e74ce8622926e20deacf916a897d5f27f7bf1f6cd112d9f2261f6891930` |
| `out/g4_leader_429359_joint_fullcollatz_multistart.json` | `f7951ff6e116de9d81411b9abd214642f2fed9d6490735ec7ddc098a8fc31675` |
| `out/g4_leader_429359_joint_fullcollatz_exact_graph.json` | `1d7b5271ca477784f9e40b4796604030f76a4777798086387187690be0b90f1e` |
| `out/g4_leader_429359_exact_graph_linear_optimized.json` | `b9141dad0e74e170822a9c3981b240ee4e1c91ad4663a5a216f21366c31d5ee2` |

The merged BSP ledger names Peach source artifacts under
`/tmp/permute-conv-g2/out/`. Their expected hashes are embedded in the
ledger. A final certificate must make every source available to the verifier
or package its authenticated content.

### Leading-profile decomposition

The best stored 320-step witness under the old adversarial-hole linear-BL
outer gives

```text
outer                         +268976.6146352452
inner                         -267070.1202656985
combined                        +1906.4943695467
uniform target                   -111.4150650164
deficit                           2017.9094345631
```

Longer power iteration on the same frozen parameters converges near

```text
inner                         -267186.9845
combined with old outer         +1789.6301
inner MGF                         30.8161
state-domination charge           25.4891
```

The exact graph-conditioned linear-BL outer reduces the outer contribution
to approximately `268485.2337363`. Reoptimizing that outer family improves
it by only about `0.005` bits. Combining the converged inner diagnostic with
that outer gives approximately

```text
combined                        +1298.2492
uniform target                   -111.4151
remaining pointwise deficit       1409.6643
```

The 20,000-step trajectory is packaged as
`out/g4_leader_429359_collatz_20000.json` (SHA-256
`7e0f634d169a84e7266670ec775b6b1459245c372fb84cf4c4aef20f92d1f76a`).
The file records the converged inner witness with the older adversarial-hole
outer. The companion
`out/g4_leader_429359_collatz_20000_exact_graph.json` (SHA-256
`330e262bde5d624ed937c3208c587c013eee158a8958ae687e5c4d5f3a2bc7c6`)
replaces that outer by the inherited exact graph-conditioned outer. The
separately optimized exact-graph outer recovers about another `0.005` bits.

### Relaxations tested at the leading profile

| Test | Observation | Interpretation |
| --- | --- | --- |
| Power iteration from 320 to 20,000 steps | About `116.9` bits recovered | Finite Collatz convergence mattered, but the trajectory is now effectively saturated. |
| Exact graph-conditioned outer | About `491.4` bits recovered relative to the adversarial-hole outer | The old hole tax was substantial. |
| Direct optimization of exact-graph linear-BL tilt | About `0.005` additional bits | This outer family is locally optimized. |
| Shared systematic-column inner bound | About `0.035` bits beyond the row-specific bound | Global column sharing does not own the remaining deficit at the current tilt. |
| Exact drive-pattern transfer at the current tilt | No material gain after convergence correction | The point-cap relaxation is not the main owner at this tilt. |
| Arbitrary 65-state test-vector improvement | The complete converged MGF is only about `30.8` bits | This component has no plausible 1400-bit gain unless the operator or tilt changes. |

These are diagnostic comparisons. They localize slack but do not prove that
the construction has the target distance.

## Active proof routes

The current witness family has reached its redesign trigger at the leading
profile. More Powell depth, more Collatz iterations, or another arbitrary
state vector is not the preferred next move.

### Primary candidate: band-pair-spectrum outer

The repository contains exact split spectra for bands `(0,1)` and `(1,2)`:

```text
out/ebch85_band01_split_spectrum.csv
out/ebch86_band12_split_spectrum.csv
```

The candidate route combines their exact membership probabilities through a
rank-one domination and degree-4 Finner inequality. The intended result is a
profile-shaped outer bound that retains more three-band structure than the
current linear-BL family.

`scripts/probe_packet_group_g4_pair_spectrum_outer.py` is an unfinished
prototype. It has not passed its mass invariant, mathematical audit, or a
hard-profile run. No value produced by that script should be cited yet.

### Fallback: retuned exact drive-pattern inner

`scripts/probe_packet_group_pattern_exact_transfer.py` retains all 4845
multisets of four packet-drive patterns. At the inherited tilt it gives no
material improvement. A joint retuning could expose a better saddle point,
but the current Python implementation is too expensive for a broad search.
Use caching or a native data-oriented kernel before launching a campaign.

### Stop condition for the current construction

Challenge the construction only after both structured routes fail at the
leading profile under a properly audited joint optimization. A negative
binary64 experiment is not a counterexample. An explicit codeword of weight
at most `188743`, or a rigorous lower bound showing the required first moment
cannot close within the permitted ensemble, would be contrary evidence.

## Theorem-facing non-Riffle rows

The structured RM/EBCH certificate uses `4096` copies of
`RM(4,9) [512,256,32]` outside the full-split EBCH inner. Exact rational and
outward artifacts prove a complete first moment at most `2^-37.27`.

The promoted field-symbol row uses `16383` independent systematic EBCH
`[128,64,22]` blocks plus one componentwise-XOR parity block. Its dimension is
`2^20-64`, its outer minimum weight is at least 44, and exact arithmetic proves
a complete first moment at most `2^-41.29`. Scalar extension uses the same
binary generator over `GF(2^128)` and requires no field multiplication.

The BCH256 and BCH512 spectrum projections remain heuristic. The BCH512 code
construction is authenticated, but its modeled spectrum is not a theorem
input. Do not promote either projection without an exact spectrum or a
rigorous BCH-specific envelope.

## Verification gates

The basic paper and theorem-row checks are

```powershell
python scripts\verify_dense_claims.py --delta 0.109
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
python scripts\build_bch512_256_candidate.py
python scripts\audit_bch512_spectrum_obstruction.py
python scripts\compare_outer_modes_fullsplit.py --delta 0.09
pdflatex -interaction=nonstopmode main_permConv.tex
pdflatex -interaction=nonstopmode main_permConv.tex
```

For the Riffle line, a successful diagnostic run is not a verification gate.
The final `g=4` gate must invoke an independent outward verifier on a complete
ledger and must reproduce a final upper endpoint at most `-40`.

Never run two benchmarks or long verification jobs simultaneously.
