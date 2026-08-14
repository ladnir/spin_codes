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
| Riffle with packet width `g=4` | Certificate needs construction binding | The outward ledger gives `2^-61.7881515553` for the global-lane puncture variant. It does not certify independent puncture lanes in the sloped layout. |
| Riffle with packet width `g=8` | In progress; same binding gap | This is the intended deployment width. The conditioned-row branch needs the global-lane puncture rule or a proof-only replacement tax. |

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

## `g=4` certificate and layout-binding condition

Let the `g=4` setup use the global-lane puncture rule described below, and
sample the remaining lane bijections, packet permutation, and recursive state
permutations specified by Riffle. For a
sampled setup, let `Z_d` count nonzero messages whose codeword has weight at
most `d=188743`. The complete outward certificate proves

```text
E[Z_d] <= 2^-61.78815155534398186100.
```

The expectation is over the declared setup randomness. Markov's inequality
therefore gives

```text
Pr[d_min <= 188743] <= 2^-61.78815155534398186100.
```

The certified margin beyond the required 40-bit threshold is
`21.7881515553` bits. In particular, at least one frozen setup defines a
binary `[2^21,2^20,d_min>=188744]` code.

The verifier establishes the bound through the following finite ledger.

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

4. **Per-cell BSP.** An exact rational binary space partition refines 899
   source cells. Each tree covers its parent exactly. Closed boundaries may be
   overcounted because every parent retains its original profile-count bound.

5. **Fixed witnesses and mixtures.** A leaf uses one fixed witness or one
   fixed rational convex mixture at every vertex. Convexity then covers the
   complete leaf. Different witnesses at different vertices do not suffice.

6. **Finite-atlas pricing.** Discovery uses dual column generation. The final
   verifier does not trust its floating-point objectives. It independently
   hardens every selected witness and evaluates every fixed rational mixture.

7. **Cell-local union bound.** Each certified cell bound is multiplied by a
   rigorous integer-profile count upper bound. The verifier performs an
   outward log-sum-exp over cells and supports.

Relevant mathematical notes are
`explorations/g2_triangular_cover_math.md`,
`explorations/g4_delaunay_cover_math.md`, and
`explorations/riffle_group_chain_proof.md`.

The independent replay reports:

```text
feasible exact-support strata                 30
exact support-mesh root cells              40902
root cells replaced by exact BSPs            899
exact BSP leaves                             2132
used fixed witnesses                         2144
component evaluations at vertices          483058
BSP leaf inequalities                        15284
outward log2 E[Z_d] upper endpoint   -61.78815155534398186100
```

Canonical artifacts:

| Artifact | SHA-256 |
| --- | --- |
| `out/g4_support_lower30_outward_sparse_grid9_mixtures.json` | `ddd57ef25a1ef815614ca3592941075a1ae12dc8652b8c76806c233c511a92d7` |
| `out/g4_support_full_asymmetric_incremental_peach.json` | `80b03f60b4b3c6c7b22d1c20cf2e1e203561cd9f111a6baf8bd84cf24e349ef8` |
| `out/g4_lower296_bsp_diagnostic.json` | `e1208a95642ea6092dbde94baea8015a6cf675b652294d6ae492e7a0a7442880` |
| `out/g4_cell_bsp_remaining603_asymmetric.json` | `4a67362e2a0b191617c93933e0d0474769adf6ffc4442bb002544c5b673c6e19` |
| `out/g4_end_to_end_bsp_outward_certificate.json` | `9a416987feda6c0edf15ea2af91b30eb84d3be8b97a002831cf85f03e7086364` |

The final report binds the evaluated inequalities with SHA-256
`512f46a7e6cf2686e7bcf88b8b8af22af11080d6741f8f60eaaba0bb785d4c6f`.
It binds the hardened witness reports with SHA-256
`adfc9e21df420c4c6dd0fc9823c0b2c457bc91a84d1e833a153fd5e3fa39f0e1`.

An independent adversarial audit re-derived the one-conditioned-row lemma,
cross-implemented its numerical branch, and reproduced the complete warm
replay bit for bit. The audit found no weakened verifier invariant. Its PASS
verdict is recorded by commits `762a9977fbd9526dad2d4845cfd3ebf3e5489098`
and `eb942354ec25897c3143b86e13d99a6768b4e7c7` on branch
`codex/g4-race-claude`. The second commit records that the optional fresh-cache
replay was stopped at the user's request. It was not needed for the verdict.

A later global-layout audit found a missing construction obligation. The
conditioned-row formula needs all 128 punctured blocks to lie in one perfect
matching of the three band-tile classes. Independent per-tile puncture lanes
do not ensure this. The numerical ledger therefore applies to the variant
that samples one global lane and punctures that lane in all selected tiles.
Uniformly sampling the global lane preserves each block's puncture marginal
and changes no hot-path operation. The exact argument and an allowed
counterexample to the independent-lane claim are in
`explorations/conditioned_row_sloped_layout_audit.md`.

The theorem uses two declared inputs beyond the finite certificate machinery.
The EBCH and graph spectrum tables are authenticated mathematical inputs. The
hole analysis also assumes the construction's graph-syndrome uniformity. A
proof of either input would strengthen the dependency chain but would not
change the finite ledger.

`G4_CERTIFICATE_MANIFEST.json` freezes the result, hashes, assumptions, and
source commits. Run the fast integrity gate before starting from this
checkpoint:

```powershell
python scripts\verify_g4_savepoint.py
```

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

The complete `g=4` gate is:

```powershell
python scripts\certify_packet_group_g4_end_to_end_bsp.py `
  --lower-ledger out\g4_support_lower30_outward_sparse_grid9_mixtures.json `
  --full-ledger out\g4_support_full_asymmetric_incremental_peach.json `
  --bsp-batch out\g4_lower296_bsp_diagnostic.json `
  --bsp-batch out\g4_cell_bsp_remaining603_asymmetric.json `
  --workers 8 --iterations 40 `
  --checkpoint-dir out\g4_race_hardened_cache `
  --output out\g4_end_to_end_bsp_outward_certificate.json
```

The command must reproduce `passed: true` and an upper endpoint at most
`-61.78815155534398186100`. A diagnostic discovery run is not a substitute
for this independent replay.

Never run two benchmarks or long verification jobs simultaneously.
