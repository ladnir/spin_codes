# BCH-256 / RM2Sub: paper integration handoff

Milestone: 2026-09-07. The selected (t,s)=(128,19) construction combines
full finite distance/setup certificates with measured no-fanout implementation
performance. No additional parameter search is needed before drafting.

## Construction and claim

The outer is the fixed binary [256,128,d>=38] intermediate BCH constituent
specified in `T64_S20_FULL_CLOSURE.md` and instantiated by the bare encoder's
generator. It is not an arbitrary code with these dimensions. Its exact
spectrum remains unknown; the proof uses certified spectrum inequalities.

For K message bits, there are L=K/128 outer rows and N=2K output bits.
Setup uses independent row permutations, independent transposed-region
permutations, and fresh independent nonzero GF(2^19) multipliers alpha_i.
The selected fixed maps A and B satisfy B=A^T and BA=0. The inner recurrence is

    q_0 = 0
    Y_i = X_i + A q_i
    q_(i+1) = alpha_i q_i + B X_i

The state persists across regions; output precedes update; there is no flush
and no parity fanout. One setup theta is shared by every message. The ledgers
pin the selected maps, field representation, parameters, and instance identities.

For each listed K, let H=floor(N/10). The exact rational bound U in its
ledger satisfies

    Pr_theta[exists x != 0: wt(E_theta(x)) <= H] <= U < 2^-40.

The reported margin is -log2(U), rounded for display. This is a lower bound
on the distance/setup margin, not an estimate of the true failure probability
or the total security level of the application. No independence between
codewords or heuristic spectrum assumption enters this statement.

## Results to use together

| log2 K | Full certified margin (bits) | Measured online time (ms) |
|---:|---:|---:|
| 16 | 53.9443672720 | 0.560 |
| 18 | 52.3463883689 | 2.322 |
| 20 | 50.4482033129 | 11.259 |
| 22 | 48.4706838718 | Not measured |
| 24 | 46.4762221823 | Not measured |

Each proof includes every occupancy Q=1,...,L, where Q counts nonzero outer
rows. These are five finite results, not uniform coverage of every K between
them. The implementation accepts message exponents 16 through 20 only.

Timings concern the transposed block encoder: 2K input blocks map to K output
blocks, with each 128-bit block carrying 128 parallel binary instances.
They are steady-state online times, excluding setup and allocation, on one
Ryzen 9 7950X thread pinned to CPU 15 under GCC 15.2.0. Each entry is the
median of three serial run medians, with 101 calls per run after warmup.
They are not full-application timings. At K=2^20 the selected (64,20)
reference takes 12.059 ms; (128,19) takes 6.6% less time.
`../bare_bch_rm2sub/PERFORMANCE.md` records setup, memory, flags, and variability.

## Proof structure and evidence map

1. Bound the bad-message expectation by occupancy, then sum all contributions.
2. Use the BCH weighted inequalities and positive inner transfer for Q1.
3. Cover the sparse remainder with the positive polynomial recurrence.
4. Cover dense occupancies with four-state interval bounds, two tilts, and
   the sharpened shuffled Poisson-binomial density comparison. The all-one
   label is counted separately and included in the complete sum.
5. Reconstruct the exact rational union after checking coverage and receipts.

| Material | Entry point |
|---|---|
| Current construction and earlier reference | `CURRENT_UNDERSTANDING.md`, `T64_S20_FULL_CLOSURE.md` |
| K16 / K18 split proofs | `T128_S19_M16_CLOSURE.md`, `T128_S19_M18_CLOSURE.md` |
| K20 / K22 / K24 closure | `T128_S19_M20_CLOSURE.md`, `T128_S19_M22_M24_CLOSURE.md` |
| Scalable dense argument | `DENSE_TWO_TILT.md`, `POISSON_DENSITY_REFINEMENT.md` |
| All-rung audit and regression gates | `audit_ladder_completion.py`, `test_audit_ladder_completion.py` |
| Measured implementation and provenance | `../bare_bch_rm2sub/README.md`, `PERFORMANCE.json` in that directory |

The retained completion audit freshly checked all 1,497 dense boxes at 768
bits, authenticated sparse 512-bit replays, reconstructed the exact unions,
and checked preservation of earlier certificate pins. The completed ladder
regression run passed 71 tests. The audit is not a formal proof-assistant
verification; the mathematical inequalities remain part of the proof argument.

## Engineering interpretation

For this selected map, the certified m20/m22/m24 margins lose almost one bit
per doubling of K. Q1 dominates all three full unions. A local planning rule is

    M(K) approximately 50.45 - log2(K / 2^20).

Present the finite certified points as such and any interpolating line as a
guide. Smaller exact-spectrum BCH studies explain the large-state plateau
through late activation by low-weight outer words. Higher occupancies still
constrain how far the epoch size can grow or the state size can shrink.

The spectrum-model curve in `CURVE_SPECTRUM_ASSESSMENT.md` uses the selected
(64,20) map, not this (128,19) map. Its roughly 20-bit gap from the certified
curve is unresolved spectral uncertainty, not demonstrated extra margin.
Either label that comparison separately or recompute the matched (128,19)
model before overlaying it. Backtests on known spectra support the model but
do not provide an error guarantee for BCH-256 or a universal BCH growth law.

## Artifact preservation and reproduction

Git retains the source, tests, proof notes, compact full ledgers, and completion
audit. The audit's `source_sha256` records the dependency hashes. Bulk search
trees, numerical worker receipts, build products, and benchmark logs remain
local. A source-only checkout is therefore not a self-contained replay bundle.
Before external artifact submission, restore or regenerate the pinned inputs
and receipts, or publish a separate reviewed artifact archive. No archive was
uploaded as part of this cleanup.

The repository attributes preserve bytes for the packaged proof and benchmark
sources. Older calibration dependencies already use explicit CRLF checkout
attributes; the checker respects those existing rules. Use a Git checkout
when restoring sources, since a raw blob export need not apply those rules.

Run the lightweight provenance check from the repository root:

    python -B workstreams/bch_rm2sub_bridge/verify_paper_milestone.py

It checks retained evidence, exact final sums, measured-source hashes, and
available dependency hashes; it does not rerun interval arithmetic. Add
`--require-local-evidence` to reject missing local dependencies, and
`--check-index` before committing to check the checkout bytes of staged files.
Full numerical reproduction commands are in the size-specific closure notes.
Use fresh output names and preserve the hash-bound sources and proof notes.
Benchmark reproduction is separate and must run serially.

## Next paper pass

Start from the current GitHub manuscript, reconciling concurrent paper edits
before integration. The local `paper/implementation.tex` has an explicit
handoff comment but still retains its earlier BA-3 text until that pass.
Do not identify that factored routing construction with this no-fanout encoder.

Present the construction and probability space before its finite theorem.
Then explain the proof decomposition, give the combined results table, and
interpret the engineering curve. Keep diagnostic search history and complete
receipt listings in the artifact documentation. Match notation and connect the
distance/setup event to the application's actual theorem hypotheses during
integration; do not silently relabel it as end-to-end protocol security.
