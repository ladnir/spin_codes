# Final theorem evidence coverage

Updated read-only snapshot: 2026-09-28, 23:43:24–23:43:35 UTC. The target is
`SpinCodes.Structured.ConcreteNativeTheoremPin`. No proof dependencies were
modified by the inventory audit. A separate isolated source recheck described
below has repaired 50 identified evidence gaps without changing
live objects: 13 concrete-outer modules, 23 older foundations, nine arithmetic
helpers, and five archival timestamp cases. Numerical certificates, final
theorem assembly, and all four distance-corollary modules have now passed.

**No current source/object pair disagreed with a recorded successful pair.**
There are archival coverage gaps: some earlier successful checks recorded
only a source hash, and other modules have only dependency/transfer snapshots.
Those records do not establish a cryptographic association
between the current source and the object produced by a particular compiler
run. File timestamps were not used as evidence of such an association.

The complete discovered graph contains 5,215 project-module nodes. All source
files now exist, and independent dry-run topological discovery succeeds.

| Strongest evidence available in the final inventory | Modules |
| --- | ---: |
| Current source/object pair matches an individual PASS compiler record | 4,320 |
| Current pair matches the initial Peach checked-dependency snapshot, without a stronger pair record | 34 |
| Prior PASS record has the current source hash, but no object hash | 795 |
| Final dependency inventory matches, but no individual compile-output pair is recorded | 66 |
| Pending numerical or final-assembly objects | 0 |
| Missing source files | 0 |

The Fourier controllers append completed entries only after both
module compilations return zero and four axiom audits pass; this was checked
in `scripts/route_fourier_box_batch.py` and the corresponding lane controllers.
Their completed entries were useful evidence while the batch status was
`RUNNING`; those Fourier batches now have whole-report `PASS` records.
The machine snapshot preserves the exact report path, JSON entry, and hashes.

The initial `scripts/map_data/peach_dependencies.json` snapshot matches all
691 of its imported module pairs. Most also have stronger later records;
34 rely on the initial snapshot as their strongest paired evidence. These
include `Prob`, `ConcreteMapData`, `ConcreteMaps`, `MapSpectrumBridge`,
`Transvection`, `Transfer`, and `VarianceBound`. Its documented role is an
archive of previously checked dependencies, not a fresh compilation of them.

The 795 source-only records include all 387 majorant numerical parts, the
refined majorant assembly, four numerical majorant helpers, and
403 sparse/polynomial modules. Their historical PASS reports and matching
source hashes were retained. All 795 also have matching transfer-cache
hashes for the current objects, but combining those facts does not turn the
old source-only compiler record into a recorded compiler-output pair.

The remaining 66 snapshot-only modules comprise 45 old cover modules,
19 majorant segment assemblies, and `Golay`/`GolayEnum`. Their hashes match
the final dependency inventory and mutable transfer cache. The final
assembly snapshot is explicitly classified as inventory evidence for those
imports; compiling the final theorem does not freshly compile its imports.
The following 13 newer semantic modules lacked a dedicated source/object pair
in the initial inventory. Their evidence gap is now repaired:

- `ConcreteOuterDefect`, `ConcreteOuterDefectEnvelope`, `ConcreteOuterDefectNumeric`;
- `ConcreteOuterEntropy`, `ConcreteOuterEntropyTransition`;
- `ConcreteOuterEnvelope`, `ConcreteOuterEnvelopeCoefficient`,
  `ConcreteOuterEnvelopeShiftEntropy`, `ConcreteOuterEnvelopeShiftSupport`,
  `ConcreteOuterEnvelopeShiftTransition`;
- `ConcreteOuterMajorantBridge`, `ConcreteOuterMajorantRegularity`,
  `ConcreteOuterMajorantSpectrum`.

`scripts/encoder-recheck-outer-isolated.py` freshly compiled all 13 current
sources on Windows in dependency order, with one Lean compiler and a
4,000 MB cap. Outputs went to a separate directory; reused dependencies were
read-only inputs. All 13 fresh `.olean` files were byte-for-byte identical to
their original live objects. A separate isolated pin audited one relevant
theorem per module and found only `propext`, `Classical.choice`, and
`Quot.sound`. All live source and object hashes were unchanged afterward.
The run took about nine minutes, with minimum observed free RAM 4.45 GiB.

The dedicated record is
`scripts/map_data/encoder_outer_isolated_verification.json` (`PASS`,
22:33:32–22:42:35 UTC). It retains each original-object hash, fresh-object
hash, source hash, exit code, command, log, and explicit byte-comparison
result. Only actual byte matches are exposed as a current paired-object
hash. The refreshed inventory recognizes all 13 repaired pairs. This is
fresh source evidence, not an inference from timestamps or transfer caches.

A second bounded run, `encoder_foundations_isolated_verification.json`,
freshly checked 12 current semantic foundations in dependency order:
`Framework`, `Distance`, `Selection`, `Cover.DenseTail`, and Structured
`Certificate`, `Regimes`, `Instantiation`, `Composition`, `Interleaver`,
`Enumerator`, `Schedule`, and `WeightedNorm`. All 12 source checks and 12
axiom audits passed, using one local compiler capped at 4,000 MB; minimum
observed free RAM was 3.91 GiB. All live files remained unchanged.

The first direct-compiler outputs differed from the cached objects because
of compiler metadata. They were not promoted on that basis. All 23 named
Framework declarations, including full types, proof bodies, and axiom lists,
printed identically under both import environments; the 348,253-byte
transcripts are recorded in `encoder_framework_comparison.json`.

Controlled recompilation on Peach resolved the differences and reproduced
**all 12 original objects byte for byte**. The original Lake setup supplies
package `spincodes`; preserving that field with an empty `importArts` map
reproduced ten objects exactly. `Instantiation` and `WeightedNorm` also embed
their original absolute source filenames. Passing the unchanged source
bytes through Lean's `--stdin` with those filename labels reproduced both
remaining objects exactly. No binary files were edited to obtain a match.

The successful records are
`encoder_foundations_package_peach_verification.json` (12 source checks and
12 audits, 23:10:30–23:12:40 UTC) and
`encoder_foundations_filename_peach_verification.json` (two controlled
rechecks and two audits, 23:14:30–23:15:00 UTC). Each records original/fresh
hashes, source hashes, setup hashes, commands, logs, and actual comparisons.
Only standard axioms occurred. They used one compiler with an 8,000 MB cap,
a 20 GiB host-memory floor, and isolated import/output trees. Live files
remained unchanged. The refreshed inventory recognizes all 12 exact pairs.

Earlier resource-limited attempts are retained: the local package test hit
the 3 GiB free-memory guard, and the first Peach test reached its explicit
4,000 MB compiler cap. The latter report and log remain under
`scripts/map_data/encoder_package_peach_20260928_160850/`. Those attempts
made no source/output pairing claim and did not change live objects.

One optional follow-up batch also completed before final assembly:
`encoder_semantic_package_peach_verification.json` freshly checked
`Accumulator`, `AccTuple`, `BASpectrum`, `BAGolay`, `BATails`, `UpperTail`,
`Majorant`, `BetaSub`, `EmptyEpoch`, `Coarse`, and `Simplex`. All 11 outputs
matched their original cached objects exactly, and all 11 axiom audits used
only standard axioms. It used the same isolated package/filename reproduction
method, one 8,000 MB compiler, and a 20 GiB host-memory floor; minimum observed
free memory was 54.43 GiB. This run changed no live proof files and did not
replay the historical Golay enumeration or numerical certificate trees.

The final optional helper batch,
`encoder_numeric_helpers_package_peach_verification.json`, rechecked all
nine remaining `Numeric` helper sources: `FixedDefs`, `BAEvalDefs`, `Fixed`,
`LogBounds`, `FixedLog`, `RatLog`, `BAConstants`, `BAExponent`, and `BAEval`.
These contain arithmetic definitions and analytic/checker soundness proofs,
not embedded box-certificate datasets. All nine reproduced their cached
objects exactly and all nine axiom audits passed (the pure definitions need
no axioms). The single compiler retained the 8,000 MB cap and 20 GiB host
reserve; minimum observed free RAM was 49.03 GiB. `Golay` was inspected but
excluded because it includes the full 4,096-word histogram enumeration.

The final timestamp preflight found 68 imports whose object timestamps precede
their source timestamps. Sixty-three `LowCancellationData` blocks already had
matching individual compiler records. Five had weaker archival coverage:
`Majorant.RefinedData`, `MapSpectrumData.Basis`, `SparseProgramDefs`,
`SparseMaxima`, and `SparseContributionSound`. All five were freshly checked
in the isolated tree using their original relative filename metadata and no
package label. Their outputs matched the original objects exactly and five
axiom audits passed. `encoder_mtime_gap_peach_verification.json` records the
checks; minimum remote free RAM was 62.24 GiB. The subsequent independent
`encoder_mtime_preflight_review.json` rehashed all 68 current source/object
pairs and confirmed matching individual PASS records for every one. The
initial review is retained separately. This repairs the evidence gap without
changing timestamps, source files, or live objects.

One historical **source-only** record is superseded: the source hash for
`MapSpectrumData.Block257` in `kernel_blocks_257_257.json` differs from the
current source. That module has a current matching PASS source/object pair,
so it is not an unresolved current mismatch. Two apparent absent Fourier
objects in the initial scan were completed during the scan; the repeated
check resolved them, and the final snapshot has zero nonmatching historical
source/object pairs.

The final graph has no pending or missing source/object files, no current
pair disagreeing with an individual successful compiler record, and no
unresolved report labels. `native_theorem_verification.json` is `PASS` for
the freshly compiled final theorem and pin. Its complete 5,215-module hash
inventory is retained, with imported snapshots distinguished from individual
compiler-output pairs as above.

Four additional corollary modules, outside this base-pin closure, were freshly
compiled after base completion. `encoder_native_distance_consequences_verification.json`
records their source/object pairs and seven actual axiom audits; all use only
`propext`, `Classical.choice`, and `Quot.sound`. The independent read-only
`encoder_native_distance_consequences_final_review.json` rehashed all four
current pairs, inspected both pin logs (five plus two audits), and confirmed
that the base verification record stayed unchanged. The final theorems are
`relative_distance_success_tendsto` and
`eventually_exists_rate_half_distance_gt_eleven_percent`, in
`Spin.Structured.ConcreteNativeFamily`. They state actual success probability
tending to one and eventual positive-probability seed realizations with
exact rate one half and relative minimum distance strictly greater than 11%.

Reproduce this review with:

```powershell
Set-Location C:/Users/peter/repo/permute_conv-github-bch/research/lean
& C:/Python314/python.exe scripts/encoder-audit-native-evidence.py
```

The machine-readable result is
`scripts/map_data/encoder_native_evidence_coverage.json`. It records the
review time window, input-report hashes, every module's current hashes,
matching evidence classes, nonmatching historical records, and missing
sources. It is an evidence inventory, not a compiler verification record.

Recommended completion wording: report the final kernel-checked theorem and
its axiom audit, while stating that previously checked imported objects were
reused. Preserve the older successful-check records and these coverage
distinctions. A uniform claim that every current project object was freshly
recompiled from its recorded source requires the optional source replay in
`FINAL_REPRODUCTION.md`; taking a new hash snapshot alone does not establish
that stronger claim.
