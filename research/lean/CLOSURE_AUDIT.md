# Concrete theorem closure audit

Updated 2026-09-28. **CLOSED: the actual native distance/rate theorem and
unconditional corollaries pass.** Final default invariants report ALL CHECKS
PASS. The consolidated audit passes all 17 gates and checks 5,220 recorded
source/object pairs, with no missing/invalid gates or issues.

## Exact completion requirement

For the actual native construction, the probability of minimum Hamming
distance at most `floor(0.11 * Nsched m)` tends to zero. The law uses one
shared Golay-BAA outer and the actual independent route permutations and IMT
transvections. Every realization has rate `1/2`. The paper, original statement
pin, construction, native schedule, and threshold are preserved. No holes,
new project axioms, or compiled-trust certificates are used.

The ordinary/transposed encoder operation bound and arbitrary requested-length
wrapper remain outside the agreed native distance/rate scope.

## Requirements and evidence

Reports below are in `scripts/map_data/`. Each PASS applies to its stated
scope. An axiom audit alone does not discharge an implication's hypotheses.

| Requirement | Concrete proof/evidence | Result |
| --- | --- | --- |
| Exact table maps and finite spectrum/count bounds | `ConcreteMaps`, `ConcreteCounts`; `concrete_counts_verification.json` | Checked |
| Actual permutation routing, serialization, IMT recurrence | `ConcreteRoutedEncoder.stream_moment_bound`; `routed_bridge_verification.json` | Checked |
| Native schedule, shared outer, actual emitted codeword | `concreteFamily`, `ConcreteNativeCodeword`; `native_codeword_verification.json` | Checked |
| Linearity, injectivity, dimension N/2, rate 1/2, pairwise minimum distance | `realizedCode_rate`, `concrete_probBad_minimumDistance`; `native_linear_code_verification.json` | Checked |
| Outer selection failure and growing sparse range | `native_selection_failure_tendsto`, `native_sparse_eventually`; `outer_tail_closure_verification.json` | Checked |
| Every fixed positive occupation | `ConcreteNativeFixedAll.native_fixed_EZ_tendsto`; `native_fixed_all_verification.json` | Checked without analytic premises |
| Closed rectangle cover, selected profile/layer counts, vanishing native prefactors | `DenseGeometryRate`, `DenseOccupationFamilyRateLayer`, `ConcreteNativeDenseRates`; `dense_native_integration_verification.json` | Checked; numerical argument now supplied |
| All 283 scalar indexed certificates and aggregate | `DenseOccupationScalarCertified.scalar_certified`; scalar geometry and mixed reports | PASS |
| All 433 occupation indexed certificates and aggregate | `dense_occupation_complete_verification.json` | PASS: 133 witnesses, 1,566 pairs, 2,132 audits |
| All 307 Fourier indexed certificates and aggregate | `dense_fourier_complete_verification.json` | PASS: 105 witnesses, 824 pairs, 1,543 audits |
| Complete mixed numerical rate | `DenseGeometry.certified_denseRates`; `dense_mixed_aggregate_verification.json` | PASS: eta=4/10000000, C=24000000000000 |
| Unconditional minimum-distance limit and rate pin | `ConcreteNativeTheorem`, `ConcreteNativeTheoremPin`; `native_theorem_verification.json` | PASS: two fresh modules, four audits |
| Success probability tends to one; eventual positive-probability good codes | `ConcreteNativeDistanceConsequencesFinal` and pin; `encoder_native_distance_consequences_verification.json` | PASS: four fresh modules, seven audits |

`ConcreteNativeFixedAll` combines actual Q=1, Q=2, and Q≥3 limits. Finite
subset insertion, gap exchangeability, order-statistic integration, Laplace
domination, matrix path expansion, all 62 numerical endpoints, and the native
profile application are supplied. The finite sum over occupations 1–4095
tends to zero without requiring uniformity over a growing occupation range.

`minimum_distance_of_fixed_and_dense_rates` is a reusable implication. The
final `native_minimum_distance_failure_tendsto` supplies both its fixed-moment
and mixed `DenseRates` arguments. Its full printed statement has no remaining
analytic or certificate premises.

## Independent reviews and trust checks

`NATIVE_SEMANTIC_REVIEW.md` checks the construction against the paper, including
all 38 map basis words, shared outer, uniform transvection law, persistent
state, and threshold. `DENSE_SEMANTIC_REVIEW.md` checks density parameters,
closed boundaries, inequality directions, profile counting, selection
averaging, and native remainders. `FIXED_NATIVE_SEMANTIC_REVIEW.md` checks the
actual first moments, continuum normalizations, and finite occupation sum.
Conservative intermediate constants differ from the paper but establish the
same limit. None of these reviews found a blocking mismatch.

`FINAL_STATEMENT_REVIEW.md` checks the actual printed final statements,
recursive axiom outputs, current source/object hashes, and report links.
Every final recursive audit uses only `propext`, `Classical.choice`, and
`Quot.sound`. `final_source_trust_scan.json` binds the current inventory of
5,299 project Lean sources and reports zero proof-hole, new-axiom, or
compiled-trust token hits after stripping nested comments.

Original `SpinCodes/Pin.lean` SHA-256:
`447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca`.
Paper `.tex` hashes match the earlier verified parallel checkpoint; this is
preservation of that working-tree baseline, not a claim of a clean Git tree.

The two occupation generator template errors found during reproduction
review were repaired before certification. B152 and B372 both passed their
individual kernel checks, and the generator audit checks all 433 associations.

## Evidence scope and final gates

The final assembly snapshot records 5,215 project imports. The independent
coverage inventory finds 4,320 current source/object pairs with individual
successful compiler records, 795 historical source-only PASS records, 34
initial checked-dependency snapshot records, and 66 final/transfer snapshot
records. There are no pending sources or objects and no mismatching paired
records. These categories must not be conflated with fresh compilation.

Fifty modules with weaker archival evidence were independently reproduced
in isolated outputs, byte-identical to their live objects. For 68 timestamp
inversions, `native_dependency_timestamp_verification.json` requires direct
current source/object compilation evidence: 63 existing low-cancellation
records and five fresh exact reproductions. No timestamps or live dependency
objects were modified to satisfy the check.

The final default check is recorded separately in
`native_final_invariants_verification.json`. The read-only consolidated check
is `final_closure_evidence.json`; it requires all mandatory reports and
revalidates current hashes. The default build does not itself import the
actual native theorem, so dedicated checks are required in addition to it.

A fresh clean replay of every transitive project source has **not** been run
in this closing phase. `FINAL_EVIDENCE_REVIEW.md` records exact coverage;
`FINAL_REPRODUCTION.md` documents the optional clean source replay. This is
the recommended next verification step after closure, not a remaining
mathematical hypothesis. No commits were made.
