# Structured SPIN formalization — status

Updated 2026-09-28. Project: `research/lean/`. Paper:
`research/paper/structured_proof.tex` and its appendices.

The approved paper presentation update is now applied. The public and anonymous
PDFs build and the edited pages pass visual review. The theorem and Lean proof
sources are unchanged. `scripts/map_data/paper_revision_verification.json`
records the new manuscript hashes; earlier closure reports retain their original
paper snapshot, so their preservation gate will detect this intentional revision.

## Native distance and rate theorem

**CLOSED: the unconditional native distance/rate theorem and its corollaries
are kernel-checked.** Final default invariants pass. The consolidated audit
passes all 17 gates, checking 5,220 recorded source/object pairs with no
invalid gates or issues. See `scripts/map_data/final_closure_evidence.json`.

For the actual native Structured SPIN construction,

```text
Pr[minimumDistance ≤ floor(0.11 * Nsched m)] → 0.
```

The probability uses the original product law: the shared Golay-BAA outer
seed and the actual independent row/region permutations and IMT
transvections. Every realization has exact rate `1/2`. The success
probability for strict relative distance greater than `11/100` tends to one;
for every sufficiently large native index, such a realization has positive
setup probability. There are no remaining selection, fixed-occupation,
dense-rate, or numerical-certificate premises in these final statements.

Lean entry points:

```lean
import SpinCodes.Native
```

`SpinCodes.NativePin` audits the four results through this public import.
`PROOF_GUIDE.md` gives the reading order; `CONTRIBUTING.md` and
`scripts/README.md` describe maintenance and verification. The public entry
modules supplement the original checked proof without changing its sources.
Both modules compile; `scripts/map_data/native_public_entry_verification.json`
records four standard-only audits, a 5,301-source trust scan, and a dry-run
replay plan for the public import. The original 17-gate closure audit still
passes after the cleanup.

`ConcreteNativeTheoremPin.actual_minimum_distance` pins the full probability
statement and exact floor threshold. `actual_rate` pins the rate for every
seed. `NATIVE_RESULT.md` describes the construction and all final declarations.

## Completed proof obligations

| Component | Result |
| --- | --- |
| Concrete table maps, spectra, finite counts, Fourier/fiber bounds | Checked |
| Actual route permutations, persistent IMT state, emitted words | Checked |
| Shared Golay-BAA outer, native schedule, linearity and injectivity | Checked |
| Outer entropy majorant, selection tail, growing sparse range | Checked |
| Every fixed positive occupation and the finite sum over 1–4095 | Checked |
| 283 scalar, 433 occupation, and 307 Fourier indexed rectangles | All 1,023 checked |
| Mixed dense rate, selected profile/layer sums, native remainders | Checked |
| Actual unconditional minimum-distance limit and rate pin | PASS |
| Success probability and eventual positive-probability existence | PASS |

The occupation certificate covers 133 matrix witnesses, 1,566 recorded
source/object pairs, and 2,132 exact standard-only axiom audits. The Fourier
certificate covers 105 witnesses, 824 pairs, and 1,543 such audits. The mixed
certificate proves `DenseRates (4/10000000) 24000000000000` and supplies the
last numerical argument to the actual theorem.

For fixed occupation, the actual finite subset marginal, routing law,
continuum order-statistic integral, matrix bound, and all 62 numerical
endpoints are supplied. The Q=1 and Q=2 first moments have eventual bounds
`1000 * exp(-b/50)` and `1000 * exp(-b/20)`. Each fixed Q≥3 has bound
`(1600/3) * exp(-Q*b/2000)`. Only a finite sum of fixed-Q limits is used;
growing sparse and dense occupations have their own estimates.

## Verification and scope

`scripts/map_data/native_theorem_verification.json` records the fresh final
theorem/pin checks, four recursive axiom audits, and a 5,215-module project
dependency snapshot. `encoder_native_distance_consequences_verification.json`
records four freshly checked modules and seven audits. All final theorem
axioms are among `propext`, `Classical.choice`, and `Quot.sound`.

`final_source_trust_scan.json` records a comment-aware scan of all 5,299
project Lean sources, with no proof holes, new axioms, or compiled-trust
constructs. The original `SpinCodes/Pin.lean` and paper source hashes are
preserved relative to the prior verified checkpoint. No commits were made.

The closing checks reuse previously built dependency objects. They are
**not a fresh clean replay of every transitive source**. `FINAL_EVIDENCE_REVIEW.md`
distinguishes direct successful source/object records from archival
source-only records and dependency snapshots. `FINAL_REPRODUCTION.md`
provides the separate full-replay workflow. The default `scripts/check.sh`
does not import the actual native theorem, so dedicated final checks remain
necessary.

The result covers native-length distance and exact rate. Encoder and
transposed-encoder operation counts, an arbitrary requested-length wrapper,
and a practical finite-length probability cutoff remain outside this result.
See `CLOSURE_AUDIT.md` for the requirement-by-requirement audit and the
independent native, fixed-regime, dense-regime, and final-statement reviews.

Earlier milestones and their evidence remain in `ROUTED_MOMENT.md`,
`CONCRETE_MAPS.md`, `CONCRETE_STEP.md`, `LOW_CANCELLATION.md`,
`SPARSE_BRIDGE.md`, `DECISIONS.md`, and the verification reports. Historical
references to open obligations in those checkpoint documents are superseded
by the final theorem and this status.

**Recommended next turn:** review the rendered formalization discussion, then
prepare a versioned proof artifact. The separate project-source replay in
`FINAL_REPRODUCTION.md` remains a verification follow-up. The distance/rate proof
needs no additional lemma.
