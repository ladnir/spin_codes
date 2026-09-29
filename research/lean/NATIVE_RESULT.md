# Native Structured SPIN result

**Status: CLOSED for native distance and exact rate.** All numerical families, mixed certificate, actual theorem, full-statement pin, and unconditional success/existence corollaries pass. Final default invariants pass, and `scripts/map_data/final_closure_evidence.json` records all 17 consolidated gates passing with 5,220 current recorded source/object pairs and no issues. Previously built dependencies were reused; a clean full source replay was not run.

## Construction and claim

The result concerns the concrete binary linear code emitted by Structured SPIN at its native lengths. Lean indexes the family by `m : ℕ` and defines

```text
L = Lsched m = 128*(m+1)
b = bsched m = least positive multiple of 24 with b ≥ (39/4)*log₂(L*b)
N = Nsched m = L*b.
```

The `m+1` shifts the paper's positive index to Lean's natural-number index. These definitions are in `SpinCodes/Structured/Schedule.lean`.

Setup samples two independent uniform permutations for one Golay-BAA outer constituent and reuses that same constituent in all `L` rows. Independently, it samples the row permutations, region permutations, and IMT transvections. The theorem's probability is exactly the product of `ConcreteOuter.nativeSeedLaw m` and `ConcreteRoutedEncoder.experimentLaw L b (rounds m)`, with `rounds m = b*L/128`. It is the original setup distribution; the outer selection event used inside the proof is not a premise of the final claim.

For each setup, `realizedCode m out inner` is the actual emitted code, and `minimumDistance m out inner` is its minimum Hamming distance. The final theorem states

```text
Pr[minimumDistance ≤ floor(0.11*N)] → 0  as m → ∞.
```

Every realization has dimension `N/2` and exact rate `1/2`. The success corollary gives probability tending to one that `minimumDistance/N > 11/100`. The existence corollary gives such a rate-half realization, with positive setup probability, for every sufficiently large native index.

This matches the distance equation `eq:structured-spin-distance` and the exact-rate assertion in the paper's Theorem `thm:structured-spin-scalable` (`../paper/structured_proof.tex`). The declarations below do not establish the theorem's encoder/transposed-encoder operation counts or the wrapper for arbitrary requested lengths. They also provide no practical finite-length cutoff for the asymptotic probability statement.

## Lean entry points

```lean
import SpinCodes.Native
```

This public import exposes the original declarations without replacing their
definitions. `SpinCodes.NativePin` prints their full types and recursive axiom
dependencies. Start with [PROOF_GUIDE.md](PROOF_GUIDE.md) for the reading order.
Both public modules have passed their dedicated compilation check; see
`scripts/map_data/native_public_entry_verification.json`. The original
closure sources and certificate objects were preserved during this cleanup.

The main declarations are:

| Declaration | Meaning |
| --- | --- |
| `Spin.Structured.ConcreteNativeFamily.native_minimum_distance_failure_tendsto` | Actual failure probability tends to zero. |
| `Spin.Structured.ConcreteNativeFamily.realizedCode_rate` | Exact rate `1/2` for every setup. |
| `Spin.Structured.ConcreteNativeFamily.relative_distance_success_tendsto` | Strict relative distance above `11/100` with probability tending to one. |
| `Spin.Structured.ConcreteNativeFamily.eventually_exists_rate_half_distance_gt_eleven_percent` | Eventual existence of a positive-probability setup with that distance and rate. |

The original explicit statement pins remain in
`SpinCodes.Structured.ConcreteNativeTheoremPin` and
`SpinCodes.Structured.ConcreteNativeDistanceConsequencesFinalPin`. They are
checked by the dedicated closure controllers described below.

## Verification

The approved manuscript presentation update is recorded separately in
`scripts/map_data/paper_revision_verification.json`. It preserves the theorem
statement and Lean sources. Earlier closure reports bind the preceding paper
snapshot; their paper-preservation check will detect the intentional revision.

Run commands from `research/lean`, using the pinned Lean/mathlib environment. On this Windows workspace, the read-only evidence check is:

```powershell
C:/Python314/python.exe scripts/route-final-closure-evidence.py
```

It requires all mandatory gates to pass and current recorded hashes to match. Its output separates individual check records from dependencies covered only by an assembly snapshot. The default `scripts/check.sh` build does not itself import the concrete native theorem; dedicated checks are required.

To repeat the final assembly, while no other controller owns these checks, run:

```powershell
C:/Python314/python.exe scripts/check-native-theorem.py
C:/Python314/python.exe scripts/encoder-check-native-distance-consequences.py
C:/Python314/python.exe scripts/check-native-final-invariants.py
C:/Python314/python.exe scripts/route-final-closure-evidence.py
```

The first two commands use the configured Peach helper. They check the final declarations and axiom closures while reusing previously checked dependency objects. The invariant command runs the default checks locally. Accepted axiom audits contain only `propext`, `Classical.choice`, and `Quot.sound`; compiled numerical evaluation is excluded.

A full project-source replay is a separate, more expensive verification mode. Inspect its plan first:

```powershell
C:/Python314/python.exe scripts/replay-native-closure.py
```

Adding `--execute --backend peach --jobs 1` recompiles every transitive project source of the native theorem pin, after the replay's idle check. Lean and mathlib remain pinned external dependencies. An assembly PASS must not be described as this full replay; the replay produces its own report.
