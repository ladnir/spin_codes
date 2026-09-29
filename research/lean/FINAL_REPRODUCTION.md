# Reproducing the native theorem

The native theorem proves relative minimum distance strictly above `11/100` with probability tending to one and exact rate `1/2` for every realization. The public import is `SpinCodes.Native`; `SpinCodes.NativePin` prints the four canonical distance, rate, success, and existence statements and their recursive axioms. See [NATIVE_RESULT.md](NATIVE_RESULT.md) for the declarations and [PROOF_GUIDE.md](PROOF_GUIDE.md) for the proof structure.

The recorded closure passed on 2026-09-28: 17 evidence gates, 5,220 current source/object pairs, and no reported mismatch. It reused previously checked dependency objects. A complete project-source replay has **not** been executed. [FINAL_EVIDENCE_REVIEW.md](FINAL_EVIDENCE_REVIEW.md) distinguishes individual compiler records, historical source-only records, and dependency snapshots. The public entry/pin is an additional import surface; the original closure records concern `ConcreteNativeTheoremPin` and `ConcreteNativeDistanceConsequencesFinalPin`.

## Environment and preserved inputs

Run commands from `research/lean`. The recorded Windows checkout is:

```powershell
Set-Location C:/Users/peter/repo/permute_conv-github-bch/research/lean
lean --version
lake --version
```

Keep `lean-toolchain`, `lakefile.toml`, and `lake-manifest.json` unchanged. The toolchain is `leanprover/lean4:v4.34.0`; the manifest pins mathlib at `5ed2965256430c3649e86755f9576b54eca72435`. Lean and Lake must be on `PATH`, including inside Git Bash for shell checks. Python commands below use the recorded installation `C:/Python314/python.exe`; substitute your compatible Python installation when working locally.

With elan installed, the dependency-installation commands are:

```powershell
elan toolchain install leanprover/lean4:v4.34.0
lake exe cache get
git -C .lake/packages/mathlib rev-parse HEAD
git diff --exit-code -- lean-toolchain lakefile.toml lake-manifest.json
```

These install external dependencies, not project certificates. The remote helper is specific to prepared Peach: it requires the configured PuTTY executable/key and pinned SSH host key, `/tmp/spin-lean-peach/project`, the toolchain under `/tmp/spin-lean-peach/toolchain/lean-4.34.0-linux/bin`, installed external dependencies there, and the local `scripts/map_data/peach_dependencies.json` inventory. It does not provision a new machine. The connect-to-peach skill describes the configured connection. `prepare-peach-counts.py` uploads a historical dependency archive; it is not a complete environment installer.

Preserve original reports and objects before an intentional replay. Controllers overwrite named reports, and compilation can replace live objects. Package or source-filename metadata can change object bytes even for identical source text. A new compiler run must not be relabelled as the old one.

## Inspect completed evidence without compiling

The manuscript was subsequently revised with user approval to present the Lean
result. `map_data/paper_revision_verification.json` records that change and the
unchanged theorem statement. The earlier closure reports retain the pre-revision
paper hashes. Running the collector below against the revised manuscript will
therefore report a paper-preservation mismatch; it does not indicate a changed
Lean theorem. Preserve both records rather than replacing historical hashes.

```powershell
& C:/Python314/python.exe scripts/route-final-closure-evidence.py
```

This starts no Lean process. It verifies required reports, current hashes, exact numerical-family counts, axiom records, original pin/paper preservation, and final invariants. Inspect `scripts/map_data/final_closure_evidence.json`: its `status` must be `PASS` with no missing or invalid gates. Process exit alone is insufficient; `INCOMPLETE` is a valid diagnostic result.

The main records, all under `scripts/map_data`, are:

| Record | Scope |
| --- | --- |
| `native_theorem_verification.json` | Two fresh base theorem/pin modules; four recursive audits; dependency snapshot |
| `encoder_native_distance_consequences_verification.json` | Four fresh helper/final corollary modules; seven audits |
| `native_final_invariants_verification.json` | Subsequent default checks and preservation of final artifacts |
| `dense_fourier_complete_verification.json` | 105 witnesses and 307 indexed boxes |
| `dense_occupation_complete_verification.json` | 133 witnesses, 433 local and 433 indexed boxes, plus aggregate |
| `dense_scalar_exact_boxes_verification.json`, `dense_scalar_geometry_verification.json` | 283 scalar boxes and their indexed statements |
| `dense_mixed_aggregate_verification.json` | All 1,023 indexed boxes combined |
| `native_dependency_timestamp_verification.json` | Paired compiler evidence for 68 imports with inverted timestamps |
| `encoder_native_evidence_coverage.json` | Historical evidence classification, not another compilation |
| `native_public_entry_verification.json` | Two new public entry/audit modules; four recursive audits; preserved original closure |
| `native_public_source_scan.json` | Separate source scan including the public modules; the original closure scan is preserved |

## Repeat assembly with checked dependencies

This requires existing checked objects and matching producer reports. On the prepared Windows/Peach setup, run sequentially and stop on a failed report:

```powershell
& C:/Python314/python.exe scripts/check-native-theorem.py
& C:/Python314/python.exe scripts/encoder-check-native-distance-consequences.py
& C:/Python314/python.exe scripts/check-native-final-invariants.py
& C:/Python314/python.exe scripts/route-final-closure-evidence.py
```

The first two controllers accept `--wait` when upstream checks are running. A `WAITING` file does not establish that its controller is alive. The base checker requires matching fixed/mixed aggregate and dependency-provenance PASS records. `check-native-dependency-provenance.py` can revalidate the existing 63 historical and five isolated compilation pairs; it does not create missing compiler evidence.

The final invariant controller invokes `C:/Program Files/Git/bin/bash.exe scripts/check.sh`. The default check builds the legacy `SpinCodes` root and audits its existing theorem and certificates. It does **not** replace the dedicated native checks or import `SpinCodes.Native`. On another platform, `bash scripts/check.sh` runs the default checks directly; the Python controller itself has the recorded Windows Bash path.

For a single module, `peach-lean.py` accepts one project source path:

```powershell
& C:/Python314/python.exe scripts/peach-lean.py SpinCodes/NativePin.lean
```

Its imports must already have compiled objects, including `SpinCodes.Native`. It compiles only the requested source with `lake env lean -j1 -M 20000` and fetches its objects. A one-module PASS does not replay dependencies or refresh every consolidated report.

## Replay the public import closure from project sources

The planner discovers source imports and detects missing files/cycles. It defaults to a dry run:

```powershell
& C:/Python314/python.exe scripts/replay-native-closure.py --target SpinCodes.NativePin
```

Without `--target`, it instead selects `SpinCodes.Structured.ConcreteNativeTheoremPin`, which excludes the public wrapper and final corollary wrappers. Each invocation writes a separate `scripts/map_data/encoder_native_replay_<time>_<pid>/report.json`; a dry run has status `PLAN` and starts no compiler.

When all other Lean runs are idle and sources are frozen, execute the selected closure:

```powershell
& C:/Python314/python.exe scripts/replay-native-closure.py --target SpinCodes.NativePin --execute --backend peach --jobs 1
```

Use `--backend local` instead when the pinned external dependencies are installed locally and `lake` is on `PATH`. Both backends recompile each selected project module in dependency order, replacing its object in the existing build tree. This is a project-source replay, not a hermetic clean checkout or a bootstrap of Lean/mathlib from source. The executor has not completed a full run; planning and bounded supporting checks have been exercised.

Supported controls are `--target MODULE`, `--jobs 1` through `--jobs 4`, and `--max-modules N`. A bounded run reports `PARTIAL` unless it covers the entire selected closure. Each compiler has one worker and a 20,000 MB cap; four may require about 80 GB plus overhead. Peach checks inherit the helper's 900-second timeout. A failed or timed-out replay is not evidence of mathematical falsity.

The executor checks for active Lean processes before starting and holds its own replay lock. It cannot prevent unrelated processes from starting later. Do not edit sources or run another build against the same object tree during replay. On failure it stops new submissions and drains running jobs. A complete `PASS` records compilation of the selected closure and a stable source/pin snapshot; inspect the final pin log for its exact declarations and recursive axioms.

Replay may invalidate hashes in historical producer/assembly reports. Preserve its separate logs and report; do not rewrite old evidence to make mismatches disappear. Re-establish matching producer evidence before expecting the archived assembly/consolidation controllers to pass again.

## Focused checks and generators

[scripts/README.md](scripts/README.md) groups the generators, family checkers, and evidence tools. Isolated rechecks repaired 50 archival evidence gaps without replacing live objects; their narrower scope is recorded in [FINAL_EVIDENCE_REVIEW.md](FINAL_EVIDENCE_REVIEW.md). A fresh compilation producing different bytes is not proof that an older cached object is the same output.

Python schedules work, generates candidate data, and checks records. Lean checks the proofs. Final recursive audits must contain only `propext`, `Classical.choice`, and `Quot.sound`; a generator result or Python PASS marker does not substitute for those audits.
