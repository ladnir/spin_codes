# Verification and certificate scripts

Run scripts from `research/lean`. Start with [FINAL_REPRODUCTION.md](../FINAL_REPRODUCTION.md) for the verification sequence and platform prerequisites. This directory includes current entry points, generators, and historical controllers; not every script is a suitable restart command.

## Current entry points

| Script and accepted arguments | Action and output |
| --- | --- |
| `route-final-closure-evidence.py` | Starts no compiler. Revalidates reports/current hashes; writes `map_data/final_closure_evidence.json`. Inspect its status, not just exit code. |
| `check-native-theorem.py [--wait]` | Peach compilation of actual native theorem/full pin after matching upstream PASS; writes `native_theorem_verification.json`. |
| `encoder-check-native-distance-consequences.py [--wait]` | Compiles four helper/final corollary modules after base PASS; writes `encoder_native_distance_consequences_verification.json`. |
| `check-native-final-invariants.py` | Windows Git Bash default checks after base/corollary PASS, followed by artifact preservation checks. |
| `check-native-dependency-provenance.py` | Revalidates successful compilation pairs for 68 timestamp-inverted imports; no compilation. |
| `replay-native-closure.py` | Dry-run import plan; `--execute` recompiles the selected project closure. Options: `--target MODULE`, `--backend peach` or `local`, `--jobs 1..4`, `--max-modules N`. |
| `peach-lean.py SOURCE.lean` | Compiles one source on prepared Peach using available project imports; fetches objects and writes `peach_<stem>.log`. |
| `check.sh` | Default legacy build, source trust scan, and legacy headline/certificate audits. Does not replace the actual-native checks. |

Use `--target SpinCodes.NativePin` for a replay of the public native entry and all four printed results. The replay script's default target is the older `ConcreteNativeTheoremPin`, which has a narrower import closure.

These tools write reports or objects even when their proof inspection is read-only. Preserve completed records before rerunning them. A wait option does not start its upstream producer. A `WAITING`/`RUNNING` record is not evidence of a live process. A transfer-cache entry is not a compiler record. [FINAL_EVIDENCE_REVIEW.md](../FINAL_EVIDENCE_REVIEW.md) explains the recorded evidence classes.

## Numerical certificate families

Generators emit candidate data or Lean source; checkers compile that source. Many batch scripts do both. These maintenance tools are not required to import the already checked theorem.

| Family | Main files | Scope |
| --- | --- | --- |
| Scalar dense | `dense_scalar_exact.py`, `dense_scalar_exact_batch.py`, `dense_scalar_geometry.py`, `dense_scalar_geometry_batch.py` | 283 boxes and indexed wrappers. Single generators take an index; the exact batch accepts `--start`, `--stop`, `--workers`. |
| Fourier dense | `route_fourier_exact.py`, `route_fourier_box.py`, `route_fourier_batch.py`, `route_fourier_box*_batch.py` | 105 witnesses, 307 boxes. Single generators take an index. The lane controller takes a configured name from `map_data/dense_fourier_box_lanes.json`. |
| Occupation dense | `dense_occupation_fixed_*.py`, `dense_occupation_admitted_resume.py` | 133 witnesses, 433 local boxes, 433 indexed wrappers. The admitted controller takes `matrices`, `boxes`, or `rates` and coordinates compilation permits/hash-gated parents. |
| Geometry/aggregation | `dense_occupation_geometry.py`, `dense_occupation_global_indices.py`, `dense_occupation_aggregate_wait.py`, `dense_occupation_mixed_aggregate_resume.py` | Global indexing and combined 1,023-box certificate. |
| Fourier evidence | `route_fourier_index_audit.py`, `route_fourier_final_verification.py [--wait]` | Index/source and recorded proof verification; no Lean compilation. |
| Occupation evidence | `dense_occupation_complete_wait.py` | Waits for six granular PASS reports, then checks exact module/audit sets and current hashes; no compilation. |
| Fixed-occupation numerics | `route_fixed_numeric_generate.py`, `route_fixed_numeric_semantics.py` | Endpoint data and its semantic connection. |
| Outer envelopes/tails | `route_defect_*.py`, `route_dense_*_generate.py` | Scalar defect and outer-tail certificate construction. |

Dense generators consume candidate replay data under `../workstreams/inner_design/imt_asymptotic/d11/`, including `DENSE_REPLAY.json`. Their external or floating-point candidate calculations are not trusted proof steps: the emitted inequalities are checked in Lean.

Do not restart completed lanes merely to inspect status. They may regenerate sources, overwrite manifests, or share objects. `dense_occupation_parallel_resume.py` is explicitly disabled: it exits before imports or I/O, retaining its old implementation only for historical inspection. Its weaker existence/old-log admission was superseded by the admitted controller. `dense_occupation_mixed_aggregate_wait.py` is also an older entry point superseded by the mixed resume controller; older waiters can encode an earlier lane layout. Review source and current manifest ownership before deliberate regeneration; [FINAL_PIPELINE_REVIEW.md](../FINAL_PIPELINE_REVIEW.md) records the handoff repairs.

## Earlier layers and development tools

| Files | Purpose |
| --- | --- |
| `emit_cover.py`, `gen_dense_tail.py`, `emit_dense_tail.py`, `check-cover.sh` (argument `small`, `full`, or `all`) | Earlier interval covers. The shell driver can skip an absent full-cover source; its success line alone is not native closure evidence. |
| `majorant*.py`, `emit_majorant.py`, `assemble_majorant.py`, `check-majorant.py`, `check-majorant.sh` | BA majorant search, generation, part checks, assembly. Some drivers default to several compilers. |
| `sparse_*.py`, `check-sparse-bridge.py`, `check-sparse-semantic.py`, `check-sparse.sh` | Sparse programs and semantic/certificate checks. The shell driver regenerates data before building. |
| `map_spectrum_*.py`, `map_fiber_*.py`, `check-map-*.py`, `concrete_*.py`, `check-concrete-*.py` | Exact map spectra, fibers, counts, and concrete transitions. |
| `emit-low-cancellation.py`, `check-low-cancellation.py`, `check-low-final.py` | Low-cancellation data and assembly. |
| `check-*-record.py`, `route_completed_verification.py`, checkpoint scripts | Narrow historical collectors. Inspect each exact module list/provenance rule; they do not replace final consolidation. |
| `fixmirror.py`, `check_mirror.py`, `majorant_why.py`, `majorant_progress.py`, preview/probe scripts | Development diagnostics; some modify generated files. |

## Evidence repair and semantic review

`encoder-audit-native-evidence.py` inventories existing source/object associations; `encoder_native_semantic_review.py` records a bounded semantic comparison; `final-source-trust-scan.py` performs a comment-aware prohibited-construct scan. These create reports without recompiling the native theorem.

`encoder-recheck-outer-isolated.py` and `encoder-recheck-foundations-isolated.py` run bounded local checks while preserving live objects. `encoder-recheck-package-peach.py` performs selected isolated remote checks and accepts `--original-filenames`, `--modules`, `--plan`, `--optional`, `--report-name`, `--no-package`, and `--relative-filenames`. A plan maps selected module names to audit theorem names. These options do not constitute a general proof of arbitrary cached artifacts. Its `encoder-peach-package-worker.py` helper is not a standalone environment installer.

`lean-local.ps1` constructs `LEAN_PATH` from installed package/project objects and forwards arguments to `lean`. `prepare-peach-counts.py` creates a historical transfer archive; `peach_dependencies.json` and `peach_extra_files.json` support transfer reuse, not proof authority. `map_data/` also contains compiler logs, producer manifests, candidate data, isolated-run directories, and consolidated reports.

## Execution constraints

Remote helpers contain machine-specific PuTTY paths and prepared Peach locations. The final invariant controller contains a Windows Git Bash path. The closure replay's local backend is the alternative when these are unavailable; pinned installed external dependencies are still required.

Use one controller owner per manifest. Never rebuild a dependency while another compiler consumes it. Keep memory headroom: the one-module helper permits 20,000 MB per compiler, and older drivers may default to multiple workers. Do not run benchmarks concurrently. Python reports describe orchestration evidence; Lean declarations and recursive axiom audits establish the proof.
