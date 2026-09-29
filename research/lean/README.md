# SpinCodes — Lean formalization of Structured SPIN

The native Structured SPIN distance and rate result is complete. Under the
actual setup distribution, the probability that minimum distance is at most
`floor (0.11 * N)` tends to zero along the native lengths. Every realization
has exact rate `1/2`.

## Use the result

```lean
import SpinCodes.Native

open Spin.Structured.ConcreteNativeFamily

#check native_minimum_distance_failure_tendsto
#check realizedCode_rate
#check relative_distance_success_tendsto
#check eventually_exists_rate_half_distance_gt_eleven_percent
```

`SpinCodes.Native` is the documented public import. It exposes the original
checked declarations without aliases or replacement definitions.
[SpinCodes/Native.lean](SpinCodes/Native.lean) explains their objects and scope.
[SpinCodes/NativePin.lean](SpinCodes/NativePin.lean) prints all four full
statements and their recursive axiom dependencies.

The law samples one shared Golay-BAA outer constituent, independently of the
actual routing permutations and IMT transvections. The final probability is
unconditioned. Relative distance strictly above `11/100` holds with probability
tending to one, and a positive-probability good setup exists at every
sufficiently large native index.

## Build and verification

Use the pinned toolchain and dependency manifest; run commands from this directory.
With elan and the pinned dependencies installed:

```sh
lake exe cache get
lake build SpinCodes.Native
lake env lean -j1 -M 20000 SpinCodes/NativePin.lean
```

The build command creates missing or stale project dependencies. On a fresh
checkout, this includes substantial numerical certificate checking. It is
not a quick smoke check. The final command imports the built public result
and prints its statements and axiom dependencies.

Bare `lake build` still selects the historical `SpinCodes` default target.
That target and `scripts/check.sh` check the earlier abstract framework;
they do not by themselves establish the concrete native theorem.

The completed closure checks report only `propext`, `Classical.choice`, and
`Quot.sound` in the final theorem dependencies. Numerical certificates are
checked by Lean's kernel. The closing run reused previously built objects;
it was not a clean replay of every transitive project source.
[FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md) gives the exact workflows,
resource limits, platform prerequisites, and evidence distinctions.

## Read and maintain the proof

| Document | Purpose |
| --- | --- |
| [NATIVE_RESULT.md](NATIVE_RESULT.md) | Exact construction, theorem declarations, and result scope |
| [PROOF_GUIDE.md](PROOF_GUIDE.md) | Reading order and dependencies from the encoder to the final theorem |
| [STATUS.md](STATUS.md) | Current completion status and recommended follow-up |
| [CLOSURE_AUDIT.md](CLOSURE_AUDIT.md) | Requirement-by-requirement closure evidence |
| [FINAL_EVIDENCE_REVIEW.md](FINAL_EVIDENCE_REVIEW.md) | Source/object provenance and archival coverage limits |
| [FINAL_REPRODUCTION.md](FINAL_REPRODUCTION.md) | Reproduction commands and full source replay |
| [CONTRIBUTING.md](CONTRIBUTING.md) | Editing proofs, generated certificates, and verification records |
| [scripts/README.md](scripts/README.md) | Verification drivers, generators, and report layout |

`ROADMAP.md`, `TASKS.md`, checkpoint notes, and `LOOP_LOG.md` retain the
history of the proof's development. Their earlier open obligations are
superseded by the current status and final checked declarations.

The verified result covers native-length distance and exact rate. It does
not include encoder/transposed-encoder operation counts or the complete
arbitrary requested-length wrapper. Some fixed-occupation intermediate
bounds are more conservative than the paper's; see
[FIXED_NATIVE_SEMANTIC_REVIEW.md](FIXED_NATIVE_SEMANTIC_REVIEW.md).
