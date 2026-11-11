# Turbo SPIN Lean handoff — 2026-11-08

## Read this first

Continue in `C:/Users/peter/.codex/worktrees/316b/permute_conv`, branch
`codex/turbo-spin-results`. The latest commit is `89b5418a`
(`Close native Scalable claims and establish Turbo proof interfaces`).
The later proof work is saved in the working tree, largely as untracked
Lean modules and scripts. It has **not** been committed or pushed. Do not
clean untracked files or assume another checkout contains this work.

The user cleared the automated goal and requested this handoff for another
session. Do not restart an autonomous goal without their instruction.
The immediate checkpoint is the authenticated local Turbo transfer family;
the end-to-end K18 and K20 distance certificates are **not closed**.

The primary planning document is `REVISED_SPIN_PROOF_PLAN.md`. Consult it
periodically and maintain it. The user's mandatory rule is: before changing
proof strategy or starting a component, audit the corresponding old proof,
then explain what is reused and what differs. The old Structured proof used
untrusted witness generation plus small Lean kernel checks, not enormous
all-in-Lean computation. Do not forget that history again.

Use the Controlled Writing for Cryptography skill for mathematical prose:
`C:/Users/peter/.codex/skills/controlled-writing-for-cryptography/SKILL.md`.
Never run concurrent benchmarks. No benchmark is needed for this task.

## Exact objective and status

Close the actual ideal sampled-setup, 22-state Turbo SPIN theorem, at rate
1/2 and 10% binary distance, for these endpoints:

| Endpoint | Input dimension | Output length | Bad-weight cutoff | Required margin |
| --- | ---: | ---: | ---: | ---: |
| K18 | 262144 | 524288 | 52428 | 58.24 bits |
| K20 | 1048576 | 2097152 | 209715 | 51.82 bits |

The cutoff is **floor(N/10)**. Do not change it to a ceiling. Do not quietly
weaken the rounded margins to 40 bits. No discrepancy currently requires
that concession. K16, C++ refinement, new performance experiments, and
paper edits are outside this task. Natural lengths are acceptable.

Scalable SPIN is already closed in
`SpinCodes/UniformMixer/ScalableCertified.lean`: original-setup failure
probability for distance at most floor(0.11 N) tends to zero; binary linear
injection, exact rate 1/2, natural schedule, and ordinary/adjoint circuits
of at most 272 N coefficient-AND/XOR gates. Setup generation, C++ refinement,
and a new clean dependency replay are not included in that claim. Do not
reopen Scalable to avoid finishing Turbo.

External outward-arithmetic diagnostics pass Turbo's targets:

| Contribution | K18 margin | K20 margin |
| --- | ---: | ---: |
| Complete q=1 calculation | 58.2456338686 | 56.2570901999 |
| Complete q=2 calculation | about 105.37 | 102.6824180241 |
| q>=3 tail | 66.4739187234 | 100.5938735669 |
| Total, including applicable rare event | 58.2408311085 | 51.8224664170 |

These are **diagnostics, not Lean endpoint certificates**. Mathematical
evidence is favorable; substantial finite arithmetic and assembly remain.
Conditional theorem interfaces do not count as closing either endpoint.

## Checkpoint reached at handoff

The local-family pipeline now connects checked integers to the actual
physical kernel. Important completed modules in `SpinCodes/Turbo/`:

- `PhysicalProfile22`, `PhysicalProfileNumeric*`: actual expansion-map
  census of 25 packet-weight profiles, fixed interval arithmetic and
  authenticated class maxima/means. The uniform-live class is index 2;
  it requires a mean, not a maximum substituted inconsistently.
- `PhysicalPacketMoment22`, `PhysicalPacketFourier22`,
  `PhysicalLayerPolynomial22`, `PhysicalInputPolynomial22`,
  `PhysicalCharacterPolynomial22`, `PhysicalFeedbackFourier22`,
  `PhysicalProfileMoment22`: exact physical input moments and zero-feedback
  Fourier formula agree with the numeric polynomial banks.
- `PhysicalSmallBirth22Census`, `PhysicalSmallBirth22Polynomial`,
  `PhysicalSmallBirth22Budget`, `PhysicalSmallBirth22BudgetSmoke`:
  actual 240 one-packet and 27000 two-packet input census, normalization,
  and checked concrete birth cones; empty input is exact as well.
- `PhysicalProfileFamily22`, `PhysicalProfileFamilyExtract22`,
  `PhysicalProfileFamilyFourier22`: cached integer checks imply maxima,
  uniform means, Cauchy return bounds, and actual zero-hold bounds.
- `PhysicalProfileRows22`, `PhysicalProfileAssembly22`,
  `PhysicalProfileCompleteRows22`: integer checks plus authenticated birth
  cones give actual physical and selected-step row soundness. The latter
  uses the already proved binomial thinning. Nonnegative matrix entries
  are also checked.
- `PhysicalProfileBirthFamily22`, `PhysicalProfileBirthFamilySmoke22`:
  all 17 birth cases, exact for j=0,1,2 and coarse for j>=3.

All these complete-module checks passed; printed declarations use only
`propext`, `Classical.choice`, and `Quot.sound` (or a subset).

The final wrapper is `PhysicalProfileFamilySmoke22.lean`. It exposes
`ProfileFamilySmoke22.physicalRows_sound`, `selectedRows_sound`, and the
corresponding nonnegativity theorems, with **no numerical or birth premises**.
Its final compilation status is recorded in the closing checkpoint below.
Despite the historical name “Smoke,” it covers every occupation and row
for one rational local family; it is not just a spot check.

That family uses z=8607079764250578/10^16 and integer endpoint denominator
2^80. The rational z is an arbitrary admissible witness, approximately
exp(-3/20); **do not claim equality with that exponential**. Its complete
Boolean certificate is `PhysicalProfileFamilySmoke22Check.familyChecked`.

Critical numerical convention: `returnCaps[j][i]/endpointDen` is already
the **final zero-column entry**, divided by 4194303. The checker establishes
min(M,sqrt(S F))/4194303 <= that value. Never divide it a second time.
For j>=3 births, the row places zero-hold Z in coordinate 0 and a total
moment cap in coordinate 1. Overlap is permitted: these bound different
pieces, not disjoint estimates of the same probability.

## Sparse composition: exact next blocker

These modules passed complete checks:

- `SparsePhysicalKernel22`: actual shared-permutation regional kernels
  agree with normalized occupancy coefficient DPs.
- `SparseRegionalEnvelope22`: propagate one continuous state through all
  64 regions. There is **no state reset at region boundaries**.
- `SparseSpectrum22`: all current RS8 coefficients Cv, including every
  shell and the zero coefficients outside 5..64. Do not import RS16 data.
- `SparseActiveProfiles22`, `SparseActiveKernel22`: remove inactive groups,
  reindex distinct active groups, preserve one shared permutation, and
  obtain exact occupancy-0/1/2 regional kernels.
- `SparseActualMoment22`, `SparseActualEnvelope22`: original shared-setup
  support measure and support-dependent Chernoff/envelope families.

**Update (later 2026-11-08): `SparseShellFold22`, `SparseActualShell22`,
`SparseMarkerFold22`, `SparseOccupation22`, and `SparseQ1Choice22` now pass
complete-module checks with standard axioms.** See the plan section
"sparse fold checkpoint". The historical failure notes below are kept for
reference.
Loop tick 1: `SparseSpectrumInteger22`, `SparseFixedDefs22`, and `SparseFixed22`
are also GREEN; symbolic q2 marker assembly already present. Next: soundness
for the remaining SparseFixedDefs22 checks (poly product, normalization,
q1 step/initial, power trace).
Loop tick 2: done — `SparseFixedSound22` GREEN. Next: binary-powering chain
for the regional occupancy polynomial, then K18 producer review.
Loop tick 3: `SparseRegionalPower22` GREEN; the three K18 families were
generated (all precision gates pass) and all 18 per-tilt modules are GREEN
(`ProfileFamilyK18Tilt{1,2,3}22`). Next: row-integer adapters for j≤2 and the
regional/q1/q2/power certificate producer.
Loop tick 4: `SparseRowBase22` GREEN. Kernel probe showed flat-Array checks
exhaust kernel memory at endpoint size; list-structured checks pass
(≈12 s per 65-shell step). New GREEN list pipeline: `SparseListDefs22`,
`SparseListSound22`, `SparseListQ1Regional22`. Use these, not the Array
checks, for endpoint data. Next: K18 certificate producer + Check modules.
Loop tick 5: old-proof memory audit recorded (small per-item kernel
decisions in many generated modules). K18 certificates generated at scale
2^384 into ignored `SpinCodes/Turbo/Generated/K18/` (builder
`.lake/turbo/build_k18_cert.ps1`; regenerate with
`scripts/emit_turbo_k18_sparse_cert22.py` then
`scripts/emit_turbo_k18_q1_assembly22.py`). All 54 certificate modules and
the three `Tilt{k}Q1` assemblies are GREEN. Producer q1 = 58.2456338686 bits.
Next: the K18 q1 endpoint theorem.
Loop ticks 6-8: `SparseRS8Counts22` GREEN (exact integer counts). K18 q1 CLOSED
(`Spin.Turbo.K18Q1.q1_le`, 2^-58.24563) and K18 q2 CLOSED
(`Spin.Turbo.K18Q2.q2_le`, 2^-104.18), both in ignored generated modules with
tracked templates/emitters under `scripts/`. Remaining for K18: q ≥ 3 tail
(needs ≤ ≈ 2^-66.25) and the final sum.

(Historical) `SparseShellFold22.lean` was failing, and its later dependents
were drafts. The existing error log is
`.lake/turbo/SpinCodes_Turbo_SparseShellFold22.lean.log`.

At last verification:

- Lines 47 and 66: default recursion limit during finite-sum rewriting.
- Lines 84, 95, 107: broad matrix/sum simplification reaches the same limit.
- Lines 148/150: `simp_rw` leaves a matrix product application that is not
  type-correct at implicit transparency; the subsequent fold rewrite fails.

No resource cap was raised. Likely repair: extract opaque scalar/matrix
sum identities, avoid broad simplification over `Finset (Fin 64)`, and use
explicit `congrArg` rather than automatic `congr`/unfolding. Check actual
instances before attributing errors to difficult mathematics.

After ShellFold, verify `SparseActualShell22`, `SparseOccupation22`, and
`SparseQ1Choice22`, plus whichever marker/spectrum-integer dependencies
their imports require. `SparseFixedDefs22`, `SparseFixed22`,
`SparseSpectrumInteger22`, and `SparseMarkerFold22` are not to be treated
as checked simply because their source exists.

Two essential q1 details, already reflected in drafts:

1. The diagnostic selects a different tilt **for each outer-weight shell**.
   Apply that choice before the unconditional support sum. Taking a minimum
   after summing different bounds is not the same calculation.
2. Shells 0..7 use the trivial probability bound 1 in the saved diagnostic.
   Cv5, Cv6, Cv7 are nonzero and must be charged. Do not force those shells
   through an unsuitable Chernoff bound or silently omit them.

Use the complete unconditional 65-shell q1 and complete q2 calculation for
both sizes. K20's full q2 bound is already adequate; no heavy-specific q1
or q2 transfer is necessary. Charge K20's rare event once, where the tail
truncation requires it.

## Numerical work queued, but not run

The profile producer is `scripts/emit_turbo_profile_family22.py`. It now
accepts `--z-num`, `--z-den`, and `--namespace`; `--emit-only` preserves the
existing witness and rejects a conflicting z. This CLI extension is saved
but the next batch has not been run.

Draft `scripts/emit_turbo_profile_k18_families22.py` prepares precisely three
K18 tilts: lambda=1/1250, 3/2500, 1/625. The q2 witness reuses 1/625.
Select dyadic rational z near exp(-lambda), denominator 2^80. Use that exact
z consistently in the family and endpoint proof. The sparse 2^512 grid
then represents it exactly as `baseInt=zNum*2^432`.

No K18 batch JSON, Lean witness data, or new numerical replay was generated.
Review the draft's exact precision diagnostics before running it. All
nonzero first/second moments at these tilts are large; endpoint rounding
should be tiny relative to them, but preserve structural zeros exactly.
Only one agent/process should own numerical experiments.

Suggested order:

1. Fix ShellFold and close the symbolic original-setup q1/q2 assembly.
2. Generate and kernel-check the three K18 local families using the now
   verified concrete-family pattern. Birth checks accompany every family.
3. Authenticate low-tilt regional DP and scalar endpoints with materialized
   integer traces. Authenticate all 65 shells and the complete q2 sum.
4. Authenticate K18 q>=3 witnesses and their full occupation cover.
5. Add the final K18 probability sum at 58.24 bits, without leftover
   numeric hypotheses; audit axioms and original-setup statement.
6. Reuse the pipeline for K20 at 51.82 bits. Do not start K20 by duplicating
   or replacing the checked K18 strategy.

## Tail arithmetic design and historical audit

The absolute 2^512 lower-power trace in `SparseFixedDefs22` is only for
low-tilt q1/q2 endpoints. A positive final endpoint is mandatory. All
retained q1/q2 lambda*d values are below 168 (below 243 binary exponent
bits); the sufficient first-three K18/first-two K20 families need fewer
than 122 bits. High-tilt tails cannot use this absolute grid.

Historical sources already audited: `Numeric/FixedLog.lean` (`flog_mem`,
`flogI_mem`, range reduction/Horner), `Numeric/FixedDefs.lean`, and
`Structured/DenseFourierExactCheck.lean` (integer Collatz checks).
Retain the split between exact integer transfer checks and checked logs.

Draft `MatrixScaledDefs22` / `MatrixScaled22` provides normalized finite
matrix squaring. **Not verified as a complete component.** The core file
was compiled before its last two Boolean definitions were appended; its
source is newer than its olean. The semantic file has not been compiled.

Use a materialized witness matrix for each stage, a one-stage containment
check, and an induction theorem. Do not ask the kernel to reduce an entire
recursive squaring chain: repeated list/function lookup can recompute
prior stages exponentially. Each positive rational normalization is
tracked symbolically; checked logs account for the accumulated scale.
No tail witness generator or numerical experiment was started.

For K18 dense occupations q=385..2048, the existing replay uses:

- M = sum_j binom(16,j) theta^j (1-theta)^(16-j) Physical[j], theta=15p/16;
- moment = e0 M^8192 1;
- beta=2^256/65535^8, d=52428;
- log prefactor = log choose(2048,q) + q log beta - d log z
  + 64*(-log choose(2048,q)-q log p-(2048-q)log(1-p)).

The exact regional calculation for q=3..384 uses paired coefficient DPs;
q=62..143 uses fractional exponent 4/5 at tilt17/200. The others use the
retained alpha=1 witnesses. Inspect `scripts/replay_turbo_coarse_k18.py`
and `research/workstreams/k18_margin60/replay_fractional.py` before
designing that numerical component. Do not substitute an untested Collatz
bound for the finite power: K18's tail has only about0.23 bits of slack
against the budget needed for its rounded final margin.

Convenient targets: q2 <=2^-80 for both sizes; K20 tail <=2^-80 suffices;
K18 tail needs roughly66.243 bits, so66.4 is a useful target. Final positive
rational P=num/den can be compared exactly to rounded margins with
num^25*2^1456 <= den^25 for K18, and num^50*2^2591 <= den^50 for K20,
or an equivalent verified log calculation. Never trust a floating bit count.

## Where the evidence lives

Ignored diagnostics in this worktree's `research/lean/.lake/turbo/`:

- `q1_coarse_complete65_diagnostic.json`
- `q2_common_marker_diagnostic.json`
- `q2_common_marker_k20_full_recombination.json`
- `k18_coarse_full_diagnostic.json`
- `k20_coarse_extended_diagnostic.json`
- `profile_family_3_20_rational.json`

Original current receipts are in
`C:/Users/peter/repo/permute_conv-github-bch/tmp/`:
`k18-rs8-gf32-s22-v1.json` and `k20-two-ms-rs8-transfer-v1.json`.
Historical files with RS16 in their names supply witness choices only;
they must not supply the current outer spectrum. The verified Lean
`SparseSpectrum22` identifies the actual current Cv.

The retained numeric source is in Lean data modules, not just raw JSON.
Do not commit raw experiment dumps, timing data, or generated PDFs.
Preserve unrelated `SOL_TO_ASTRA_HANDOFF.md`, `bch_spectrum_work/`,
`research/paper/$paperCheckDir/`, and `workstreams/`; they are not cleanup
targets or automatically part of this proof checkpoint.

## Build mechanics and common pitfalls

Local Lean is4.34. Dependency objects live in
`C:/Users/peter/repo/permute_conv-github-bch/research/lean`; local objects
must come first in LEAN_PATH. The ignored helper performs that setup:

```powershell
powershell -NoProfile -File research/lean/.lake/turbo/compile.ps1 -RelativeFile SpinCodes/Turbo/PhysicalProfileFamilySmoke22.lean -MemoryMB 6000
```

Logs are `.lake/turbo/SpinCodes_Turbo_<Module>.lean.log`. Objects go under
`.lake/build/lib/lean/SpinCodes/Turbo/`. Preserve the helper if moving work
to another session/worktree; do not assume `lake build` sees these local
objects or rebuilds a changed dependency automatically.

The helper serializes heavy checks with mutex
`Local\CodexSpinTurbo316bCompile0`. One audited core-only light check can
use `Local\CodexSpinTurbo316bLightCompile0`. Heavy imports can use about
12GB private memory even when RSS is much smaller. Do not start multiple
heavy compilers or raise memory/recursion limits blindly. An elan shim and
toolchain Lean process for one command are not two independent jobs.
Do not touch unrelated Lean4.35/CertifiedDelivery processes.

Use a live process/session to determine whether a queued check is still
running; a stale log is not a new failure. Wait for the complete module's
exit0 and inspect printed axioms. A file with one successful theorem and
later errors is NOT green. A failed declaration can print `sorryAx` even
without a literal `sorry` in the source.

Recent repairs worth reusing:

- Equivalent finite-enumeration or decision instances can make definitional
  equality unfold sums over2^64 inputs. Freeze the polynomial/probability
  and prove an explicit transport at that level, before applying coeff or
  another outer function. See `PhysicalInputPolynomial22` and the final
  `probability_decidable_eq22` helper in `PhysicalProfileRows22`.
- To recover the Fintype carried by P, use a non-instance implicit binder:
  `{fi : Fintype A} (_P : @FinPMF A fi)`. An instance binder `[fi]` can
  eagerly pick the wrong enumeration.
- Avoid `congr 1/2` on large concrete finite objects; use explicit congrArg.
- `&&` is left-associated. Normalize with `and_assoc` or project the actual
  tuple shape. Do not assume h.2.1 after three left-associated checks.
- `ofInt_mem 1` needs `simpa only [Int.cast_one]` for a real1 target.
- Partial functions can resist simp; the small birth smoke branches closed
  with explicit `change` to the intended function and existing theorem.

## Stopping point

All three Sol agents stopped. No new numerical batch was launched.

**Final verified checkpoint:** `PhysicalProfileFamilySmoke22.lean` compiled
with exit0. Its `physicalRows_sound`, `selectedRows_sound`,
`physicalRows_nonneg`, and `selectedRows_nonneg` declarations all print only
`propext`, `Classical.choice`, and `Quot.sound`. The complete-row and
all-birth wrappers also passed exit0 before it. No task-owned Lean process
remained in the final process inventory. All root compilation sessions
were terminal; no further numerical or compilation work is queued by this
session. The working tree and this handoff are saved but uncommitted.

Recommended next-session request: “Read TURBO_K18_K20_HANDOFF.md and the
primary plan. Resume K18: audit and fix SparseShellFold22, finish the sparse
assembly, then authenticate the three retained low-tilt families and their
endpoint arithmetic. Keep K20 and the full goal in scope.”

## 2026-11-09 checkpoint: K18 endpoint closed on Peach

`SpinCodes/Turbo/Generated/K18/Endpoint.lean` (generated by
`scripts/emit_turbo_k18_endpoint22.py`) proves

```lean
theorem Spin.Turbo.K18Final.k18_endpoint (points : Fin 8 → Field16)
    (hp : Function.Injective points) :
    (turboSetup22 128 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 128 points).Z 52428 ω) ≤ (2 : ℝ) ^ (-(58.24 : ℝ))
```

i.e. for the ideal sampled 22-state setup at R = 128 (N = 524288,
k = 262144), the probability that some nonzero message has output weight
≤ floor(N/10) = 52428 is at most 2^-58.24. Axioms: `[propext,
Classical.choice, Quot.sound]` (Peach build log). No numeric hypotheses:
q=1 (`K18Q1.q1_le`), q=2 (`K18Q2.q2_le`), the q∈[3,2048] tail
(18 term modules, all standard axioms) and the final integer inequality
`(2048·aSum + b2 + Σ A q)^25 · 2^1456 ≤ 2^6400` are all kernel-checked.

How it was built (see REVISED_SPIN_PROOF_PLAN.md for the audit trail):
- Heavy builds on Peach (`~/turbo/pbuild.py`, 14 jobs). The 1450-module
  closure's sources are byte-identical to this worktree modulo CRLF
  (sha256 comparison, 2026-11-09).
- Tail squaring chains use the packed Kronecker check
  (`SparseListPackedDefs22` core-only defs, `SparseListPacked22`
  soundness). Lean keeps kernel memory across declarations of one process,
  so each squaring step is split into 11 per-row modules
  (`Tail{k}ExStep{s}R{i}`), each with 11 per-entry declarations; the
  fractional family's pair/base checks are split into `Tail7FrRpow{a}` and
  `Tail7FrBaseC{c}` (`scripts/patch_turbo_k18_frac7_split22.py`, still to be
  folded into the main generator).
- Marker bounds are the constant 1 (`aMk`), avoiding a 730-arm `match`.

Remaining: local olean audit of K18 (in progress), then K20 with the same
pipeline (R = 512, d = 209715, target 51.82 = 2591/50 bits).

K18 LOCAL AUDIT PASSED (2026-11-09): the 1450-module olean closure was
pulled from Peach into this worktree's `.lake/build`; `SpinCodes/Turbo/K18Audit.lean`
compiles locally and prints
`'Spin.Turbo.K18Final.k18_endpoint' depends on axioms: [propext, Classical.choice, Quot.sound]`.
K18 is closed. Remaining: K20.

## 2026-11-09: K20 endpoint closed; both endpoints audited locally

```lean
theorem Spin.Turbo.K20Final.k20_endpoint (points : Fin 8 → Field16)
    (hp : Function.Injective points) :
    (turboSetup22 512 points).joint.prob
      (fun ω => 1 ≤ (turboSetup22 512 points).Z 209715 ω) ≤ (2 : ℝ) ^ (-(51.82 : ℝ))
```

N = 2097152, k = 1048576, cutoff floor(N/10) = 209715. Axioms:
`[propext, Classical.choice, Quot.sound]`, printed both on Peach and by the
local audit `SpinCodes/Turbo/K20Audit.lean` (which also re-prints K18).
The local audit uses oleans pulled from Peach; the 3094-module K20 closure
sources are identical to this worktree modulo CRLF (sha256 comparison).

Certified K20 value: about 2^-56.257 (q1 2^-56.25709, q2 2^-101.59,
tail 2^-100.594), well inside 51.82 bits. Design notes:
- No rare-event term: the K20 tail uses the same untruncated physical model
  as K18, so BadSplit22 (q1 + q2 + tail) suffices (strictly sound).
- General Lean machinery added: `SparseRpowAB*`, `PairedFracAB*`,
  `TailTermFracAB22` (`tail_term_fracAB_le`, `..._shift_le` with the
  domination shift), `TailMarkerFracAB22` (fractional markers).
- Tail = 32 pieces: exact q=3..315 (8 pieces), markers q=316..8192
  (24 pieces; q=316..615 replanned via `emit_turbo_k20_tail22.py plan`
  because the diagnostic's degree-615 exact witness had only 20 bits and hit
  a precision floor).
- Kernel-depth workarounds without raising limits: split base equalities and
  lengths proven by `decide +kernel`; the 8190-term tail sum split into
  17 chunks of 500 joined by `Finset.sum_Ico_consecutive`.
Producers: `emit_turbo_profile_k20_families22.py`, `emit_turbo_k20_*.py`
(sparse cert, q1 assembly/endpoint, q2, tail, tail_lean, endpoint).

Status: BOTH ENDPOINTS CLOSED. Nothing committed. Remaining housekeeping
(optional): fold `patch_turbo_k18_frac7_split22.py` and the K18 per-row
split into the K18 generators; generated trees under
`SpinCodes/Turbo/Generated/` are git-ignored and must be regenerated from
the scripts + `.lake/turbo` caches.

## 2026-11-09: margins tightened to the certified values

The endpoint statements now carry the full certified margin:
- `K18Final.k18_endpoint`: ≤ 2^(-58.24083)  (certified total 2^-58.2408311)
- `K20Final.k20_endpoint`: ≤ 2^(-56.25709)  (certified total 2^-56.2570902)
Both GREEN on Peach and in the local audit (`K18Audit.lean`, `K20Audit.lean`),
standard axioms only. The final step uses `nat_div_two_pow_le_two_rpow_neg`
(EndpointBits22) with small-denominator exponents 84158/1445 (K18) and
91249/1622 (K20), then rpow monotonicity to the decimal target: a 10^5
denominator exceeded the kernel's maximum numeral size, which we did not raise.
Both margins are q=1-limited; better tail witnesses could add at most
≈0.005 bits at K18 and nothing measurable at K20.

Paper updated (research/paper): theorem λ = (61.88, 58.24, 56.25); K20
rare-event paragraph removed (q=1,2 sums include 5..8-packet supports);
appendix explains exact-coefficient vs marker witnesses and cites the Lean
development as the authority for K=2^18, 2^20 (K=2^16 remains numerical);
implementation table 58.240 / 56.257; abstract and intro "above 56 bits";
AI disclosure mentions the Turbo Lean development. Paper builds cleanly.

## 2026-11-10: s22 ladder G14..G20 certified and audited

`Spin.Turbo.G{14,16,18,20}Final.g{e}_endpoint` (turboSetup22 with L/16 =
8/32/128/512, cutoff floor(N/10)) prove 2^-64.25152, 2^-62.34042,
2^-60.36266, 2^-58.36785; local audit `SpinCodes/Turbo/LadderAudit.lean`
prints standard axioms for all four. G18/G20 supersede the K18/K20 endpoints.
Producers: `scripts/ladder22.py` (q1/q2) and `scripts/ladder_tail22.py`
(plan/families/pairs/exact/marker/emit/endpoint), using the marker range lemma
`TailMarkerRange22.lean`. G22 is planning; G24/G26 need compact certificates
(see REVISED_SPIN_PROOF_PLAN.md). Nothing committed.

## 2026-11-10: G22 certified and audited

`Spin.Turbo.G22Final.g22_endpoint` (turboSetup22 2048, cutoff 838860) proves
2^-56.3691; standard axioms on Peach and locally (`SpinCodes/Turbo/G22Audit.lean`,
kept separate from `LadderAudit.lean` for memory). Tail: 13 exact pieces (q 3..545)
+ 9 marker ranges (546..32768); exact term modules emitted with
`TURBO_QCHUNK_EXACT=6` (30 per module exhausts memory at G22). Generated trees
(~0.6 GB source, 9.5 GB oleans) are git-ignored. Next: G24 (exact to 1299, markers
from 1300) and G26 (exact to 4999, markers from 5000); see REVISED_SPIN_PROOF_PLAN.md.

## 2026-11-10: G24 certified (Peach), audit pending

`Spin.Turbo.G24Final.g24_endpoint` (turboSetup22 8192, cutoff 3355443) proves
2^-54.36978, standard axioms. Key fixes for K = 2^24: planner marker normalizer
(C^-63, was C^-62 — exact pieces now only to q 123), `BinomFast22.lean`
(kernel-cheap C(N,q); N!/(q!(N-q)!) cost ~18 GB per term), marker MBase modules
`noncomputable section`, exact term modules TURBO_QCHUNK_EXACT=4, top-band marker
module with LEAN_NAT_MAX_SIZE=2^30 (≈134M-bit numerals; value/8 bits). Generated
source 265 MB (git-ignored).

## 2026-11-10: G26 certified — all seven rungs proven (s22)

`Spin.Turbo.G26Final.g26_endpoint` (turboSetup22 32768, cutoff 13421772) proves
2^-52.32161 with standard axioms. Ladder: 64.25152 / 62.34042 / 60.36266 / 58.36785 /
56.3691 / 54.36978 / 52.32161 bits at K = 2^14..2^26. Scoped build limits beyond the
standard -M 14000 are listed in REVISED_SPIN_PROOF_PLAN.md (G24 TailP15MQ0; G26
TailP18-21MQ0). Remaining: Peach audit (leanchecker + fresh regeneration), regenerate
entry points, paper update. Nothing committed.
