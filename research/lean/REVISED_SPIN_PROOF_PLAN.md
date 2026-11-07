# Revised SPIN proof plan

This is the primary working plan for Lean verification of the paper's revised
Scalable SPIN theorem and finite Turbo SPIN distance claims. Update it when a
gate closes, a proposed method fails, or the proof target changes. The paper,
the actual setup definitions, and checked Lean statements determine the claim;
research scripts and this plan do not.

## Working rule: audit before designing

Before starting or redesigning a proof component, inspect the corresponding
part of the completed Structured SPIN development. Record:

1. its theorem statement and probability space;
2. the mathematical reduction and which hypotheses it discharges;
3. how candidate data were generated, represented, and checked in Lean;
4. the soundness theorem and the final theorem connection;
5. observed build-time or memory limits.

Then state what is reused, what changes because the construction changes, and
why the new method is preferable. Consult this plan and the old proof again
when a new obstacle appears. A faster-looking standalone check is not a
substitute for a connection to the sampled encoder or complete parameter
coverage.

## Target and present status

The Scalable target is the paper's one-bit route, shared growing Golay--BAA
outer, and independent uniform `GL(22,2)` updates with the actual 22-bit state
map. The desired result is a rate-half family whose probability of minimum
distance at most `floor(0.11*N)` tends to zero under its original sampled
setup. The finite Turbo target is the separate four-bit packet construction
at the lengths and setup-failure margins claimed in the paper. Do not transfer
the old BCH/IMT theorem or the Scalable theorem to Turbo by name similarity.

For Scalable SPIN, the fixed-occupation and growing-sparse results have Lean
proofs. The dense analytic and selected-outer first-moment bridges are proved
conditional on `ScalableDenseRates` in
`SpinCodes/UniformMixer/ScalableDenseRateDefs.lean` and
`SpinCodes/UniformMixer/ScalableDenseConditionalFinal.lean`. The sound
fixed-point vertex checker is in `ScalableDenseVertex.lean`. The
`ScalableDensePilot.lean` theorem checks one rectangle, including its transfer
witness, four vertices, outer support, and pointwise rate conclusion. It does
not establish the full dense cover. The final unconditioned distance, rate,
and finite Turbo theorems remain open in this revised development.

## Audit of the completed proof

| Previous component | Checked method and lesson | Revised counterpart |
| --- | --- | --- |
| `SpinCodes/Cover/` and `scripts/emit_cover.py` | An untrusted search emits an exact split tree and rational data. Lean proves the checker sound and checks every leaf. A monolithic 13,513-leaf check exhausted 16 GB; 44 independent parts worked. | Use chunked checked data for the 46,317-box Scalable dense cover. Never emit one giant term or one theorem per vertex without a resource trial. |
| `SpinCodes/Structured/DenseOccupationGeometry*.lean` | Rational rectangles and a split tree are checked independently of the numerical family. `global_cover` proves coverage of all real parameters. | Reuse the `Rect`/`Tree` definitions and soundness argument if their performance scales; partition into independently checked subtrees if needed. |
| `SpinCodes/Structured/DenseFourierExactW*`, `DenseFourierExactCheck.lean`, and `DenseFourierExactB*` | The transfer/Collatz inequalities use exact integer cross multiplication with rational parameters and integer vectors; only the logarithm/KL vertices use fixed-point intervals. Witnesses are shared across boxes; a uniform prefactor is established before assembly. | Use the three-state `paperSpectrumTransfer` and its actual production spectrum. Preserve exact transfer checking, rather than forcing tiny matrix entries onto the `10^30` fixed-point grid. The old seven-coordinate IMT matrix and divisor 128 do not apply; the new exponent has divisor 64. |
| `SpinCodes/Structured/DenseGeometryRate.lean` and `DenseOccupationMixedCertified.lean` | Checked boxes, exact geometry, and per-family bounds become a full `DenseRates` theorem. | Prove `ScalableDenseRates η C` with a positive uniform `C` and a checked `η`, then feed the existing conditional dense assembly. |
| `SpinCodes/Structured/ConcreteNativeTheorem.lean` and `SpinCodes.Native` | The proof pays the selected-outer failure once and connects first moments to the original product setup law and minimum distance. Final declarations are audited for recursive axioms. | Adapt the generic first-moment and selection argument to `scalableSetupLaw`; do not import the old IMT encoder conclusion as the revised result. |

The old `FINAL_REPRODUCTION.md` records the build and axiom standard. Its
Python evidence is diagnostic; Lean's checked statements are authoritative.

## Dense certificate plan

The current external source is
`../workstreams/turbo_spin_migration/shared_inner_dense_certificate.py`.
It reports 46,317 accepted rational boxes over
`[1/10000,1] × [13/125,112/125]`, using the paper's refined outer
supports. It can report its hardest box and a rational candidate with
`--report-hard-box`; it does not yet emit the complete certificate. Its
floating-point optimizer only proposes witnesses; its Arb arithmetic checks
them externally. Neither result is a Lean proof.

### Gate 1: hard pilot and exact representation

Select a box near the reported worst outward bound, below
`-1.00147e-7`. Emit its exact outer line, endpoints, `p`, `y`, rational
`z`, rational Collatz radius `ρ`, and positive three-state vector. The
external search uses `z = exp(-λ)`; replace this with a rational `z` and
recompute all bounds. Check this hard box in Lean at target exponent
`η = 1e-7`, including support membership and the actual transfer. If
rational rounding erases the margin, increase precision or refine that box.

**Closed locally:** `ScalableDenseHardPilot.lean` checks the reported hardest
box at `η = 1e-7`. The first rational radius failed the fixed-point transfer
check by two units at scale `10^30`. Raising its numerator by three units
passes all three rows and all four exponent vertices. The box also reaches
the selected-outer point-rate interface with the actual support line. The
recursive axiom audit lists only Lean's `propext`, `Classical.choice`, and
`Quot.sound` for that theorem.

### Gate 2: reusable sound leaf checker

Generalize the pilot's proof into one leaf interface. The Boolean checks
must imply: (a) the production-spectrum transfer is bounded by `ρ` using
the supplied positive vector; (b) all four exponent vertices are at most
`-η`; (c) the indicated affine line belongs to the refined outer
majorant; and (d) the vector gives a specified finite moment prefactor.
Use the existing fixed-point log and KL soundness theorems. Test a shared
witness on several boxes before generating all leaves. The generator is
untrusted: Lean must check every numerical and semantic condition.

**Soundness bridge proved:** `ScalableDenseTransferEvalDefs.lean` computes
outward intervals for the exact production-spectrum three-state matrix;
`ScalableDenseTransferEval.lean` proves each interval and Boolean row check
sound; `ScalableDenseLeaf.lean` combines the transfer, four vertices, and
uniform prefactor into a parameterized box theorem. The hardest box has been
instantiated through this generic theorem. Support membership is checked
separately against `Majorant.refined.supports`, as in the old proof. A second
pilot with two adjacent boxes from the same witness bank now checks its
shared transfer and eight vertices in `ScalableDenseSamplePilot.lean`.
That module takes about 37 seconds to compile, versus about 33 seconds for
an import-only module on this machine. This is a small resource trial, not
an estimate for the full cover. The generator should cache transfer checks
by witness and vertex checks by the full witness/coordinate tuple; adjacent
boxes can share corners. Test a larger chunk before fixing module sizes.

**Correction after a mixed-witness test:** the first fixed-point transfer
checker works on the hard and first shared-bank witnesses, but a fallback
witness has a coordinate about `8.6e9`. Rounding a tiny transfer entry at
`10^30` creates a row error about `5.2e-21`, too large for this box's
`1e-7` exponent margin. All four exponent vertices still pass. A closer
inspection of `DenseFourierExactCheck.lean` showed that the old completed
proof avoided this issue by exact integer cross multiplication for Collatz.
The initial audit recorded rational witnesses but missed the importance of
that *exact* arithmetic. The next checker should copy that separation of
responsibilities: exact rational transfer, fixed-point transcendental
vertices. Do not compensate for this failure by widening `ρ` until the
exponent fails.

**Exact replacement proved:** `ScalableDenseTransferExactDefs.lean` evaluates
the production three-state matrix in `ℚ`, and
`ScalableDenseTransferExact.lean` proves its three-row Boolean check sound
against `paperSpectrumTransfer`. The generic `ScalableDenseLeaf.lean` now
uses this checker, while leaving KL/log vertex checks on the old fixed-point
path. `ScalableDenseMixedPilot.lean` checks a fallback witness with the large
coordinate and reaches `ScalableDensePointRate` with the actual outer support.
The hard pilot also passes through the exact generic leaf. Both point-rate
theorems list only standard Lean axioms.

### Gate 3: complete geometry and resource control

Have the generator retain the exact split decisions and accepted boxes.
Use the previous `check_sound`/`global_cover` pattern so every real
parameter in the full rectangle reaches a checked leaf. Separate numerical
leaves from geometry where that reduces duplication. Split the output into
bounded, independently compiled modules; test one representative chunk's
time, memory, and `.olean` size before selecting a chunk limit. A checked
box count alone does not prove coverage.

### Gate 4: uniform rate and final Scalable theorem

Assemble all leaves with a single `η = 1e-7` and a finite positive
prefactor `C`, proving `ScalableDenseRates η C`. Apply
`scalableGoodLayer_total_tendsto_of_dense_rates`. Connect the selected
first moment to the original product law in `ScalableSetup.lean`, charge
the bad outer event once, and prove the minimum-distance statement for the
realized code. Check injectivity and exact rate for every setup. Pin the
full theorem statement and recursive axioms; only `propext`,
`Classical.choice`, and `Quot.sound` may remain.

## Finite Turbo gate

Before designing its Lean checker, audit the old finite certificate and
the current Turbo replay code against the actual RS8 outer, packet route,
and uniform state-map law. Record the exact first-moment summands, one- and
two-group exceptions, reference-outer domination for larger occupations,
directed arithmetic, and how the sum is connected to a setup event. Then
choose a checked-data format and prove its soundness before generating all
three length-specific certificates. The Scalable dense checker may supply
fixed-point primitives, but its asymptotic cover does not establish a finite
Turbo margin.

## Progress and decision log

- 2026-11-07: Audited the old split-tree cover, generated Fourier witnesses,
  geometry checker, dense-rate assembly, final native theorem, and resource
  failure recorded by `emit_cover.py`. The revised dense pilot compiles
  without `sorry` or custom axioms. The complete cover and final theorem do
  not yet compile because they have not been constructed.
- 2026-11-07: Added an optional hard-box report to the existing Arb search.
  The 256-bit replay retained its 46,317 accepted boxes and original worst
  bound. It located a box near `α = 1`, `x = 0.5035`. A rational `z`, `ρ`,
  and positive three-state vector still give an outward Arb exponent below
  `-1e-7` for its four corners. This is candidate data, not yet a Lean check.
- 2026-11-07: Built and compiled the exact transfer evaluator, its soundness
  theorem, a generic leaf theorem, and the hard-box instance. The external
  candidate was adjusted by `3e-30` in spectral radius to account for
  directed fixed-point rounding. The hard box is Lean-checked with the
  `1e-7` exponent margin and reaches `ScalableDensePointRate`.
- 2026-11-07: The first two accepted boxes use one shared witness; their
  transfer and eight vertices pass Lean. A small timing control showed about
  four seconds beyond module import. The 98 initial witness-bank radii range
  from about `1e-16` to `1`, within the `10^30` fixed-point resolution.
- 2026-11-07: A second, locally optimized witness exposed a fixed-point
  transfer failure from a large Collatz coordinate despite passing vertices.
  Re-audited the old exact Fourier checker and found its integer
  cross-multiplication strategy. This is now the preferred transfer path.
- 2026-11-07: The new exact rational checker and its soundness theorem
  compile. The mixed witness now gives a fully checked point-rate rectangle;
  the generic leaf and hard pilot also compile against the exact predicate.
- Next action: make the untrusted generator emit shared exact witnesses,
  vertex data, and the exact split tree; compile a larger mixed-witness
  chunk to choose resource limits. Then assemble the full Scalable dense
  cover and final sampled-setup theorem. Keep finite Turbo as its separate
  gate.

At each future checkpoint, update the status and next action here before
reporting completion or changing strategy. Preserve failed candidates and
their reason in research notes when they inform a later decision.
