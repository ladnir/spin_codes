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
not establish the full dense cover. `ScalableDenseLeaf.lean` supplies the
general leaf theorem; `ScalableDenseCoverCore.lean` supplies exact rectangle
composition. Generated 8-, 32-, 64-, and 128-leaf research chunks compile,
but no theorem covers the full 46,317-box domain. The final unconditioned
distance, rate, and finite Turbo theorems remain open in this revised
development.

For Turbo SPIN, Lean now checks the value-carrying packet route, its
22-state inner moment identity, the generic field-scalar action, and the
RS8 information-set, shortening, and exact-support identities. In
particular, zeroing `m ≤ 4` selected columns leaves
`|F|^(8*(4-m))` outer groups, and inclusion--exclusion now gives a
checked exact count for every specified output-column support. Lean also
checks the paper's closed `A_h` expression for the actual number of
messages with `h ≥ 5` active columns, assuming the base field has 16
elements and eight distinct evaluation points. The proof handles the
zero message's empty support separately. It does not yet identify the
concrete GF16/GF(2^32) implementation, establish the complete expected
packet measure, or close any finite numerical endpoint.

## Audit of the completed proof

| Previous component | Checked method and lesson | Revised counterpart |
| --- | --- | --- |
| `SpinCodes/Cover/`, `SpinCodes/Numeric/BAEvalDefs.lean`, `SpinCodes/Numeric/BAEval.lean`, and `scripts/emit_cover.py` | An untrusted search emits an exact split tree and rational data. Lean proves the checker sound and checks every leaf. A monolithic 13,513-leaf check exhausted 16 GB; 44 independent parts worked. The old `LogTree.check` verifies each fixed-point logarithm once; `treeLog_oracle` proves lookup sound even if the search tree is unordered (a bad tree only fails to find a key). | Use chunked checked data for the 46,317-box Scalable dense cover. Share expensive KL and log evaluations before choosing a bulk module size. The old table has one integer key because all its logarithm denominators equal `scale`; the new checker needs rational pairs or, more simply, cached KL terms and witness logs. A new cache must carry its own soundness theorem. |
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
`--report-hard-box`; `--cover-output` now emits the full proposed split tree,
boxes, and witnesses. Its
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

### Audit of the finite proof and correspondence to the old development

The old **Lean** `SpinCodes/Framework.lean` proves the first-moment identity
for a product setup law and Markov's bound for the event that some nonzero
message has low output weight. `SpinCodes/FiniteMoment.lean` also gives the
generating-function tail inequality. The completed Structured SPIN theorem
shows how to connect a checked numerical bound to an actual sampled encoder,
not merely to an abstract matrix. Reuse these probability and conclusion
interfaces, while defining the Turbo outer, packet route, and state maps
anew. The old Structured BA/Golay enumerator and BCH--IMT finite appendix
describe different distributions; neither is a Turbo certificate. The
retained **Turbo reference** certificates for the RS[16,8] outer are external
Python/Arb receipts and independent audits, not Lean proofs. Authentication
of their source hashes and exact dyadic sums does not by itself prove their
numeric endpoints or link them to Lean's sampled setup.
The K20 independent audit explicitly says it authenticates the full-state
census and matrix endpoints without regenerating them. Its local receipts
contain 17 outward physical operators, each `11×11`, derived from all
`2^22` states. K18's `replay_fractional.py` thins these operators, takes
entrywise fractional powers only after the chronological two-step product,
then raises a regional operator to 64 regions. K16's fresh sparse producer
similarly computes state operators and route averages with Arb. Therefore
checking only receipt hashes, `C_v`, and the final dyadic union would leave
the main finite probability inequality as an unsupported premise. A Lean
endpoint checker needs a sound path from the actual map/update law to
outward local operators, then through placement and matrix powers. This
is a distinct gate from exact arithmetic on imported endpoints.

The target is the ideal setup in `research/paper/turbo_spin.tex` and
`turbo_appendix.tex`: independent nonzero GF(2^32) column multipliers,
uniform packet and regional permutations, and independent uniform GL state
updates, sampled once per code. `spin/RS8.md` and
`research/workstreams/turbo_spin_migration/SETUP_AUDIT.md` compare this law
with the executable seed expansion. The executable's SplitMix choices are
not literally the independent uniform law, so the information-theoretic
theorem must not be stated for those finite seeds without another premise.

The local expected packet counts are coefficients of
`Σ_(h=5)^8 A_h (((1+15z)^8-1)/(2^32-1))^h`, with `A_h` from the eight
parallel GF16 RS[8,4] rows. These are expected **message counts**. Each
group has 128 message bits, 256 output bits, and 64 four-bit packets. A
group shuffle makes a support of size `v` uniform over 64 positions, with
uniform nonzero packet labels. The eight random columns are independent;
different messages sharing a setup need not be. Prove this measure and its
normalization in Lean before importing its numerical coefficients.
The existing `UniformMixer/ScalableSetup.lean`, `ScalableRealization.lean`,
and `ScalableFullMoment.lean` already define the uniform 22-state map law,
chronological encoder, and moment bridge. Their one-bit route and shared
Golay--BAA outer cannot be reused for Turbo. No paired 15-state map or
GF(2^32) RS8 outer is currently modeled in Lean. Start the new outer with
an abstract 32-bit field action and an explicit additive-coordinate
equivalence; Mathlib has `GaloisField 2 32`, but the proof must relate this
model to the paper's chosen binary adjoint map. The 22-state Turbo profiles
can reuse the checked state-map law after their 64-bit packet stream is
defined. K16 needs a separate 15-state map and law.

| Length | Sparse part | Larger occupations | Receipt source |
| --- | --- | --- | --- |
| `K=2^16`, paired 15-state map | Fresh one- and two-group support/label folds, including all `v=5..64` and all support pairs | For each `q=3..512`, transfer the authenticated RS16 endpoint by the exact factor `((65535/65537)^4)^q` | `k16_outer_routing/certify_rs8x8_s15.py`, retained `paired-s15-drop10-whole-p256.json` |
| `K=2^18`, exact 22-state map | Full 65 one-group endpoints and 2,145 triangular two-group endpoints, with off-diagonal multiplicity two | Retain every RS16 endpoint for `q=3..2048` using the no-worse complete atom envelope | `rs8_sizes/k18_rs8_certificate.py`, retained `k18-margin60-fractional-whole-p256-v2.json` and sparse receipts |
| `K=2^20`, exact 22-state map | Full one-group fold; globally charge the event that **any** group has a local word of packet support 5..8; compare the remaining two-group measure with the RS16 measure using the exact maximum shell ratio squared | Retain every RS16 endpoint for `q=3..8192` using the no-worse complete atom envelope | `k20_two_ms/rs8_transfer.py`, retained `k20-small-outer-whole-v1.json` and sparse receipt |

For the K20 exceptional event, the two-group comparison is an
**unconditional truncated first moment**, not a setup law conditioned on
its complement. This keeps the direction of the union bound. For all three
lengths, the complete atom bound is `1/(2^32-1)^4`; comparing it with the
RS16 bound `1/(2^16-1)^8` is a pointwise statement about the whole group
measure. Matching only single-column marginals would not transfer the
regional or fractional-moment endpoints. The reference route, inner map,
state continuity, output cutoff, and ordering of the active-group union
must be authenticated in the mathematical bridge, not only in metadata.

The exact bad cutoffs are 13107, 52428, and 209715; the respective
minimum-distance lower bounds outside failure are 13108, 52429, and
209716. The advertised decimal margins 61.884, 58.245, and 51.822 bits
are displays of exact upward dyadic unions. Lean should check a rational
threshold (initially `< 2^-40` at each length), not reason from the
displayed logarithms. The original receipts exist in the canonical
checkout but are ignored and absent in this worktree. Import only the
needed proof data into ignored research build storage, recording hashes;
do not commit raw receipts or treat Python replay as a Lean theorem.

**Finite implementation order.** (1) Pin Lean's actual Turbo setup and
bad event, then prove the RS8 expected outer measure and packet-route
symmetry. (2) Prove a finite first-moment decomposition by active group
count, with the exact one/two-group folds and K20 exceptional event.
(3) Prove the physical 64-bit state operator bound for one selected
`s=22` packet profile before bulk data generation. Audit whether the
existing Scalable 22-state map spectrum/census proofs can supply its
entries; if not, use a checked integer census in the style of the old
Structured `ConcreteCounts` and `FiberCertificate` developments.
The initial inspection finds `ScalableStepLaw.uniform_state_map_update`
reusable for the uniform GL(22,2) transition and the checked `Map22`
coordinates reusable for the fixed expansion. The Scalable image-weight
spectrum is for Bernoulli bit inputs; it does not itself prove the 17
packet-activity operators or their 11-state quotient. Those quotient and
birth-class transitions need an explicit soundness proof.
`SpinCodes/Turbo/PacketRoute.lean` now carries complete packet values under
the old independently sampled row and regional permutations. Lean proves
that taking nonzero-packet support after routing gives exactly
`ConcreteRoute.routeEval` on the original support, hence the old checked
`ConcreteRoute.law` applies. It also proves preservation of active-packet
count. These results use only the standard Lean axioms. They establish the
route interface, not the RS8 outer measure or finite distance bound.
`SpinCodes/Turbo/PacketInner22.lean` serializes 16 four-bit packets per
64-bit step and 64 regions in chronological order. Its ideal law is the
product of the checked packet-route seed law and independent uniform
GL(22,2) updates. `routeInnerMoment22` applies the old
`ScalableEncoderKernel.sampledInputMoment_kernel` to the **entire** stream,
with one initial zero state and no regional reset. Lean checked this
route-and-inner moment identity with only standard axioms. The outer
sampler and its expected measure still need formalization; this theorem
does not establish any finite distance margin.
`SpinCodes/Turbo/ColumnMixer.lean` proves that a uniform nonzero field
scalar sends each fixed nonzero column uniformly over nonzero values,
that independently sampled column maps have the corresponding product
law, and that they preserve the exact nonzero-column support. This is
the first part of the outer-measure argument. The remaining work is to
instantiate `F` as the 32-bit field, identify the binary adjoint action,
prove the eight-row RS[8,4] column-support enumerator, and then derive
the expected packet counts `C_v` under the actual local setup law.
(4) Formalize complete-measure domination and a theorem transferring
reference endpoints only under explicit numerical and route/inner premises.
(5) Use the old checker pattern: untrusted code emits rational/dyadic data;
Lean checks soundness from local operators through placement and powers,
then complete occupation coverage and the exact positive sum in bounded
modules. Hashes remain provenance, never the proof. (6) Connect the sum
to the sampled encoder, prove injectivity for every setup, and audit the
final theorem's recursive axioms. A first pilot should check one sparse
support and one `q>=3` transfer for K18, where the existing 22-state map
can be reused; this tests semantics before thousands of endpoints are
generated. Port the separate K16 15-state map after that interface works.

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
- 2026-11-07: Re-audited `DenseOccupationGeometryDefs.lean`,
  `DenseOccupationGeometry.lean`, `DenseOccupationGeometryCover.lean`, and
  `scripts/emit_cover.py` before designing the full cover. The existing
  two-coordinate `Rect`/`Tree` and `check_sound` match the revised domain.
  The generated global table/tree and assembly must be new; split-tree
  subclaims should be compiled independently, following the old generator's
  response to its 16 GB monolithic-check failure.
- 2026-11-07: The untrusted search now emits a 22.7 MB candidate containing
  the exact split tree, 46,317 rectangles, and 113 shared rational witnesses.
  Independent diagnostic scripts checked its tree bookkeeping and all
  proposed rational transfer, prefactor, and four-corner exponent bounds.
  There are 81,660 distinct witness/vertex checks after deduplication.
  These diagnostics do not establish Lean coverage.
- 2026-11-07: Reused the old `check_sound` interface in
  `ScalableDenseCoverCore.lean`. Lean checked a generated 32-leaf subtree and
  a 256-leaf geometry subtree. The 256-leaf geometry build took about 25 s,
  largely import time. A numerical module with 24 distinct vertices passed
  in about 55 s; one with 85 vertices hit kernel memory limits. Numerical
  chunks need to be much smaller than geometric chunks. Lowering the search
  margin from `1e-7` to `1e-8` reduced the box count only to 45,706 and is
  not a useful shortcut.
- 2026-11-07: A generated 8-leaf module checks 3 shared transfer witnesses,
  24 distinct vertices, all leaf point rates, and its root `DenseClaim` in
  about 63 s at a 4 GB Lean limit. A 12-leaf full claim failed at 4 GB but
  passed at 8 GB. A 32-leaf full claim with 85 distinct vertices passed at
  12 GB in about 116 s. All final claim audits list only standard axioms.
  The 32-leaf size is a provisional resource ceiling, not yet validated
  across regions with harder arithmetic. Generated bulk sources remain in
  ignored `.lake` storage until the complete proof and artifact layout are
  decided.
- 2026-11-07: Two separately checked 8-leaf numerical chunks compose across
  their actual split into a 16-leaf `DenseClaim`; the assembly theorem has
  only standard Lean axioms. The cover candidate has 39 top-level segments.
  A simple 32-leaf cut would produce 3,891 modules; 32-leaf full claims
  compile at a 12 GB Lean limit but may vary by region. This is a resource
  pilot, not a full-cover theorem.
- 2026-11-07: Re-audited the old cover's `LogTree`, `treeLog_oracle`, and
  generator before changing the revised numerical checker. The old checker
  stores `flogQ key scale n` values and proves every lookup sound; its
  single-integer key is insufficient for the revised KL arguments with
  varying denominators. The new cover's 81,660 distinct witness/vertex
  checks use only 1,291 distinct `(witness, α)` KL terms, 32,056 distinct
  `(witness, x)` KL terms, and 113 witness pairs `(ρ,z)`. These counts suggest
  a factorized checked vertex: prove each KL/log interval once per chunk,
  then check each vertex with cheap fixed-point arithmetic. Cross-chunk
  duplication remains; compare actual compile time and memory before
  considering a general rational-key log table.
- 2026-11-07: `ScalableDenseVertexCached.lean` now proves a factorized vertex
  check sound by reducing it to the existing checked vertex theorem; the
  resulting theorem has only standard axioms. The generator emits and Lean
  checks the distinct KL intervals and witness logarithms, then reuses them
  at each vertex. On the same 32-leaf, 85-vertex subtree, its complete root
  claim compiled in about 83 s at a 12 GB Lean limit, versus about 116 s for
  direct checks. Cached 64-leaf / 157-vertex and 128-leaf / 235-vertex chunks
  also passed at that limit, in about 125 s and 189 s respectively; the
  latter `.olean` is 13.1 MB. The 128-leaf check used roughly 4.5 GB working
  memory during its run. All three root claims have only standard axioms.
  These are build-resource measurements, not evidence that the entire cover
  has been checked. Before bulk generation, test a hard region of comparable
  size and choose a chunk cap by both distinct vertices and leaves.
- 2026-11-07: The 101-leaf region containing the hardest candidate box
  (index 16606) also passed the cached Lean checker at a 12 GB limit, in
  about 156 s, with four distinct witnesses and 192 vertices. Its root
  theorem has only standard axioms. The full candidate JSON is retained in
  ignored `.lake` storage and can be regenerated with
  `shared_inner_dense_certificate.py --max-boxes 100000 --precision 256
  --cover-output PATH`; it is proposed data, not a certificate.
- 2026-11-07: Revisited the old `checkTree`/`checkTree_sound`/`leafOK_sound`
  chain after the larger pilots. The old certificate encodes each chunk as
  a compact tree and data table, then proves one Boolean checker sound.
  The current generated pilots instead emit a theorem for every KL interval,
  vertex, and leaf. They validate the arithmetic and assembly, but copying
  that format across 46,317 leaves would create excessive proof source and
  repeated elaboration. Before bulk generation, define a compact revised
  leaf record and checker whose success implies `DenseTarget` over its box;
  combine it with the existing geometry `check_sound` and bounded chunk
  assembly. The checker must fail on a missing witness or support index and
  must check positivity, exact transfer, four vertices, and uniform prefactor.
- 2026-11-07: `ScalableDensePackedLeaf.lean` now has a compact, decidable
  line/witness/box checker and a soundness theorem that reaches the actual
  `DenseTarget` through `scalable_dense_checked_leaf`. Its indexed witness
  lookup returns `false` on a missing entry. The chunk checker combines
  numerical checks with the old rational-rectangle `check_sound`; no generated
  assertion of coverage is trusted. A one-leaf pilot and 8-, 32-, 128-leaf
  chunks compiled with only standard axioms. A single `decide` over the whole
  numerical list did not reduce, but each leaf check did; combining the
  independently checked leaf equalities with `List.all` fixed this without
  weakening the theorem. The 128-leaf compact source was about 94 KB versus
  about 621 KB for the earlier theorem-per-vertex pilot.
- 2026-11-07: A second compact checker verifies the line and each witness
  transfer once per chunk, then verifies only box conditions at each leaf.
  Its soundness theorem also has only standard axioms. The 128-leaf chunk
  passed at an 8 GB Lean limit in about 150 s, versus about 198 s when those
  checks were repeated per leaf. A 101-leaf region containing the hardest
  candidate box passed in about 133 s with four witnesses. This shared
  compact format is the current preferred bulk format; it still recomputes
  KL/log terms at shared vertices, so a rational-key log table remains a
  possible later optimization if bulk resources demand it.
- 2026-11-07: A 256-leaf, 21-witness compact chunk passed at an 8 GB Lean
  limit in about 238 s. A second chunk at the selected vertex cap (200 leaves,
  400 distinct vertices, 20 witnesses) passed in about 194 s. The 256-leaf
  and 400-vertex cut yields 611 chunks over the 46,317 candidate leaves.
  This budget has two distinct stress tests; it is still possible that a
  later chunk needs a smaller cut.
- 2026-11-07: `emit_scalable_dense_packed_cover.py` emitted a complete
  proposed source tree into ignored `.lake/fullcover`: 611 part files (about
  34 MB), 39 segment assembly modules, and a global `ScalableDenseRates`
  module. Its manifest records candidate SHA-256
  `164e3a944752a4a38003bfd5bca3917a3b69c79249375b678d3ead667de6f694`.
  All part bounds and split choices are checked in Lean when compiled. Parts
  17 and 18 passed independently, and their generated Segment 9 assembly
  passed with only standard axioms. A sequential, resumable build script
  compiled Part 0 and checked its axiom output. The manifest's
  `source_emitted` status does **not** claim a full certificate; most parts,
  segments, and the global theorem have not been compiled.
- 2026-11-07: Began the sequential checkpointed build of the full Scalable
  cover at an 8 GB Lean limit and below-normal process priority. The build
  has checked Parts 0--2 and continues; 611 parts, 39 segment assemblies,
  and the global rate theorem are still required. Separately audited the
  finite Turbo proof against the old Lean first-moment framework and the
  three RS8 transfer scripts. The length-dependent proof split and the
  difference between authenticated Python receipts and Lean-checked
  endpoints are recorded in the finite gate above. The retained receipts
  are present in the canonical checkout but not in this worktree.
- 2026-11-07: Added and checked the value-carrying Turbo packet route and
  its support-law and count bridges. Audited the K18/K20 endpoint producers:
  their outward physical operators are 17 `11×11` matrices originating in
  a `2^22` state census, and the independent Python audit authenticates
  rather than recomputes these inequalities. The finite Lean pilot must
  therefore prove a local operator bound before checking endpoint sums.
  The Scalable bulk build remains active and has reached at least Part 12.
- 2026-11-07: Added the actual four-bit packet serialization and uniform
  22-state sampled inner law. The full route-and-inner moment identity
  compiles with only standard axioms. This reuses the old step-kernel proof
  after auditing its input and state-continuity semantics; the next finite
  gate is the RS8 outer measure and its product setup with this law.
- 2026-11-07: Checked the generic nonzero field-scalar action and the
  product law across independently randomized columns. The output
  support agrees with the pre-randomization support for every setup.
  All new Turbo route, inner, and scalar-law declarations have only the
  standard Lean axioms; the full RS8 enumerator and finite endpoints
  remain open.
- 2026-11-07: Audited the old `Structured/BASpectrum.lean` outer argument
  before starting the Turbo RS8 outer. Its local-code injectivity and
  weight-enumerator interface are the relevant pattern; the Golay--BAA
  polynomial and its numerical bounds cannot be reused for an RS row.
  `Turbo/RSRow.lean`, `Turbo/RSRowCode.lean`, and `Turbo/RS8Code.lean` now
  prove that four distinct evaluations determine a degree-below-four row,
  that the row message space has cardinality `|F|^4`, and that eight such
  rows form `|F|^32` messages with at least five active columns for every
  nonzero group. All these declarations compile with only standard Lean
  axioms. The argument is currently generic over a field and an injective
  eight-point evaluation map. It has not yet identified the concrete GF16
  basis (modulus `0x13`) and evaluation points `0..7` used by the encoder.
  The support floor is only a first semantic gate: it does not provide the
  exact RS8 support enumerator needed by the finite first moment. Next prove
  the shortened-code cardinality from the four-position information-set
  property, derive the exact support counts by inclusion--exclusion, and
  connect those counts to the old first-moment framework. Keep the actual
  field representation and binary adjoint law as separate explicit gates.
- 2026-11-07: Closed the next RS8 algebraic gate without importing a
  numerical receipt. Evaluation at any four distinct points is now a
  bijection from an RS row to four symbols, and evaluation at any four
  distinct columns is a bijection from the eight-row group to four
  32-bit columns (`RSRowCode.lean`, `RS8Code.lean`).
  `RS8Shortening.lean` proves that zeroing any three selected columns
  leaves exactly `|F|^8` groups, of which `|F|^8-1` are nonzero; every
  nonzero survivor has exactly five active columns. At `|F|=16`, this is
  `2^32-1` words for each specified five-column support. These are the
  semantic ingredients of the `A_5` term, though the theorem summing over
  all 56 supports is not yet stated. The `A_6,A_7,A_8` terms still require
  shorter zero constraints and inclusion--exclusion. The three-zero
  theorem assumes an injective four-position extension; proving that
  every three-position selection admits one over the concrete eight
  positions, and connecting to the actual GF16 implementation, remain
  separate gates. All new declarations have only standard Lean axioms.
- 2026-11-07: Generalized the shortening pilot in
  `Turbo/RS8ShorteningGeneral.lean`. For any `m ≤ 4` distinct selected
  zero columns, the surviving group messages are in bijection with the
  `4-m` free 32-bit columns; their count is
  `|F|^(8*(4-m))`. This gives the exact inclusive counts needed for
  supports of sizes four through eight. The generic proof and its
  recursive axiom audit pass. It still stops short of the *exact*
  `A_5,...,A_8` spectrum: we must formalize extension of arbitrary
  zero-column sets to four distinct positions and use checked
  inclusion--exclusion to require the remaining columns to be nonzero.
  The present proof is the algebraic side of the old enumerator interface;
  the expected packet measure and first-moment bridge are not yet proved.
- 2026-11-07: Following the old `BASpectrum` weight-class bookkeeping,
  used Mathlib's checked inclusion--exclusion theorem to count RS8 words
  with **exactly** a specified set of active columns
  (`Turbo/RS8Spectrum.lean`). The intersection in each term is a
  zero-column constraint. `Turbo/RS8PositionSets.lean` builds a bijective
  ordering with any chosen zero set first and proves the shortening
  cardinality for arbitrary sets of size at most four.
  `Turbo/RS8SpectrumCount.lean` proves that larger zero sets leave only
  the zero group, identifies every intersection with the corresponding
  zero-set subtype, and combines the cases into an exact finite-set-sum
  formula for each support. No codeword enumeration or external receipt
  is trusted in this derivation. The new declarations compile with only
  standard Lean axioms. The paper's closed binomial `A_h` and its
  expected packet-weight polynomial still need to be derived from this
  checked formula; in particular, this formula is not yet the full
  sampled-setup first moment.
- 2026-11-07: `Turbo/RS8SpectrumNumeric.lean` now kernel-checks the
  simplification of the exact-support finite sum on all 256 subsets of
  the eight positions. The first attempted equality failed at the empty
  support; this is the unique zero message and contributes one, whereas
  the paper's `A_h` applies only at `h ≥ 5`. With that case explicit,
  Lean verifies the per-support alternating binomial sum, sums over
  `binom(8,h)` supports, and proves the result equals the actual count
  of eight-row messages of column weight `h`. The final `A_h` theorem
  has only standard Lean axioms and assumes `|F|=16` plus distinct
  evaluation points. The remaining outer gate is the *randomized*
  expected packet measure: model the 32-bit field coordinates and
  binary adjoint multiplication, then combine independent column laws
  with this deterministic message count.
- 2026-11-07: Audited the old `ScalableStepLaw.uniform_state_map_update`
  transitive-action proof before adding the column-coordinate bridge.
  `Turbo/ColumnCoordinates.lean` now proves that any zero-preserving
  bijection from the field to a finite coordinate type transports the
  uniform nonzero scalar output to a uniform nonzero coordinate word.
  This theorem has only standard Lean axioms. It is conditional on the
  coordinate bijection and treats forward multiplication; the actual
  GF(2^32) tower basis, its binary adjoint map in the transposed kernel,
  and the count of four-bit packet supports are still open. Do not use
  this conditional law as the finished sampled outer measure.
- Next action: run the checkpointed part build, adapting any resource-heavy
  chunk by a finer cut rather than weakening the target. Compile all 39
  segment claims and the global `ScalableDenseRates` declaration, audit its
  recursive axioms, then complete the actual sampled-setup theorem. The
  checkpointed build is currently running and had checked Parts 0--51 at
  the last inspection; do not restart it merely to check status. In
  parallel, prove the independent-column expected packet measure and pin the
  actual Turbo setup/event, and prove its expected outer measure. Use the
  finite audit above before choosing numerical checker data structures.

At each future checkpoint, update the status and next action here before
reporting completion or changing strategy. Preserve failed candidates and
their reason in research notes when they inform a later decision.
