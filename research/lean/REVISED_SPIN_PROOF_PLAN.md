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

Treat this file as the decision log for the active proof. At each completed
lemma, failed strategy, or change of interface, update the relevant gate and
the next action. Before starting a new component, identify its old-proof
counterpart by declaration and check its assumptions and conclusion in the
source. If no counterpart exists, record that fact and the new obligation.
The audit is a design input, not a requirement to copy the old code.

## Target and present status

The Scalable target is the paper's one-bit route, shared growing Golay--BAA
outer, and independent uniform `GL(22,2)` updates with the actual 22-bit state
map. The desired result is a rate-half family whose probability of minimum
distance at most `floor(0.11*N)` tends to zero under its original sampled
setup. The finite Turbo target is the separate four-bit packet construction
at the lengths and setup-failure margins claimed in the paper. Do not transfer
the old BCH/IMT theorem or the Scalable theorem to Turbo by name similarity.

For Scalable SPIN, the native-length distance, exact-rate, schedule, and
ordinary/transposed linear circuit-size claims are now
checked in `SpinCodes/UniformMixer/ScalableCertified.lean`. The complete
generated cover proves
`ScalableDenseFullCover.rates : ScalableDenseRates (1/10000000) 10000000`.
It discharges the numerical hypothesis in the fixed/sparse/dense assembly;
the original-setup theorem charges bad outer selection once. The final
theorem proves that the actual realized code's probability of minimum
binary Hamming distance at most `floor(0.11*Nsched m)` tends to zero under
`scalableSetupLaw`. Every realized encoder is linear and injective, has
exactly `Nsched m` output coordinates, dimension `Nsched m/2`, and rate 1/2.
The final aggregate compiled successfully and all eight printed declarations
depend only on `propext`, `Classical.choice`, and `Quot.sound`.

The actual encoder and its adjoint have explicit coefficient-AND/XOR circuits
using at most `272*Nsched m` gates. The constant is a conservative formal bound,
not an implementation timing estimate. This covers fixed sampled setup and
free permutation wiring; setup generation and C++ refinement are not covered.
The schedule additionally satisfies the eventual rounding bound
`0 <= bsched m - (39/4)*log2(Nsched m) < 24`.
The user accepts natural lengths and has removed the
requested-length wrapper from the current scope. A clean dependency
replay remains a delivery check; this pass reused existing checked objects
and the completed generated cover. Finite Turbo remains open as described
below. The earlier status entries are retained as historical checkpoints.

For Turbo SPIN, Lean checks the RS8 support spectrum, exact expected
packet measure, complete value law after the sampled row shuffle, and
the actual 22-state setup. The vector encoder is connected to that setup's
weight experiment and is proved binary linear, injective, and exactly
rate 1/2 at every positive natural length in this 22-state model.
Its first moment is now partitioned by active
groups and expressed exactly as a sum of packet-array counts multiplied
by one joint regional-route and continuous-inner failure probability.
The symbolic bridge currently covers the 22-state K18/K20 construction.
The exact 16-slot physical input law, full-state weighted step kernel, and
zero-inclusive selected-slot thinning identity are also checked. The last
identity justifies the producer's `15/16` mixture at the complete-word level;
authentication of its 11-coordinate numerical entries remains open. The revised
RS8 active-group expected count is now pointwise bounded by the retained
uniform-reference atom factor, before and after the sampled row shuffle.
The exact reference route is now disintegrated into hypergeometric
physical-step count profiles, with complete four-bit values retained.
The checked full-state kernel and lower-tail inequality bound the actual
shared-outer active-set first moment by a chronological product of 64
averaged reference-region kernels. This is a semantic bound, not yet a
finite numeric certificate.
An 11-coordinate state-measure cone, a complete symbolic selected-slot
operator, and preservation of complete-row soundness under exact binomial
thinning are checked. One complete retained rational local family is now
authenticated for the exact-nonzero physical kernel and its selected-step
mixture; the other required families remain to be checked. Conditional probability
capping, actual S-specific routing fibers, conditional reference domination,
and chronological paired-matrix fractional subadditivity are checked.
The conditional reference mean is now connected to the full-state product
of selected-step kernels. Exact physical Cauchy rows, coarse birth rows,
and chronological eleven-coordinate cone propagation are also checked,
with their numerical hypotheses explicit. The paired hypergeometric
coefficient identity, its ordered matrix DP, and the positive-marker bound
are now checked. The fractional tail assembly retains the conditioning and
the outer multiplicities in their correct order; an unconditional reference
mean alone would not justify it.

The actual expansion map's 25-profile census and the integer interval
coefficient/aggregate checkers are checked. The exact one- and two-packet
birth histograms are authenticated against the actual input layers, and
their normalized birth cones (including the empty input) are checked.
The exact packet character formula, zero-feedback character average, and
the class-maximum/uniform-mean weighted adapter are also checked.
The polynomial-to-actual-moment joins and sparse regional product joins
are now checked. The sparse analysis also permits choosing a Chernoff
parameter separately for each support before summing over supports.
The first complete 17-row rational profile family passed its integer/Fix
kernel checks and actual physical/selected-row soundness assembly in
`PhysicalProfileFamilySmoke22`. All four exported soundness/nonnegativity
theorems have only the standard axioms, with no remaining numerical or
birth premises. This is one local family, not a finite distance certificate.
Only successful complete-module replays count as delivery.

No end-to-end Turbo finite certificate is closed yet. Complete external
diagnostics retain 10% distance and total margins of 58.2408311085 bits at
K18 and 51.8224664170 bits at K20. These diagnostics do not discharge the
Lean numerical premises. The next gates are authenticated physical-row
data, the remaining original-setup sparse composition, and exact numeric
checks covering every required shell and occupation, followed by the final
probability sum. Finish K18 at the paper's rounded 58.24 bits before K20 at
51.82 bits. The separate K16 inner and a C++ refinement proof are outside
this goal.

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

### Sampled column adjoint and physical input law

The immediate task is to justify the packet law for the *transposed* encoder.
For a fixed nonzero 32-bit column `x`, the paper samples a uniform nonzero
field coefficient `c` and applies the binary adjoint of multiplication by
`c`. We need a checked theorem that this image is uniform over all nonzero
32-bit columns. The forward multiplication theorem in `Turbo/ColumnMixer.lean`
does not establish this: its output has the wrong orientation.

Old-proof audit: `UniformMixer/ScalableStepLaw.uniform_state_map_update`
turns a uniform sampled `GL(22,2)` state map into an exact transition law by
counting the event fiber; `uniform_state_map_expect` then transports the law
to observables. `Structured/ConcretePairing.chi_adjoint_of_basis` proves a
concrete binary transpose identity from basis incidences. The first supplies
the probability-space and observable pattern; the second supplies the
implementation-identity pattern. Neither proves the Turbo field adjoint
action, so the new algebraic bridge is necessary.

First prove a field-generic lemma under an explicit nondegenerate binary
pairing and an additive field-to-column equivalence. Its adjoint identity
must imply that, for fixed nonzero `x`, the map from nonzero coefficients to
nonzero adjoint images is a bijection: a zero image would annihilate the
pairing against every field element, and equal images for two coefficients
would do the same for their nonzero difference. Transport the uniform law
through that bijection and then through the product sampler for independent
columns. This abstract lemma is an intermediate theorem. A separate checked
bridge must show that the concrete GF(2^32) multiplication, binary basis,
and transposed kernel satisfy its pairing identity and nondegeneracy. Only
then may the local packet generating function be attached to the actual
sampled Turbo setup.

**Algebraic gate closed:** `Turbo/ColumnAdjoint.lean` proves the nonzero
image, coefficient injectivity, resulting nonzero-column equivalence,
uniform image law, and product law for independently sampled coefficients.
Its recursive axiom audit reports only `propext`, `Classical.choice`, and
`Quot.sound`. Its hypotheses include a nondegenerate pairing, an additive
coordinate equivalence, and the adjoint identity. These hypotheses are not
yet proved for `spin/src/rs8/Field.h` and `SetupScalar.cpp`; the concrete
field-tower and transpose bridge is the next semantic gate. The old
`ConcretePairing` basis-incidence proof is the model for that check. The C++
scalar tests are implementation checks, not a Lean proof of this bridge.

**Ideal GF(2^32) gate closed:** `Turbo/ColumnAdjointField.lean` instantiates
the abstract theorem with Mathlib's `GaloisField 2 32` and a binary linear
basis. It defines multiplication as a linear map, takes the actual matrix
transpose as its adjoint, and proves the dot-product pairing identity using
`Matrix.dotProduct_transpose_mulVec`. This is the precise counterpart of the
old `ConcretePairing` incidence argument. It then regroups the 32
coordinates into eight four-bit packets and proves both the sampled
one-column generating function and its product over `h` independent active
columns. These declarations have only the standard Lean axioms. The chosen
basis is abstract; the probability law is basis-independent. The file does
not verify that the C++ tower arithmetic in `Field.h` realizes a field or
that its generated reverse columns implement this particular matrix
transpose. That implementation bridge remains separate from the paper's
ideal sampled-law proof and must not be conflated with it.

After this gate, combine the checked `A_h` count with the independent-column
packet generating function to derive the full expected group measure. Then
prove packet-route support symmetry and connect that measure to the existing
route/inner moment law. The finite numerical endpoint checker remains a
separate gate; no present local lemma implies a distance margin.

**Expected group measure closed:** `Turbo/RS8FieldMoment.lean` takes eight
actual RS[8,4] polynomial rows over `GaloisField 2 4`, evaluates their
eight output columns, identifies each column with 32 binary coordinates,
and samples eight independent nonzero `GaloisField 2 32` coefficients.
`rs8_fixed_group_packet_gf` proves the exact moment for each fixed group,
including zero columns. `rs8_nonzero_group_packet_gf_closed` sums over
every nonzero group message and rewrites the result using the checked
`A_h` MDS formula. This follows the old `Framework.expect_Z_eq` and
`FiniteProductLaw.prod_sum_eq_sum_prod` separation: one sampled setup is
shared by all messages; independence is invoked only among its eight
column coefficients, and linearity of expectation handles the message
sum. The theorem assumes eight distinct evaluation points. Its arbitrary
binary coordinate basis is sufficient for the ideal distributional law,
whose value is basis-independent; it has not yet checked the exact C++
GF16 polynomial basis and point encoding. Recursive axioms are the three
standard Lean axioms. `rs8_sampled_packet_weight_range` now proves that every
nonzero message has between 5 and 64 active packets under every sampled
nonzero coefficient tuple, so the generating function has exactly the
paper's claimed degree range. This uses the old support floor and the new
adjoint nonzero theorem. `rs8ExpectedPacketCount` now defines each `C_v`
as the expected number of nonzero messages of packet weight `v`.
`rs8_expected_packet_count_gf` proves the paper's coefficientwise
generating-function equation, and
`rs8_expected_packet_count_total_numeric` checks that their sum is
`2^128-1`. These use the old `Prob.expect_comp` fiber decomposition,
with the every-setup `5..64` range as its coverage premise. No routed
support/label law or finite distance certificate follows yet. The next
gate is the value-carrying packet shuffle and its placement law, then
composition with the inner moment.

### Next gate: value-carrying packet shuffle

The old `Structured/ConcreteShufflePermutation.lean` proves
`uniform_permutation_shuffleLaw` by first identifying a Boolean support
fiber and then proving `shuffle_fiber_total`: for fixed source support
`S` and target support `T` of the same size, every `T` has the same
number of uniform permutations. `Turbo/PacketRoute.lean` already proves
that routing packet values induces this same support action. What the old
proof does **not** establish is the law of the nonzero four-bit labels.

Extend the old fiber count rather than replacing it. For a fixed source
support `S`, choose a uniform labeling of its active positions, then a
uniform packet permutation. For a fixed target packet word `y`, each
permutation sending `S` to `support(y)` determines exactly one source
labeling, namely the pullback of `y`. The resulting fiber is therefore
equivalent to the old support-permutation fiber. Count source labelings
as `15^|S|` and use `shuffle_fiber_total` to obtain a uniform target
law among all 64-packet words of weight `|S|`. First prove the
deterministic support identity and the fiber equivalence; then the
finite probability formula. Finally connect this law to the sampled
RS8 outer: independent column adjoints make labels uniform conditional
on a fixed pre-shuffle support, and summing these support classes with
the checked `C_v` gives the required local atom probabilities. Keep the
actual group and regional permutations sampled independently as in
`ConcreteRoute.seedLaw`. The support-only theorem and a labeled-word
theorem have different scopes and must remain distinct.

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

### Audit before the Turbo physical-operator pilot

The old `Structured/ConcreteFixedRegionComposition.lean` proves that
chronological state kernels compose as matrix products and that
independently shuffled regions can be averaged one at a time. Its state
and input are binary support sets, so its *matrix entries* do not apply
to Turbo's four-bit packet labels. Reuse its composition pattern and
finite-product argument after proving a new packet-valued local operator.
`Structured/ConcreteCounts.lean` identifies frozen tables with the
actual expansion/feedback maps through `MapSpectrumBridge` and
`FiberNumericsAll`; `FiberCertificate.lean` proves an integer quadratic
cap check. They illustrate checked data attached to explicit maps, but
neither proves the Turbo GL(22,2) operator or its 11-coordinate quotient.

The retained K18 producer is
`permutation_locality/packed_mixing/s16_closure/sparse_kernel.py` plus
`kernel_birth_density.py`, invoked by `k18_margin60/tune.py` and then
`replay_fractional.py`. Its 17 outward physical matrices are `11×11`:
coordinate 0 is zero-state mass, coordinates 1 and 2 are arbitrary and
uniform-dominated nonzero mass, and the remaining eight coordinates
retain expansion-weight classes. `sparse_kernel.outward_at_z` uses
weighted first and second moments, a Cauchy bound on returning to zero,
and uniform-GL refresh of any nonzero state. The optional capped-density
row uses an exact Walsh inversion of the weighted feedback law, a
uniform-density cap, and nonnegative residual mass by expansion class.
Its floating continuation potential only *selects one complete candidate
row*; taking entrywise minima is invalid. `replay_fractional.py` thins
for zero packet labels, composes two physical steps chronologically,
then takes outward fractional powers and regional products. The receipt
verifier authenticates the resulting numbers and their union but does
not prove any matrix majorizes the actual encoder kernel.

The next proof interface must define a sound representation cone from
exact nonnegative measures on `State22` to these eleven coordinates,
prove that each local matrix preserves this representation for one
physical step with `j=0..16` active packet slots, and identify the
weight-class and feedback records with the explicit 22-state map.
`UniformMixer/ScalableStepLaw.uniform_state_map_update` proves the
uniform-GL refresh law and is reusable; its Bernoulli-bit input spectrum
does not prove the packet/feedback moments. For a first K18 pilot, check
one complete local occupancy matrix and its effect on an entering
state representation, including the capped newborn row only if the
baseline operator is too weak. Then check one sparse endpoint and its
placement average against `turboRemainingFailure22`. Use rational or
dyadic outward intervals for the `exp(-tilt)` weight. The old checked
numeric-table design (small sound checker plus generated chunks) is the
model for scaling the 17 operators and all endpoint data. A source hash
is provenance and cannot replace the local operator inequality.
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
- 2026-11-07: `Turbo/PacketSupportCount.lean` now proves that a specified
  support of `v` four-bit packets inside a 32-bit column has exactly
  `15^v` possible labels. Summing over supports gives
  `binom(8,v)·15^v` binary columns, and under the uniform nonzero-column
  law the packet-weight probability for `v ≥ 1` is exactly
  `binom(8,v)·15^v/(2^32-1)`. The zero column is handled separately.
  These theorems compile with only standard Lean axioms. A direct
  specialization of the nonzero-cardinality proof to the concrete
  32-bit function type exhausted the 6 GB Lean kernel limit; proving
  that cardinal identity for an abstract finite type first avoids
  expanding the `2^32`-element fintype. This mirrors the old proof's
  practice of keeping symbolic counting separate from checked finite
  arithmetic. Next derive the local generating polynomial and its
  independent product across active RS columns. The field-coordinate
  and binary-adjoint bridges remain open.
- 2026-11-07: Reused the old `FiniteProductLaw.prod_sum_eq_sum_prod`
  identity in `Turbo/PacketProductMoment.lean` to prove a general
  product-expectation theorem for independently sampled columns.
  `Turbo/PacketLocalMoment.lean` now proves the exact local generating
  function `((1+15z)^8-1)/(2^32-1)` under a uniform nonzero 32-bit
  column and its `h`th power for `h` independent active columns. These
  moment theorems compile with only standard Lean axioms. This closes
  the probability calculation **after** active columns have been
  identified. The sampled RS8 outer still needs a proof that its
  independent field multipliers (including the transposed adjoints)
  induce this product law, followed by summing these moments over the
  checked `A_h` classes. No finite distance endpoint follows yet.
- Next action: run the checkpointed part build, adapting any resource-heavy
  chunk by a finer cut rather than weakening the target. Compile all 39
  segment claims and the global `ScalableDenseRates` declaration, audit its
  recursive axioms, then complete the actual sampled-setup theorem. The
  resumed checkpointed build had checked Parts 0--65 at the last
  inspection; do not restart it merely to check status. In
  parallel, prove the independent-column expected packet measure and pin the
  actual Turbo setup/event, and prove its expected outer measure. Use the
  finite audit above before choosing numerical checker data structures.
- 2026-11-07: The first full-cover process stopped after checking Part65;
  its PID was absent, no Lean process was running, and the output/error
  files had not changed for hours. All 66 part `.olean`s remained present.
  Resumed the same manifest with the checkpoint-aware build script (PID
  39240, `research/lean/.lake/fullcover/resume.out` and `resume.err`).
  Inspect this handle and these files before any further restart. This is
  still a partial cover; no `ScalableDenseRates` theorem has been claimed.
  The next Turbo action at this checkpoint was the ideal GF(2^32) adjoint
  bridge; the subsequent entry records its result.
- 2026-11-07: The ideal GF(2^32) matrix-adjoint bridge and exact sampled
  packet generating functions now compile in `Turbo/ColumnAdjointField.lean`.
  The next proof step is to identify the eight RS8 output columns with
  `Bits32` under the eight-row coordinate transpose, then sum the
  independently sampled `h`-column law over the checked `A_h` classes.
  Only after that should the packet-route and inner first moments be
  combined. At last inspection the resumed Scalable cover process was
  live and had checked Parts 0--72; its error file was empty.
- 2026-11-07: `Turbo/RS8FieldMoment.lean` now closes the ideal sampled
  RS8 expected packet polynomial through the exact `A_h` classes.
  The old `Framework` and `FiniteProductLaw` interfaces were inspected
  before the group-sum proof. No coefficient receipt or independent
  per-message setup assumption was used. At last inspection the resumed
  full-cover process remained live through Part82 with an empty error
  file. The next finite gate is coefficientwise `C_v` and its uniform
  support/label placement after packet shuffling; the numerical finite
  operator checks remain open.
- 2026-11-07: Added an every-setup packet-weight range proof for each
  nonzero RS8 group message: 5 through 64 packets. It depends on the
  checked MDS support floor and the proved nonzero action of each sampled
  field adjoint. The full-cover process remained live through Part86 with
  an empty error file. The next finite step remains coefficientwise
  extraction of the expected `C_v` and the value-carrying shuffle law.
- 2026-11-07: Defined `C_v` as actual expected message counts and checked
  both its generating-function identity and `Σ_{v=5}^{64} C_v=2^128-1`
  in `Turbo/RS8FieldMoment.lean`. This reuses the old finite-valued
  expectation decomposition instead of trusting coefficients emitted by
  Python. The Scalable cover process was live through Part87 with an
  empty error file. Next audit the old uniform-shuffle fiber count and
  extend it to complete four-bit packet values before deriving Turbo's
  routed local measure.
- 2026-11-07: `Turbo/PacketShuffleValues.lean` now proves the exact
  value-carrying shuffle law. A labeled source support has `15^|S|`
  possible words, and `packetShuffleFiberEquiv` identifies the complete-word
  shuffle fiber with the support-permutation fiber of the old
  `Structured/ConcreteShufflePermutation.lean`. Its `shuffle_fiber_total`
  yields probability `1/(binom(64,|S|) 15^|S|)` for every target word of
  packet weight `|S|`, zero otherwise. This is an admission-free theorem
  whose recursive axioms are only Lean's standard three. It assumes a
  uniform labeling conditional on source support. Next derive that
  conditional labeling law from the independently sampled RS8 column
  adjoints, then compose the row shuffle, transpose, and regional shuffle
  with the actual inner event. The Scalable cover build was live through
  Part102 with an empty error log at this checkpoint.
- 2026-11-07: Extracted `field32_uniform_packet_column` from the existing
  sampled-adjoint moment proof. It proves an exact **value-level**
  pushforward: a uniform nonzero GF(2^32) multiplier acting on any fixed
  nonzero 32-bit column gives a uniform nonzero eight-packet column. The
  old support-only route had no counterpart for this coefficient-to-label
  statement; it is the missing hypothesis for the new value-carrying
  shuffle. Next prove the uniform labeling conditional on each nonempty
  pre-shuffle packet support, then use the independent column product to
  state it for one RS8 group. The subsequent theorem must retain shared
  setup across messages and independent row/regional seeds. The Scalable
  cover was live through Part103 with an empty error log at inspection.
- 2026-11-07: `Turbo/PacketLabelConditioning.lean` now checks the missing
  local conditional law: a uniform nonzero eight-packet column has support
  `s` with probability `15^|s|/(2^32-1)`, and, for nonempty `s`, every
  compatible value has conditional probability `15^(-|s|)`. The support
  fiber is proved equivalent to the previously counted unrestricted fiber;
  the new field-adjoint pushforward supplies the actual sampled source.
  Recursive axioms are only Lean's standard three. Next combine the
  value-level field pushforward and this conditioning theorem for the eight
  independent RS8 columns, giving a group-level uniform-label law. Then
  compose with the complete-word packet shuffle theorem; keep the
  coefficient seed shared between different messages in the first-moment
  count.
- 2026-11-07: Reused the old `FiniteProductLaw.piPMF_map` (the same product
  pushforward interface used by `ColumnAdjoint.uniform_active_adjoint_columns`)
  to prove `field32_uniform_active_packet_columns`: independently sampled
  coefficients on any fixed family of nonzero columns yield independent
  uniform nonzero packet columns. The theorem compiles without admissions
  and has only the standard three Lean axioms. The remaining RS8 group step
  must restrict the eight coefficient coordinates to the actual active
  column support, show zero columns stay zero, and then combine these
  independent column laws with the conditional label theorem. This
  restriction/product bridge is new for the Turbo RS8 outer; inspect the
  old finite product and route pushforwards before implementing it.
- 2026-11-07: `Turbo/SampledRS8PacketLaw.lean` now gives the complete
  eight-column law for an actual fixed RS8 group message and one sampled
  coefficient tuple. It reuses `FiniteProductLaw.piPMF_map` exactly as
  the old route and outer product proofs do. The point formula says each
  inactive RS column is zero with probability one and each active column
  outputs any specified nonzero 32-bit packet value with probability
  `1/(2^32-1)`; it also proves equal mass for any two outputs with the
  same packet support in every column. This is stronger than the earlier
  moment calculation, is admission-free, and has only the standard Lean
  axioms. The next gate is to flatten columns to one 64-packet row,
  transfer support and weight through that exact bijection, and combine
  this value-level fiber symmetry with the old uniform row permutation
  count. The group coefficient tuple remains shared across all messages.
  This establishes the paper's single-group placement law only after the
  row-shuffle composition itself has been checked.
- 2026-11-07: `Turbo/SampledRS8Row.lean` now provides the exact
  eight-columns-to-64-packets equivalence. It proves that the sampled
  row's packet support has the `rs8SampledPacketWeight` used by the
  previously checked `C_v`, and transfers equal mass within every
  fixed packet-support fiber to the 64-packet row. The old
  `Structured/BASpectrum.wtF_concatBlocks` supplied the appropriate
  `Fintype.sum_equiv` pattern for the weight proof; the group law is new.
  No independent per-message setup has been introduced. Next prove the
  independent row-shuffle pushforward as a mixture over packet weights:
  after shuffling, each complete target word of weight `v` has mass
  `Pr[source weight=v]/(binom(64,v) 15^v)`. The already-proved
  `packetShuffleLaw_apply` is the fixed-support case. Connect that exact
  placement law to the `C_v` expected message counts, then add regional
  routing and the inner event.
- 2026-11-07: Closed that single-group placement bridge in
  `Turbo/PacketShuffleMixture.lean` and
  `Turbo/RS8ShuffledPlacement.lean`. The former proves a generic exact
  disintegration of any finite packet-word law uniform on fixed support
  fibers, using the old `FiniteLaw.map_prob`/`Prob.prob` interface; it
  composes this law with an independent row permutation and obtains
  `Pr[source packet weight=v]/(binom(64,v) 15^v)` for every complete target
  word of weight `v`. The latter pins the source weight event to the
  previously proved `rs8SampledPacketWeight` and sums over nonzero RS8
  messages, yielding exactly `C_v/(binom(64,v) 15^v)` as an expected
  complete-word count. All stated theorems compile without admissions and
  have only the three standard Lean axioms. This proves the paper's
  one-group value placement **before regional shuffling and the inner**;
  it does not yet certify a distance endpoint. Next audit the old
  `ConcreteRoutePermutation` and `ConcreteFixedRegionComposition`
  product/transpose steps, then lift the actual value-carrying route and
  inner event for one and two active groups. The finite numerical
  operators and complete occupation coverage remain open. The Scalable
  dense full-cover build remains an independent checkpointed process.
- 2026-11-07: Audited the next old-proof component before creating the
  Turbo setup. `Structured/ConcreteRoutePermutation.permutation_law_eq`
  uses independent row and region permutation products and matches the
  support projection of `Turbo/PacketRoute.routePackets_support_law`.
  `Structured/ConcreteFixedRegionComposition` factors region endpoint
  matrices under an independent product law; its binary support-only
  kernels cannot simply be reused for four-bit packet values. The new
  `Turbo/PacketInner22.routeInnerMoment22` already gives the exact
  value-carrying route and shared-state kernel for the 22-state inner,
  without a per-region reset. `Framework.Setup.expect_Z_eq` and
  `FiniteMoment.prob_weight_le` provide the correct actual-setup first
  moment and lower-tail interfaces.
- 2026-11-07: `Turbo/ActualSetup22.lean` now instantiates that framework
  with one sampled RS8 multiplier tuple per group shared across all
  messages, independently sampled route permutations, and independent
  GL(22,2) state maps. It checks the exact expected bad-message count,
  Markov step, and a lower-tail bound by the actual route/inner kernel
  moment. Importing this module exposed two independent declarations of
  `Packet4`; `PacketInner22` now imports the existing packet type instead.
  Both modules compile admission-free with only standard Lean axioms.
  This establishes the 22-state probability-space bridge, not the
  numerical upper bound or finite distance endpoint. The 15-state
  K=2^16 packet inner has no corresponding Lean module yet and needs a
  separate actual-setup connection. Next check the finite checker against
  the 22-state value kernel and build the 15-state analog only after
  auditing the shared inner proof infrastructure and implementation map.

- 2026-11-07: Before the Turbo occupation component, audited the old
  `Framework.Setup.Z_eq_sum_ZQ` and `expect_Z_eq_sum_expect_ZQ`,
  `FiniteMoment.expect_ZQ_eq`, and `Structured/ConcreteOuter.occupation_pos`
  / `occupation_le`. The generic first-moment partition is reusable;
  the old local-message and Bernoulli support counting are not. The
  Turbo occupation must count nonzero RS8 group inputs, prove it lies in
  `1..L` for every nonzero full message, and instantiate the actual
  setup's `ZQ` identity before any numerical occupation sum is claimed.
- 2026-11-07: `Turbo/Occupation22.lean` now proves that nonzero Turbo
  messages have occupation in `1..R*16`, partitions the actual bad-message
  count and its expectation by occupation, identifies each occupation
  term with the shared-outer-setup first moment, and partitions that term
  by its exact active group set. The proof reuses the old generic `ZQ`
  and `expect_ZQ_eq` statements; no finite numerical bound is inferred.
  All printed theorems have only Lean's standard axioms. Before a sparse
  fold, audited `Structured/ConcreteOuterCountingProduct.tupleCount_eq_prod`
  and `ConcreteOuterCountingProbability.sum_occupation_eq_activeSets`:
  the product factorization is valid only after fixing one active set and
  summing local messages under independent *group* setup coordinates.
  The old Bernoulli support measure cannot replace Turbo's complete
  packet values. Next prove the marginal of `turboOuterLaw22` for a fixed
  active set and push it through the complete `sampledRS8RowLaw`; then
  incorporate the common regional route and inner law for the `U₁/U₂`
  folds. Keep the local physical operator soundness gate separate.
- 2026-11-07: `Turbo/GroupArrayLaw22.lean` closes the first complete-value
  product bridge. For each fixed full message, the actual independent
  outer setup maps to the product of its sampled RS8 row laws. Summing
  over all messages with a *fixed* active group set gives an exact
  product of local expected complete-row counts. This matches the old
  `tupleCount_eq_prod` algebra while using four-bit packet values and
  the true RS8 setup. It makes no independence claim across different
  messages sharing a seed. The proof compiles with standard Lean axioms.
  The next sparse gate is to incorporate the *actual independent row
  permutations* in this product law, so the previously proved
  `rs8_expected_shuffled_word_count` supplies `C_v` at each active
  group. The regional permutations and one continuous inner state then
  enter one common nonnegative test function; do not multiply their
  event probabilities across groups or regions.
- 2026-11-07: `Turbo/RowShuffledArray22.lean` now incorporates the *actual*
  independent row permutations. The fixed-message pushforward is the
  product of the independently shuffled RS8 complete-word laws; its
  active-set expected count again factors by group. The active local
  factor is exactly `C_v/(choose(64,v)*15^v)` by the earlier
  `rs8_expected_shuffled_word_count` theorem. The generic pairing lemma
  proves that products of finite setup arrays can be paired without
  replacing the shared setup by per-message randomness. Lean checks
  these facts with only standard axioms. This closes the outer and row
  shuffle part of the sparse fold. Next connect the row-shuffled array
  law to the remaining sampled *regional* permutation and continuous
  22-state inner event, keeping that latter event as a single nonnegative
  test function. Then derive exact `U₁/U₂` expectation formulas before
  importing any numerical endpoint bounds.

- 2026-11-07: `Turbo/RemainingRoute22.lean` connects the row-shuffled
  complete-value counts to the actual `routeInnerLaw22`. It proves the
  route is exactly row shuffle followed by regional shuffle, disintegrates
  `qd` over row permutations, and keeps the regional permutations plus
  all 22-state map updates in one `turboRemainingFailure22` event. The
  resulting `turboSetup22_occupation_complete_formula` is an exact
  identity under the sampled setup for every positive occupation; its
  local factors include the proved `C_v` law. All printed theorems compile
  with only standard Lean axioms. This closes the symbolic first-moment
  bridge for K18/K20. It does **not** establish a finite distance bound:
  no checked local operator majorizes `turboRemainingFailure22`, no
  one/two-group numerical endpoints or full larger-occupation transfer
  have been proved in Lean, and K16 still needs its separate 15-state
  inner. Next audit the old `ConcreteFixedRegionComposition` and checked
  `ConcreteCounts`/`FiberCertificate` numerical soundness in detail,
  together with the retained Turbo operator producer, before selecting
  the first K18 sparse operator pilot. Check whether the producer's
  11-state quotient retains packet labels and chronological state; if
  not, fix its mathematical interface before importing any endpoint.

- 2026-11-07: Before a new Turbo physical-step component, inspected the
  completed `UniformMixer/ScalableStepLaw.uniform_state_map_update` and
  `RegionalEncoder.stepKernel`, the old structured region matrix product,
  and the retained `sparse_kernel.py` physical producer. The uniform
  GL(22,2) transition for each fixed input and entering state is already
  proved in Lean. The missing component is the **packet-valued local input
  law**: exactly `j` of sixteen four-bit slots are active, with uniform
  support and independent uniform nonzero labels, followed by one sampled
  GL map. Define that finite law and its exact weighted full-state kernel
  first. Only then introduce an 11-coordinate cone and a checked outward
  operator. This separates a semantic identity from the numerical upper
  bound and makes the old `stepKernel` reusable without misusing its
  Bernoulli-bit input distribution.
- 2026-11-07: `Turbo/PhysicalStep22.lean` now defines the exact local
  input law for `j=0..16` active four-bit packet slots. A checked
  support-fiber equivalence proves `15^j` labelings for each `j`-subset;
  the full layer has `choose(16,j)*15^j` words, and its point mass is
  exactly the reciprocal of that count. The 16 slots are definitionally
  the actual `PacketInner22.stepInput` slots. The new full-state
  weighted kernel equals sampling a fresh uniform `GL(22,2)` map for
  that step, is nonnegative when `z≥0`, and sums to one at `z=1`.
  Lean checked the printed theorems with only standard axioms. This
  fixes the semantic target for the retained 17 physical operators;
  it is not yet an 11-coordinate upper bound or a proof that the
  regional routing conditioned on slot counts has this product law.
  Next prove that conditional route law and define the sound
  11-coordinate measure representation, reusing the already proved
  `QuadraticSpectrumSupport` weight classes and the old finite checker
  only after matching their exact map and law.

- 2026-11-07: Before modeling the producer's zero-label thinning,
  rechecked `Turbo/PacketLabelConditioning.lean`,
  `Structured/ConcreteShufflePermutation.uniform_permutation_shuffleLaw`,
  and `k18_margin60/replay_fractional.py::thin`. The old support theorem
  establishes the uniform selected-position law, and the Turbo label
  theorem establishes uniform nonzero labels after conditioning on
  support. The producer's `15/16` thinning is another exact law: start
  with `j` selected slots whose four-bit values are uniform *including*
  zero, then condition on the number `k` of nonzero values. The selected
  set and packet values must remain one finite experiment. This audit
  established the complete-word mass and binomial-mixture obligations
  closed in the following checkpoint; a numerical `thin` matrix cannot
  stand in for that identity.

- 2026-11-07: `Turbo/PhysicalThinning22.lean` now proves that experiment.
  `selectedStepLaw22 j` samples a uniform `j`-subset of sixteen slots and
  independent uniform four-bit values, including zero, on the selected
  slots. For a complete word `x` with `k` nonzero slots, Lean proves its
  point mass is
  `choose(16-k,j-k)/(choose(16,j)*16^j)` when `k≤j`, and zero otherwise.
  `selectedStepLaw22_mixture` then proves the equality of complete-word
  laws with the sum of `choose(j,k)*15^k/16^j` times the exact-occupancy
  laws from `PhysicalStep22.lean`. `turboSelectedKernel22_mixture` lifts
  that equality to every entry of the full `2^22`-state weighted kernel.
  This matches the producer's `thin` weights, including their
  zero-label semantics. The recursive axiom audit reports only Lean's
  standard axioms.
  **Old-proof correspondence:** `Structured/ConcretePlacementLaw` supplies
  the uniform selected-set inclusion count; `Turbo/PacketLabelConditioning`
  supplies the support/label separation pattern. Neither old theorem
  models selected values that may vanish, so this finite experiment and
  mixture proof are new.
- 2026-11-07: Rechecked the *interface* of this thinning identity against
  `Turbo/RowShuffledArray22.lean` and the retained K18
  `SLACK_REVIEW.md`. The actual RS8 group law has uniform **nonzero**
  packet labels conditional on its support; it is not itself the
  zero-inclusive `selectedStepLaw22`. The retained fractional calculation
  first compares its outer measure pointwise with a uniform-value
  reference, charging an outer factor, and only then conditions on the
  number of routed group slots and thins zero labels. The new Lean
  mixture proves the latter operation *within that reference model*.
  The next gate must identify the exact reference measure and prove the
  complete-word domination with its correct factor for the revised RS8
  construction, then prove the reference route's conditional slot law.
  `Structured/Domination.Dominates.expect_le` and
  `ConcreteScalarRouted.routed_scalar_moment` are the old checked pattern:
  transfer any nonnegative route-and-inner test through a pointwise
  measure bound, preserving the same sampled route. The retained RS16
  factor cannot be imported by name. After this bridge, prove that the
  11-coordinate outward operator bounds the full-state kernel. A
  matching matrix shape alone does not discharge either obligation.
- 2026-11-08: Before coding the RS8 outer comparison, audited
  `Structured/Domination.Dominates.expect_le` and
  `ConcreteScalarRouted.routed_scalar_moment`: their setup law is a finite
  PMF, the comparison is pointwise on each complete outer word, and the
  same route-and-inner observable is then carried through expectation.
  The search data enter only after that symbolic reduction. Audited
  `Turbo/RS8ShorteningGeneral.rs8ShorteningEquiv`,
  `RS8PositionSets.rs8ZeroOnSet_card`,
  `RS8SpectrumNumeric.rs8_exact_support_card_closed`, and
  `SampledRS8PacketLaw.sampledRS8PacketColumnsLaw_apply` against the
  revised group sampler. A target atom with `h` nonzero 32-bit columns
  can arise only from messages with that exact column support; each such
  message contributes `(2^32-1)^(-h)`. The shortened-code theorem bounds
  their count by `(2^32)^(h-4)` for `h≥5`. This is slightly weaker than
  the retained `(2^32-1)^(h-4)` count but may still suffice for the
  required comparison with the much larger RS16 reference factor.
  First prove the exact atom formula and this support-count bound in
  Lean. Then verify the factor inequality with exact integers and lift
  the bound through the independent row shuffle and active-group product.
  Do not start the numerical operator import until that reference-law
  connection is checked.
- 2026-11-08: `Turbo/RS8AtomEnvelope.lean` now closes the complete-atom
  comparison for one active RS8 group, including the actual independent
  sampled GF(2^32) column multipliers and the sampled row permutation.
  `sampledRS8PacketColumnsLaw_apply_support` gives exact mass
  `(2^32-1)^(-h)` when the message and target have the same nonzero
  column support, and zero otherwise. `rs8_all_messages_atom_sum`
  identifies the expected atom count with the exact support-fiber count
  divided by `(2^32-1)^h`. `rs8_exact_support_count_le` uses the old
  shortening equivalence to bound the fiber by `2^(32(h-4))` for `h≥5`;
  the RS8 distance floor kills nonzero-message atoms with `h<5`.
  Lean checks that this weaker fiber bound still gives at most
  `1/65535^8 = beta_old/2^256` at every complete 256-bit atom.
  `rs8ExpectedShuffledWordCount_le_old_uniform_atom` carries the
  inequality through the *actual* sampled row shuffle by averaging
  inverse images of its bijective permutations. No numerical receipt
  or admitted fact enters these theorems; the recursive axiom audit
  lists only standard axioms. The old `Dominates.expect_le` pattern is
  now available at the complete group-word boundary. Next multiply
  this pointwise bound over independent active groups, identify the
  inactive zero words, and connect the resulting uniform-value
  reference to the conditional regional-slot law and full inner
  moment. The 11-coordinate operator and finite endpoints remain open.
- 2026-11-08: Before the active-group product and expectation bridge,
  rechecked the completed `Structured/ConcreteOuterCountingProduct.lean`:
  `tupleCount_le_fairRows` multiplies local pointwise count bounds, then
  `activeMessages_expect_le` multiplies by the same nonnegative conditional
  route failure and sums complete row arrays. Also rechecked
  `Structured/Domination.Dominates.expect_le`, which transfers a pointwise
  measure comparison to any nonnegative test. The Turbo counterpart uses
  complete four-bit packet values and the shared regional route plus
  continuous inner state, so only this algebraic structure transfers.
  An initial monolithic Turbo theorem exhausted Lean's kernel memory
  limit at `-M 6000`. Splitting it into local nonnegativity and comparison
  lemmas, a generic finite-product lemma, and a separate weighted-sum
  lemma resolved that resource failure without changing the statement.
  `Turbo/RS8UniformReference22.lean` now proves the active-set product
  bound by `beta_old^|S|` times the reference product density.
  `Turbo/RS8ReferenceMoment22.lean` transfers it through the actual
  `turboRemainingFailure22` event and back to the original shared-outer
  setup expectation. Both modules compile admission-free and their final
  theorems have only Lean's standard axioms. The current reference is an
  explicit atom density: active groups use uniform values on all 256
  bits; inactive groups retain the exact zero-message law.
  Next identify this product density with a normalized reference law and
  derive the conditional regional-slot distribution that justifies the
  producer's `selectedStepLaw22` and `15/16` thinning. Audit the completed
  `ConcretePlacementLaw` and `ConcreteFixedRegionComposition` at that
  interface. The checked thinning identity alone is not yet a bound on
  the original regional route or the 11-coordinate operator.
- 2026-11-08: `Turbo/RS8ReferenceLaw22.lean` normalizes that explicit
  density. Lean proves the complete packet-word cardinality is `2^256`,
  defines a uniform-value law for each active group and the actual
  zero-message law for each inactive group, and proves the independent
  product has exactly the density in the pointwise comparison. It then
  rewrites the actual active-set first-moment bound as `beta_old^|S|`
  times the reference law's expectation of the *same* remaining failure
  event. The theorem is checked with only Lean's standard axioms.
  For the next route component, audited the old
  `Structured/ConcretePlacementLaw` and
  `ConcreteFixedRegionComposition`: the former derives uniform shuffled
  support from the actual uniform permutation and gives exact inclusion
  probabilities; the latter composes independent regional kernels while
  retaining chronological state. Turbo must additionally retain four-bit
  labels and values that can be zero. Prove first that inactive reference
  rows are identically zero and that an active uniform packet row has
  independent uniform four-bit coordinates. Then establish, for a fixed
  active set and one regional permutation, the uniform selected-set and
  value law on complete regional rows. Restrict to a 16-slot physical
  step, condition on the number of selected slots, and identify its law
  with `selectedStepLaw22`; only then apply the checked thinning identity
  and begin the 11-coordinate operator soundness proof. Check that the
  regional shuffles remain independent across the 64 regions while the
  state is continuous across their chronological composition.
- 2026-11-08: `Turbo/RS8ReferenceRoute22.lean` starts that route bridge.
  It proves directly from the RS8 evaluation and sampled adjoint map that
  the zero group produces the zero complete packet word for *every*
  sampled column-multiplier tuple, and that any row permutation preserves
  this word. It lifts the fact through both sampled probability maps:
  the zero-message row law, before and after its row shuffle, has atom
  one at the zero complete word and atom zero everywhere else. All three
  printed theorems have only standard Lean axioms. Next combine this
  deterministic inactive row with the independent uniform active rows
  under the actual regional permutation, retaining all four-bit values,
  before conditioning on the physical-step slot counts.
- 2026-11-08: Before transposing Turbo's reference product, inspected the
  completed `Structured/ConcreteRouteTranspose.transpose_density` and the
  existing `Turbo/RowShuffledArray22.piPMF_prod_pair`/
  `FiniteProductLaw.piPMF_map` proofs. The old transpose merely commutes
  a finite product of independent cells; it is valid for four-bit values
  once both the active and inactive cell laws are identified. The new
  `Turbo/RS8ReferenceRegion22.lean` proves an active reference row is a
  product of 64 independent uniform four-bit packets and an inactive row
  is a product of 64 point masses at zero. It then proves the exact
  transposed array law and the law after the *actual independent regional
  permutations* is a product of 64 complete regional-value laws. The
  theorem retains one shared permutation within each region and all
  packet labels. All printed theorems compile with standard Lean axioms.
  Next derive the one-region law as a uniform selected group-position
  set of size `|S|` with independent uniform four-bit values on that set.
  Use the old `ConcretePlacementLaw` permutation-fiber theorem and the
  generalized version of `PhysicalThinning22.selectedGiven22`; then
  condition on the number of selected positions in each physical step.
  The independent-region law is checked, but independence across physical
  steps still requires that conditioning argument and must not be assumed.
- 2026-11-08: Before the one-region placement component, rechecked the old
  `Structured/ConcretePlacementLaw.shuffle_prob_subset` and
  `ConcreteShufflePermutation.uniform_permutation_shuffleLaw`, together
  with the existing `Turbo/PhysicalThinning22.selectedGiven22_apply`.
  `Turbo/RS8ReferencePlacement22.lean` now proves the complete pre-shuffle
  regional word has point mass `16^{-|S|}` exactly when all its nonzero
  packets lie in the active set. It rewrites the actual one-region law as
  an average over inverse images of the uniform permutation, then proves
  its point mass is the uniform selected-set inclusion probability times
  `16^{-|S|}`. The resulting closed formula is
  `choose(L-k, |S|-k)/(choose(L, |S|) 16^{|S|})` for a target word with
  `k≤|S|` nonzero packets and zero otherwise, where `L=16R`. Thus zero
  values at selected slots remain accounted for. Theorems compile
  admission-free with standard Lean axioms. This is the exact complete
  regional law produced by the reference and actual uniform route.
  Next identify its conditional distribution on a partition into 16-slot
  physical steps: given all selected-slot counts, each step must have a
  uniform selected subset and independent uniform four-bit values, with
  the selected subsets independent across steps. Use the old finite
  product-conditioning machinery only after matching this exact law;
  do not infer step independence from the unconditional one-region mass.
- 2026-11-08: Audited `Turbo/PhysicalThinning22.selectedGiven22_apply` and
  the old `Structured/ConcreteMarkedConditioning`/
  `FiniteProductLaw.piPMF_bind` before introducing the latent selected set.
  `Turbo/RS8SelectedValues22.lean` generalizes the checked 16-slot
  selected-value experiment to an arbitrary finite region length. It
  proves that, given a selected set `T`, every complete word supported
  inside `T` has mass `16^{-|T|}`, including words with zeros at selected
  positions. A fixed regional permutation maps this law exactly to the
  same law on its permuted selected set. Combined with the old checked
  uniform-permutation support theorem, this proves the actual one-region
  reference law equals: first choose a uniform `|S|`-subset, then sample
  independent uniform four-bit values on that subset. The printed Lean
  theorems have only standard axioms. `Turbo/RS8ReferencePhysical22.lean`
  pins the generic selected-value law at length 16 definitionally to
  `PhysicalThinning22.selectedGiven22`.
  The remaining placement gate is the joint conditional law of the
  selected sets across ordered physical steps or 32-slot macro blocks.
  The producer conditions on their counts, uses within-step uniform
  subsets and values, and averages count profiles with exact
  hypergeometric weights. Prove that disintegration from the uniform
  global selected set before applying the existing thinning identity or
  11-coordinate numerical envelopes. The outer reference factor and
  one-region product are already checked; do not reintroduce separate
  inner states at region boundaries.
- 2026-11-08: Before partitioning the selected set into steps, rechecked
  the old `ConcretePlacementSimplexLaw`/`ConcretePlacementCoordinates`
  finite-fiber arguments and the physical `PacketInner22.stepInput`
  indexing. The actual step order uses
  `finProdFinEquiv : Fin R × Fin 16 ≃ Fin (16R)`.
  `Turbo/RS8StepGeometry22.lean` now checks the equivalence between a
  global selected set and its ordered array of 16-slot subsets. It
  proves the global cardinality equals the sum of step cardinalities,
  and that the actual uniform regional selected-set law gives equal
  mass `1/choose(16R,|S|)` to every step-subset array with that total
  count, zero otherwise. All printed theorems have only standard axioms.
  An initial symbolic `Fintype.sum_equiv` cardinality proof timed out
  even after increasing Lean's heartbeat limit. The checked replacement
  uses the old `Finset.card_eq_sum_card_fiberwise` pattern and a separate
  explicit equivalence of each step fiber; it compiles at the normal
  resource limit. Next disintegrate this uniform array by its vector of
  step counts: each count profile receives the product-of-binomial
  hypergeometric weight, and conditional on that profile the step
  subsets are independent uniform subsets of their stated sizes. Then
  factor the selected four-bit values across the same steps and identify
  the resulting conditional complete-word law with the existing
  `selectedStepLaw22`.
- 2026-11-08: Audited the old `ConcreteMarkedConditioning` and
  `SpinCodes/FiniteProductLaw.piPMF_bind` for the next disintegration.
  The old proof conditions independent binary marks on their total count;
  Turbo instead starts from an exactly uniform regional selected set and
  keeps its four-bit labels. `Turbo/RS8StepCounts22.lean` gives the exact
  independent uniform subset law at a specified vector of 16-slot counts.
  `Turbo/RS8StepProfileMixture22.lean` proves the actual global uniform
  selected set is the hypergeometric mixture of those product laws, with
  weight `prod_step choose(16,j_step)/choose(16R,q)` at total count `q`.
  Normalization is derived from the equality of the two finite laws,
  rather than postulated. `Turbo/RS8StepValueFactor22.lean` proves that,
  conditional on the full selected set, complete four-bit regional values
  factor across the physical 16-slot steps. The proof uses the actual
  step indexing and retains selected positions whose value is zero.
  These checked theorems use only Lean's standard axioms. The first
  attempt to simplify the reciprocal product in the value factor left a
  goal; applying the explicit product-of-powers and product-of-inverses
  lemmas closed it. Next assemble the set-profile and value-factor
  identities into a law of complete physical step words conditional on
  their counts, identifying each factor with `selectedStepLaw22`.
  Continue to carry one inner state chronologically across steps and
  regions. The existing 11-coordinate numerical envelope and finite
  endpoint certificates have not yet been connected to this law.
- 2026-11-08: Checked `Turbo/RS8ReferenceStepLaw22.lean`. Its theorem
  identifies the exact normalized-reference one-region route, after the
  actual physical-step coordinate change, with a hypergeometric mixture
  of independent `selectedStepLaw22` inputs. It assembles the checked
  global selected-set law, conditional value factor, and the old
  `FiniteProductLaw.piPMF_bind` theorem; the 16-slot factor agrees with
  the preexisting `PhysicalThinning22.selectedStepLaw22` definition.
  The printed theorem uses only standard Lean axioms. Conditional on the
  count vector, physical-step words are independent; the common count
  vector is still integrated at the region level. Before the next
  component, audited old `ConcreteFixedRegionComposition`:
  `matrixExpect_pi_product` composes independent region kernels while
  retaining one state, and `pathKernel_regionStream` identifies the
  chronological stream with the matrix product. Turbo's `PacketInner22`
  likewise serializes region then step and uses one 22-bit state, but
  the new inner step law includes full four-bit packet values and a
  shared regional count profile. Next prove the array-level route
  identity across all 64 independent regions and connect its ordered
  word stream to `stepInput`/`flatInput`. Then state the exact full-state
  transfer-matrix expectation before attempting the 11-coordinate
  numerical compression or fractional two-step bounds. Audit the old
  kernel soundness and checked numerical producer at that boundary.
- 2026-11-08: `Turbo/RS8ReferenceStepArray22.lean` now checks the full
  64-region input law: independent regional profile mixtures, with
  independent complete physical-step words conditional on each region's
  count vector. It applies the already checked independent-region route
  factorization via `FiniteProductLaw.piPMF_map`. The same file proves
  `stepInput` and `flatInput` consume precisely these words in the
  encoder's region-major chronological order. Both printed statements
  use only standard Lean axioms. No step or region resets the 22-bit
  state. Next audit the old `ConcreteFixedRegionComposition` kernel
  product and the retained Turbo 32-slot fractional producer, then
  build the exact full-state weighted moment using this law. In
  particular, do not multiply stepwise scalar failure probabilities:
  matrix products must retain the shared state, and the producer's
  two-step fractional inequality needs its own soundness theorem.
- 2026-11-08: Audited the completed Scalable implementation of the old
  region-composition idea in `UniformMixer/ScalableRegionComposition` and
  `ScalableRouteAverage`. Its `inputListKernel_flatten`,
  `kernelListProduct_matrix`, and `matrixExpect22_pi_product` are already
  checked for the same 22-state space and preserve chronological state.
  Also audited `FiniteMoment.prob_weight_le_of_moment` and the retained
  Turbo producer `workstreams/k18_margin60/replay_fractional.py`: the
  producer thins before taking fractional powers of chronological
  two-step products, averages with exact 32-slot hypergeometric weights,
  and powers the regional operator 64 times. These are further proof
  obligations, not consequences of the local input-law identity.
  `Turbo/RS8ReferenceKernel22.lean` now checks the complete semantic
  bridge: the actual flat input is a concatenation of region streams;
  each region stream is an ordered product of full-state physical step
  kernels; conditioning on the regional count profile gives a product
  of `turboSelectedKernel22` matrices; and averaging the 64 independent
  regions gives their chronological matrix product. It then integrates
  the sampled state maps, applies the lower-tail moment inequality, and
  bounds the actual shared-outer active-set first moment by this
  reference kernel with the existing RS8 outer factor. All printed
  theorems compile without admissions and use only standard Lean axioms.
  This identifies the exact finite-certificate target, but the full
  `2^22` matrix is not a practical numeric checker. Next inspect the
  retained 11-coordinate map and its source invariants, compare its
  domination proof with the old Structured/Scalable quotient or
  Collatz arguments, and prove a sound compression from each full-state
  `turboSelectedKernel22` (including selected zero labels) to the
  certified 11-coordinate operator. Only then formalize the producer's
  two-step fractional inequality and check the finite endpoints and
  whole-code occupation union.
- 2026-11-08: Audited the actual `s16_closure/sparse_kernel.py` row
  semantics and `kernel_birth_density.py` capped-density split before
  defining a Lean representation. Coordinates 0, 1, and 2 represent
  zero-state mass, arbitrary nonzero mass, and nonzero mass dominated
  pointwise by uniform density; coordinates 3..10 represent arbitrary
  measures supported in the eight positive expansion-weight classes.
  The producer's selected newborn row must be checked as a *complete*
  decomposition; entrywise minima across different candidate rows are
  not justified. `Turbo/PhysicalCone22.lean` formalizes this exact
  measure cone, checks that the eight class indices cover every nonzero
  state under the production map, and proves nonnegativity and the
  terminal total-mass bound. It defines a local complete-row soundness
  condition and proves that any matrix meeting it propagates the cone
  through the full `2^22`-state kernel. The initial zero state has a
  checked cone representation, and full-state and envelope actions
  both compose in chronological matrix order. All printed theorems use
  only standard Lean axioms. **No retained numerical row has yet been
  shown to satisfy the soundness condition.** Next derive symbolic row
  witnesses from the exact physical input law: first the zero-state
  birth row and a coarse arbitrary-nonzero row, then the uniform and
  weight-class rows using the checked uniform-GL refresh law. Compare
  each bound with `sparse_kernel.outward_at_z`; only after the symbolic
  inequalities close should a finite checker authenticate its Arb
  endpoints or capped-density candidates.
- 2026-11-08: Before the first local row, rechecked the completed
  `UniformMixer/ScalableWeightedRows.bernoulliStepKernel_zero_entry` and
  `Kernel.update_zero`: a zero entering state is never refreshed by
  uniform GL, so its next state is exactly the feedback of the input.
  `Turbo/PhysicalBirth22.lean` identifies the corresponding exact
  zero-inclusive selected-packet birth measure. It proves that the
  zero-source row of `turboSelectedKernel22` equals this feedback law,
  that every birth mass is nonnegative for `z≥0`, and that its total is
  the expected tilted output weight. The theorems compile with only
  standard axioms. Next split this exact birth measure into the
  producer's zero, uniform-density, and eight weight-class pieces;
  prove that *one complete selected row* bounds the resulting pieces.
  A coarse arbitrary-nonzero row may be a pilot, but the retained
  numeric certificate needs the sharper complete row and its source
  arithmetic checked.
- 2026-11-08: Extended `Turbo/PhysicalBirth22.lean` from the pointwise
  birth law to a complete symbolic zero-source row. A deterministic
  class index uses the exact production-map weight spectrum; Lean proves
  that every state receives one class, that the arbitrary and uniform
  coordinates (1 and 2) receive zero birth mass, and that the class
  masses sum to the exact tilted birth moment. For every nonnegative
  tilt and every zero-supported source piece, the exact full-state
  `turboSelectedKernel22` output has a checked decomposition into the
  eleven cone coordinates, with the complete row of exact class masses
  as its bound. This is admission-free and uses only standard axioms.
  It matches the structural form of `sparse_kernel.outward_at_z`'s
  source row; it does **not** authenticate the producer's reported
  numerical values. Next prove sound symbolic live-state rows (arbitrary,
  uniform-dominated, and each weight class), using the checked
  uniform-GL refresh and exact tilted feedback moments. Then connect
  all symbolic row expressions to the producer's histogram/Walsh data
  through a checker, including any selected capped-density birth row.

- 2026-11-08: Audited `UniformMixer/ScalableWeightedRows` and
  `ScalableEnvelopeRows` before defining the Turbo live rows. The old
  Bernoulli-input law cannot be reused: `Turbo/PhysicalLive22.lean` instead
  retains the exact selected four-bit packet input in both the tilted
  output moment and the feedback-conditioned return-to-zero moment.
  Lean checks the uniform-GL live entry, its exact zero return, the
  pointwise live-atom upper bound, and exact row sum. It also checks the
  weighted source bound for each of the eight expansion classes and the
  arbitrary-nonzero coordinate by a finite class maximum, while the
  uniform-density coordinate uses the live-state average. These lemmas
  use only standard Lean axioms.
- 2026-11-08: `Turbo/PhysicalLiveRow22.lean` now checks the full live-source
  output measure, including its zero atom, total mass, and pointwise
  atom bound. It splits the output into the exact zero atom and a
  uniformly dominated nonzero remainder. The complete row witness is
  sound for every nonzero source coordinate. `PhysicalSymbolicRow22.lean`
  combines this with the earlier zero-source birth row and proves
  `turboSymbolicRow22_sound`, a complete 11-coordinate domination of the
  actual selected full-state physical kernel. Its axiom audit lists only
  `propext`, `Classical.choice`, and `Quot.sound`. This closes the
  **symbolic local operator**, not the finite numerical certificate.
  Before adding the next checker component, audit the old numerical
  row/census bridge and the exact coefficient definitions in the Turbo
  producer. Determine which retained rows dominate the symbolic row
  entrywise, and which use a different complete decomposition. In
  particular, the producer's feedback-return upper bound and the
  selected capped-density birth row require separate sound proofs.
  Check only complete rows, then connect the resulting operator to the
  chronological two-step and 64-region bound and numeric endpoints.
- 2026-11-08: Audited the completed Structured numerical bridge before
  designing the Turbo checker. `Structured/FiberNumericsBridge.lean`
  derives actual kernel/fiber counts from a concrete map spectrum and
  checked integer Fourier sums; `ConcreteCounts.lean` replaces those
  premises with the exhaustive map-spectrum replay and frozen data.
  `ConcreteFixedNumericSound.lean` separately proves fixed-point interval
  checks imply real rate inequalities. The producer remains untrusted;
  Lean checks its exact arithmetic and the semantic bridge. For Turbo,
  `s16_closure/sparse_kernel.py:outward_at_z` sets the birth row from
  complete class masses and live rows from histogram moments. The live
  return entry is the minimum of the first moment and a Cauchy bound
  using the second moment and feedback-nonzero probability; the uniform
  source uses histogram multiplicities divided by `2^22-1`. The existing
  symbolic row bounds the exact return expectation, so the Cauchy
  inequality and a checked identification of all producer histograms
  are necessary before asserting domination by the retained row.
  `kernel_birth_density.py` supplies a different complete birth-row
  decomposition with a capped uniform-density piece; it cannot be
  justified by entrywise minima of unrelated rows. Next prove a generic
  monotone-row transfer from `turboSymbolicRow22_sound`, then construct a
  small checked numeric pilot from exact selected-law data. Do not treat
  receipt hashes or Arb outputs as proof premises.
- 2026-11-08: `Turbo/PhysicalRowDominance22.lean` checks the generic
  numeric interface: if every entry of one outward matrix dominates
  the corresponding entry of the complete symbolic row, then that
  matrix satisfies `physicalRowSound22` for the actual full-state
  selected kernel. Both componentwise cone monotonicity and row
  transfer use only standard Lean axioms. This interface covers one
  complete row choice at a time; a capped-density birth row requires
  its own symbolic decomposition. Next build a pilot that authenticates
  exact birth masses, first/second live moments, and feedback-return
  bounds at one packet occupancy and one rational tilt, using the
  actual production maps and selected input law.
- 2026-11-08: Corrected the numeric interface after tracing the retained
  K18 replay through `k18_margin60/tune.py`,
  `k16_design/packet_inner_quadratic_extension.py`,
  `s16_closure/sparse_kernel.py`, and
  `k18_margin60/replay_fractional.py`. `sparse_kernel.outward_at_z`
  produces 17 local matrices conditioned on **exactly j nonzero
  packets** (`packetStepLaw22 j`); the replay thins them with weights
  `choose(j,k) 15^k / 16^j` to model selected slots that can be zero.
  Our checked `turboSymbolicRow22_sound` concerns the latter selected
  law (`selectedStepLaw22 j`), so comparing its entries directly with
  the unthinned retained rows would be invalid. The checked
  `turboSelectedKernel22_mixture` already supplies the full-state
  probabilistic identity. For the numerical bridge, first establish a
  symbolic complete row for the exact-nonzero law and compare its
  entries with each retained physical matrix. Then prove that
  nonnegative convex thinning preserves complete-row soundness and
  connects those matrices to `turboSelectedKernel22`. The existing
  selected symbolic row remains an independent semantic cross-check.
  For the pilot use one nontrivial occupancy, such as `j=1`, and a
  rational tilt. This is a scope correction, not evidence of a false
  finite certificate.
- 2026-11-08: Audited the old expectation and chronological product
  interface in `UniformMixer/ScalableRouteAverage.lean` before proving
  the Turbo thinning bridge. `Turbo/PhysicalThinningRows22.lean` now
  checks that nonnegative weighted sums preserve complete-row
  soundness, including support and the capped-density constraint.
  Combining this with the already checked full-state identity
  `turboSelectedKernel22_mixture` proves that sound exact-nonzero
  occupancy rows thin to a sound selected-slot row. This matches the
  retained replay's order and uses only standard Lean axioms. The
  premise that the retained 17 exact-nonzero numeric matrices are sound
  remains open; the existing selected symbolic row is a cross-check,
  not a substitute for that premise. Next derive exact-nonzero symbolic
  rows from `packetStepLaw22`, then test one nontrivial numerical row.

## Review checkpoint: next work in dependency order

The 2026-11-08 review found no contradiction in the new local row
statements. It did find that the progress reports understated the remaining
mathematical work by calling it only numerical. The checked selected-slot
row is useful, but the retained producer starts with exact-nonzero
occupancy; the late discovery of this distinction caused avoidable work.
The old-proof audit must trace the actual producer and consumer together
before choosing the next theorem interface.

1. Pin the retained fractional certificate's conditional inequality on
   paper before extending the numeric checker. `replay_fractional.py`
   thins exact-nonzero rows, multiplies two chronological steps, takes
   entrywise powers, averages placements, and composes 64 regions. The
   justification must condition on the appropriate routing information,
   cap the fixed-active-set bad-event probability at one before applying
   a fractional power, and keep the active-set union factor outside that
   power. An unconditional first-moment bound alone is insufficient.
   Audit the old route-conditioned/fractional proof for this exact order.
   The reviewed source is `research/workstreams/k18_margin60/SLACK_REVIEW.md`
   (conditional count through the regional coefficient formula), with
   witness scope in that directory's `README.md`. It already gives the
   argument: independence permits outer domination after fixing routing;
   conditional Markov and `min(1,x) <= x^alpha` cap each active-set event;
   subadditivity on scalar matrix paths permits entrywise powers after
   chronological two-step multiplication. Hypergeometric averaging and
   the active-set union stay outside those powers. Inspection agrees
   with `replay_fractional.py`; the Lean connection remains to be built.
2. Reuse the local row proof through a generic finite packet input law,
   then instantiate the exact-nonzero and selected-slot laws. Avoid
   duplicating the live/birth/cone proofs merely to replace one sampler.
   Preserve the checked selected theorem as a regression check.
3. Prove the capped-density birth decomposition and the Cauchy return
   bound used by the retained matrices. Check a small nontrivial pilot
   against the actual map census. The symbolic class-only birth row does
   not authenticate a chosen capped-density row by entrywise comparison.
4. Complete numerical authentication and composition for K18/K20, then
   sparse one/two-group coverage, the K20 exceptional event, and the final
   positive occupation sum. K16's 15-state map is a separate remaining
   formalization. Do not describe the full finite theorem as nearly
   closed merely because the generic local matrix interface compiles.

The paper currently states the finite margins using externally checked
certificates; the Lean work has not refuted them. Neither compilation of
the local lemmas nor this review is an independent validation of those
external numerical certificates. Closure confidence should be reported
separately for the analytic reductions, census authentication, composition,
and endpoint checks.
No reviewed derivation currently forces weaker bounds or a changed
construction. Proving the retained inequalities should retain their bounds;
possible margin loss is a contingency if an actual defect or a changed
numeric representation demands it, not a finding of this review.

### Retained K18 argument audit (2026-11-08)

Scope: read `k18_margin60/SLACK_REVIEW.md`, the actual
`replay_fractional.py`, `verify_fractional.py`,
`k16_design/packet_regional_power.py`, and the shared `sparse_kernel.py`,
`kernel_birth_density.py`, and `occupancy_birth_classes.py` producers.
Compared them with the current cone, symbolic row, and thinning interfaces.
This audit did not regenerate the 22-bit census or all retained endpoints,
and does not cover all K16/K20 sparse producers.

- **Routing conditioning:** fixing an active-group set and its potential
  step counts leaves the outer field maps and independent GL updates
  unchanged. Under the dominating uniform reference, packet labels are
  independent uniform four-bit values, including zero. Conditional on
  step counts, support subsets are independent across disjoint steps;
  the counts themselves have the exact fixed-population law. This
  supports the conditional count `beta^q z^-D e0 product(P_J) 1`.
  The proof must condition the event for one fixed active set; the
  number of active sets remains outside the fractional power.
- **Fractional composition:** expanding the nonnegative matrix product
  into scalar paths and using `(sum a)^alpha <= sum(a^alpha)` for
  `0<alpha<=1` gives the stated entrywise-power bound. No disjointness
  of the eleven envelope coordinates is required. Replay multiplies
  each chronological pair first, then powers, then averages its exact
  hypergeometric split. Its ordered matrix-polynomial coefficient and
  64-region product match the notes. Truncation cannot alter lower
  coefficients since all polynomial degrees are nonnegative.
- **Live return rows:** for each entering state the exact return is
  `E[z^weight 1{Cx != 0}]/(2^22-1)`. Cauchy bounds its numerator by
  `sqrt(E[z^(2 weight)] Pr[Cx != 0])`; the first moment is another bound.
  Class maxima and the uniform-source mean in the producer are valid.
  For the mean row, apply Cauchy jointly over a uniform nonzero state and
  input (or Jensen after the statewise inequality). No independence
  between feedback and output weight is assumed.
- **Capped birth row:** write the exact nonzero birth mass as
  `min(m_s,d) + max(m_s-d,0)`. The first term is uniformly dominated
  with coefficient `(2^22-1)d`; the residual stays in its expansion
  class. If `U_s >= m_s`, its class mass is bounded both by the baseline
  class mass and by `sum max(U_s-d,0)`. Therefore the per-class `min`
  inside `outward_rows` is valid for this *same* decomposition. It is
  not the invalid operation of mixing entries from different cap choices.
  The signed Walsh transform's common error allowance bounds either
  sign of every coefficient error before the nonnegative clipping.
- **Numerical trust boundary:** `verify_fractional.py` checks receipt
  identities, geometry, source hashes, coverage, and exact endpoint sums;
  it does not independently derive local census data or recompute all
  fractional endpoints. Those remain actual obligations for Lean.
  The paper's finite appendix currently describes transferred endpoints
  without spelling out the conditional fractional calculation; include
  that argument in the later paper-completeness pass.

No invalid inequality was found in these inspected paths. Ran the existing
`test_replay_fractional.py` (6 tests), `test_sparse_kernel.py` (4 tests), and
`test_kernel_birth_density.py` (3 tests): all 13 passed. They check exact
rational comparisons, chronology, tiny positive values, small-state rows,
and capped-density decompositions. They are implementation evidence, not
a replacement for the full-state certificate. The next proof work should
follow these exact interfaces; there is no audit finding requiring new
parameters, weaker distance, or new numerical optimization.

### Lean architecture audit (2026-11-08)

This review inspected the Lean statements and definitions, not only the
retained Python argument. Scope included actual sampled setup, reference
kernel, physical cone and complete rows, thinning, dense leaf/cover
soundness, and the selected-outer asymptotic conclusion. It is not a fresh
clean rebuild of every module or a completed check of all finite margins.

Assessment: **7/10 overall direction**, with **8/10 for the reviewed proof
architecture** and **5/10 for execution order and efficiency**. These are
engineering judgments, not percentages complete or correctness probabilities.
Keep the architecture; change the order of completion.

1. **Substantial checked progress:** the full generated Scalable dense
   domain now compiles. `ScalableDenseFullCover.rates` prints only the
   standard `propext`, `Classical.choice`, and `Quot.sound` dependencies.
   The exact rational leaf and rectangle-gluing approach follows the old
   hybrid proof: external search proposes data and Lean checks it.
   An audit wrapper applying it to the existing total-layer theorem also
   compiled successfully, with only those same standard axioms. This
   proves convergence of the complete selected-good-outer first moment;
   it does not yet establish the original-setup distance/rate theorem.
2. **Good Turbo semantics:** the sampled code is shared across messages;
   complete full-state kernels retain actual input/state correlations.
   The eleven-coordinate cone describes a decomposition of measures,
   including a uniformly capped nonzero component. Its coordinates need
   not be disjoint events. Complete-row soundness composes without
   inventing independence. The full-state matrix is a semantic object,
   not a matrix to materialize in the numeric checker.
3. **High-priority missing connection:** the retained finite certificate
   needs route-conditioned probability capping, chronological pair
   multiplication, fractional powers, and exact regional averaging.
   `turboSetup22_activeSet_le_reference_kernel` alone bounds an
   unconditional first moment. It does not prove that fractional
   calculation. Formalize the already audited retained argument and
   expose its remaining obligations in a conditional final theorem.
4. **Avoidable interface detour:** the symbolic complete row currently
   treats selected slots including zero, whereas the retained numerical
   producer starts with exactly j nonzero slots and then thins. The
   thinning theorem is sound but still assumes those physical rows.
   Generalize the existing row reasoning over its input law and reuse
   it; do not build a second parallel collection of local lemmas.
5. **Numerical authentication still matters:** symbolic maxima and
   moments do not certify the retained matrix entries. The capped birth
   decomposition must be proved directly; entrywise domination of a
   different, class-only birth decomposition is not a substitute.
   Authenticate the census and outward bounds, then the endpoints and
   complete occupation coverage. K16's separate state map remains a
   separate obligation.
6. **Delivery is unfinished:** the generated full cover and audit wrapper
   currently live under ignored `.lake`. Preserve reproducible generation
   inputs/instructions and a normal named final theorem before declaring
   the development delivered. This audit used existing dependency builds;
   a clean replay remains a final release check.

No reviewed Lean statement revealed a construction counterexample or a
reason to weaken the claimed parameters. Nor does this audit establish
all finite numerical claims. The strongest next action is to finish the
Scalable original-setup assembly, then close one K18 certificate from the
actual setup through the final failure bound. Expand to other profiles
only after that complete path is checked.

Audit execution: `.lake/fullcover/AuditClosure.lean` imports the full
cover, the conditional Scalable conclusion, and the Turbo thinning and
reference-kernel modules. Its final run exited successfully and printed
only the three standard axioms for `audited_scalable_total`,
`turboSymbolicRow22_sound`, `turboThinnedRow22_sound`, and
`turboSetup22_activeSet_le_reference_kernel`. The first wrapper attempt
had an unresolved `Lsched` namespace; the corrected wrapper opens
`Spin.Structured`. Only the successful corrected run counts as evidence.

### Scalable closure delegation (2026-11-08)

The user requested Sol subagents for separable proof work, with the primary
agent retaining strategy, integration, and audit. Scalable is the immediate
focus; Turbo is deferred until this assembly is settled. Three Sol workers
are assigned non-overlapping ownership:

- `scalable_final`: audit the old conditioning/final-probability argument
  and assemble the revised original-setup failure theorem, conditional on
  the already checked dense-rate interface during local iterations.
- `scalable_rate`: audit the old code/rate connection against the actual
  revised encoder; establish necessary code properties without importing
  the old IMT encoder by accident.
- `scalable_delivery`: audit generation inputs and build reproducibility,
  and prepare minimal delivery fixes without rebuilding the dense cover.

The primary agent reviews theorem scopes against the paper and joins the
results to the checked dense cover. Each worker must identify the old
counterpart before a new component, preserve complete proofs and standard
axiom dependencies, and avoid editing this shared plan. No construction,
parameter, or paper changes are part of this pass.

Checked progress in this delegation pass:

- `ScalableOriginalSetupFinal.lean` compiled. Its original-law failure
  convergence theorem has only `ScalableDenseRates eta C`, `eta>0`, and
  `C>0` as numerical hypotheses. It reuses the old
  `Setup.prob_bad_le_cond_sum` and selected outer tail, charging outer
  selection failure once. Independent read-only review found the joint
  setup and shared-message semantics match the revised construction.
- `ScalableCodeProperties.lean` compiled. The actual output vector has
  the existing weight recurrence, is additive, and is causally injective
  for each fixed entering state. The old `ConcreteNativeCodeword` pattern
  is reused, but no IMT-specific encoder theorem is transferred.
- `ScalableInputEmbedding.lean` compiled. It reuses old outer/route
  injectivity and proves the revised 64-bit reshape injective and additive
  under support symmetric difference.
- The above printed proof dependencies contain only `propext`,
  `Classical.choice`, and `Quot.sound`. Full linear-code/rate assembly,
  minimum-distance identification, and instantiation with the complete
  dense cover are the next checks. Circuit-cost claims are not implied
  by these algebraic statements.
- Delivery audit found a tracked search recipe for regenerating the
  ignored candidate. Exact historical bytes are not guaranteed across
  unpinned numerical-library versions; every regenerated proposal still
  requires Lean checking. The builder supports explicit checker refresh
  and accounts for generated dependencies in incremental assembly. Root
  review caught and the worker fixed absolute-path handling in the manifest
  consumer. Root made checker refresh opt-in: a resumed shared-dependency
  worktree must not unexpectedly fetch local Lake packages or rebuild old
  proofs. Existing-object mode requires current checker dependencies.
  No full rebuild or raw-data commit was made.

Final closure checkpoint:

- `ScalableLinearCode.lean` compiled and establishes the actual complete
  encoder's vector/weight agreement, linearity, injectivity, output
  dimension, code dimension, and exact rate. The proof follows old
  `ConcreteNativeLinearCode`, reusing only the unchanged outer and route.
- `ScalableDistance.lean` follows old `ConcreteNativeLinearCodeDistance`:
  it constructs the actual finite codeword set, equates low nonzero weight
  with close distinct codewords, defines the attained minimum distance,
  and identifies its failure event with `scalableBad`.
- `ScalableCertified.lean` compiled successfully and emitted its normal
  `.olean`. It instantiates the original-law and minimum-distance theorems
  with the complete checked cover, expressing the threshold as
  `floor(0.11*Nsched m)`. Its three axiom prints (original-law failure,
  minimum-distance failure, and exact rate) are all exactly
  `[propext, Classical.choice, Quot.sound]`. The successful log is
  `.lake/scalable_certified_audit.log`; failed preliminary compiles do not
  count as evidence.
- Independent scope review confirms native-length distance/rate closure.
  Ordinary/transposed circuit costs, the requested-length wrapper, C++
  equivalence, and finite Turbo certificates remain separate. Next:
  finish the remaining natural-length claims and prepare a controlled clean
  replay, without reopening completed numerical search or changing the
  construction. The later user instruction drops wrapper formalization.

### Remaining-claims audit and circuit plan (2026-11-08)

User instruction: audit and finish Scalable's remaining claims, including
linear ordinary/transposed work, then proceed to Turbo. Natural lengths
alone are sufficient; do not build the arbitrary-length wrapper.

Three Sol audits found no existing formal encoder-cost model or circuit
transposition theorem in the old development. Existing route transposition
and Fourier adjoints are not encoder complexity proofs. Reusable old
components are the concrete Golay fold, accumulator recurrence, outer
linearity, native reshaping, and the already checked revised encoder.
The missing cost proof therefore requires a small new circuit layer.

Chosen approach: finite-coordinate binary linear circuits with explicit
dense dot-product macros, permutation wiring, composition, disjoint parallel
composition, and tiled products. A dense macro is charged at least twice its
input/output dimension product, accounting for coefficient bit products and
XORs; arbitrary large linear maps are never assumed free. Wiring is free in
the arithmetic circuit model. A structural transpose must preserve the cost
and satisfy the actual dot-product adjoint identity. Build a chronological
state scan from those constructors, with summed local step costs, and use
it both for prefix accumulators and the fixed-width recursive inner.
Tie the composed circuit's evaluation to `scalableEncoderLinearMap`.
Do not substitute a disconnected arithmetic cost formula for this witness.

Assignments: `scalable_rate` owns generic circuit/transpose infrastructure;
root owns the shared scan and actual inner assembly; `scalable_delivery`
owns Golay/outer circuit agreement; `scalable_final` closes the missing
schedule upper error bound. Existing `Schedule` already supplies growth,
admissibility and vanishing relative gaps. Minimality should additionally
give eventual `0 <= b - (39/4) log2 N < 24`.

Complexity scope: the theorem concerns fixed sampled/precomputed setup and
linear arithmetic circuits (plus linear indexed bit movement). It does not
claim linear setup generation or bit-addressed Turing-machine complexity:
the random route itself requires order `N log N` seed bits. Seed cardinality
must not be confused with distinct encoder cardinality. C++ equivalence and
measured performance are separate from this mathematical circuit witness.

At each future checkpoint, update the status and next action here before
reporting completion or changing strategy. Preserve failed candidates and
their reason in research notes when they inform a later decision.

### Completed native-length claims and Turbo handoff (2026-11-08)

The final aggregate `ScalableCertified.lean` compiled after all component
checks. Its eight printed declarations have only the three standard axioms;
the aggregate log is `.lake/scalable_certified_audit.log`. This closes the
native Scalable main-theorem claims under the stated circuit model.

- `LinearCircuit.lean` gives explicit matrix, wiring, parallel, tiled, and
  composition constructors, with a cost-preserving adjoint.
- `LinearCircuitPrimitive.lean` expands matrix macros into scalar
  coefficient-AND/XOR expressions and proves Boolean semantics and exact
  gate counts. This closes the audit gap between an assigned macro charge
  and actual primitive operations.
- `LinearCircuitScan.lean` constructs a chronological state scan from those
  constructors, proving its semantics and sum of step costs. Composed stages
  share wires; no expression-tree substitution is charged as a circuit.
- `ScalableOuterAlgebra` and `ScalableOuterCircuit` connect the actual
  Golay--BAA outer to a `40*N` circuit. The route is an explicit permutation
  of wires. `ScalableInnerCircuit` connects the actual fixed-width update
  to a 14,792-gate step.
- `ScalableComplexity.lean` connects the composite circuit to the actual
  `scalableCodeword`, proves a conservative `272*N` bound, and proves the
  lowered transpose computes its adjoint at the same bound.
- `ScheduleUpper.lean` uses the old schedule's minimality to close the
  eventual upper rounding bound; the old schedule already proves growth
  and vanishing relative gaps.

Independent review checked the actual encoder connections and cost algebra.
The paper now states native lengths and fixed-setup circuit work. No numerical
cover was regenerated, and no claim about C++ or setup time was added. A fresh
dependency replay remains a release/reproducibility check, distinct from the
successful aggregate check using the existing checked objects.

Next is one complete K18 Turbo certificate. Before new lemmas, consult the
retained `k18_margin60/SLACK_REVIEW.md` and `replay_fractional.py` alongside
the existing Turbo modules. Two separable tasks are the conditional fractional
route argument and exact-nonzero physical-row soundness. Keep the capped-birth
decomposition joint, and retain explicit unproved numerical hypotheses until
they are authenticated. Expand to K20 and the separate K16 inner only after
this first complete path is checked.

### Turbo continuation after Scalable closure (2026-11-08)

History was re-audited before these components. `ConcreteNativeLinearCode`
supplies the actual-vector/linearity/rate architecture. The retained
`k18_margin60/SLACK_REVIEW.md` supplies the fractional argument, but describes
the earlier RS[16,8] outer: only its probability algebra and unchanged inner
analysis transfer. The current RS[8,4] atom bounds and current certificate
inputs must be used explicitly. The earlier Scalable proof has no fractional
route theorem to reuse by name.

Checked progress:

- `CodeProperties22.lean` identifies the actual emitted vector with
  `ActualSetup22`'s output weight, proves injectivity for every fixed setup
  with distinct RS evaluation points, and constructs its binary linear map.
  It reuses the continuous inner's checked recurrence, the old RS
  information-set lemma, and actual packet permutations. It does not claim
  equivalence with a concrete C++ field basis.
- `CodeRate22.lean` follows the old range-finrank proof: each RS group
  contains 128 binary message dimensions; the complete encoder has
  `2048*R` input dimensions and `4096*R` output dimensions. Its actual
  injective linear range has exact rate 1/2 for every `R>0`. This covers
  the 22-state K18/K20 model, not the separate K16 inner.
- `ConditionalFractionalMoment.lean` proves probability capping and
  fractional powers within conditional fibers and keeps the active-set union
  outside the power. Its whole-encoder theorem explicitly assumes the
  original-law disintegration and fiberwise bad-count estimates. Those are
  not consequences of the unconditional reference first moment.
  The final checked `turboSetup22_bad_prob_le_activeSet_conditional_fractional`
  permits a separate `P(S)` and `Q(S,j)` for every active set, as required
  by the retained routing-count argument.
- `ConditionalFractionalMatrix.lean` proves finite-sum fractional
  subadditivity and its chronological matrix-path version. In the paired
  specialization, the physical pair is multiplied before the entrywise
  fractional power. The summed terminal row keeps the state across all
  pairs. No routing average or group union is moved inside that power.
- `PhysicalCappedBirth22.lean` proves one joint decomposition into zero
  mass, `min(mu,d)` uniform mass, and `max(mu-d,0)` residual class masses.
  Both class bounds concern that same residual, so their minimum is sound.
  It connects the decomposition to `packetStepLaw22`, the actual exactly-j-
  nonzero input distribution, and the physical kernel's zero-source row.
  The pointwise upper census and numerical sums remain explicit obligations.
- `PhysicalGenericLive22.lean` reuses the law-independent cone reasoning
  for an arbitrary complete packet-input distribution and specializes it to
  the exact-nonzero physical kernel. Weighted first-moment and feedback-
  conditioned return caps `b*M` and `b*R` give one complete row with budgets
  `b*M` in the uniform-live coordinate and `b*R/(2^22-1)` at zero.
  This assumes no independence between output weight and feedback. It does
  not yet authenticate those caps. The old
  `ConcreteEpochVarianceBounds.expect_abs_le_sqrt_square` supplies the
  centered-square method for a later Cauchy argument, but lacks the retained
  feedback-event factor, which still needs its own proof.

These completed modules compiled with standard axioms only. Independent
read-only review found no mathematical mismatch. Review additionally caught
an application-interface issue: potential-slot counts depend on the active
group set S, so the retained application needs S-specific disintegrations.
Conditioning jointly on every S could reveal slot positions and invalidate
the intended conditional uniformity. The initial common-law lemma is valid;
the separate S-indexed endpoint is now checked as well.

Remaining tasks are precise: derive the actual routing
fibers and their conditional first-moment domination; authenticate the
exact-nonzero physical rows, then thin them using the existing checked
theorem; finally check all numerical occupations and sum the endpoint.
No new distance margin is claimed by these symbolic lemmas.

All six new Turbo modules above compiled successfully and emitted normal
objects, with the audited declarations depending only on the standard three
axioms. The final conditional endpoint includes separate disintegrations for
each active set. No proof compiler or benchmark remains running for this
checkpoint. The author paper rebuilt successfully with no undefined references
or citations; the changes describe Scalable's completed scope and do not
upgrade Turbo's numerical verification claim.

### Active goal: K18 and K20 closure (2026-11-08)

The user requested a committed checkpoint and continued Sol delegation until
both finite Turbo endpoints are checked. Close K18 first, then K20, for the
actual ideal setup at 10% distance and the paper's reported margins. A
conditional theorem or an externally reported numeric success is not the
completion criterion. Keep K16 outside this goal and preserve Scalable's
completed theorem. If a claimed margin fails audit, report the discrepancy
explicitly rather than silently weakening it.

The pre-commit review found a UTF-8 corruption in the cover emitter's Lean
tuple delimiters. Restored the intended characters and checked the emitter's
import and delimiter code points. This affected regeneration source, not the
already compiled cover; no cover search was restarted.

The checkpoint was committed as `89b5418a`. Root now owns chronological
physical-product integration; Sol workers own (1) actual routing fibers and
conditional domination, (2) exact-input live/Cauchy rows, and (3) retained
numerical provenance and exact arithmetic checkers. The first product module
reuses `ScalableEnvelopeIterate`'s induction architecture, replacing its
homogeneous three-coordinate power with the existing eleven-coordinate cone
and a chronological list of complete physical rows.

New numerical audit finding: current K18/K20 RS8 endpoints are the receipts
`k18-rs8-gf32-s22-v1.json` and `k20-two-ms-rs8-transfer-v1.json` in the canonical
checkout's ignored `tmp`, not the old RS16 60-bit union. Existing K18 physical
exports contain outward matrix entries and row-choice indices, but omit the
integer census values and syndrome-wise Walsh upper witnesses needed to
authenticate capped birth rows. The compact profile/class census proves the
baseline class masses; it does not by itself prove capped residual sums.
The retained K18 witness selects capped candidates for occupancies 6 through
16. A complete uncompressed set of eleven signed-Walsh upper arrays would be
about 352 MiB, before intermediate certificates. Do not emit that blindly.

Before engineering a large transform certificate, diagnose whether baseline
birth rows suffice with suitable proof witnesses. This keeps the actual code
fixed. Sparse one/two-group calculations use occupancies at most two, whereas
the troublesome intermediate group counts rarely produce six active packets
in one physical step. This is motivation for a bounded numerical diagnostic,
not an established bound. Compare against the exact remaining union budget
for the paper's rounded margins, and do not silently weaken the result.
If caps are essential, investigate compressed exact transform evidence after
auditing old `MapSpectrumBridge` split checks and `PolyPacked` identities.

### Turbo routing and physical-row checkpoint (2026-11-08)

The following new modules compiled with only the standard three axioms:
`RoutingCountDisintegration22`, `RoutingCountFibers22`,
`RoutingConditionedMean22`, `RoutingSelectedFiber22`,
`RoutingFixedValues22`, `PhysicalCauchy22`, `PhysicalCoarseBirth22`,
`PhysicalProduct22`, and `NumericRpow`.
The actual original setup now disintegrates over each active set's own
potential-slot count profile. Conditioning leaves the outer, row shuffle,
and GL updates independent with their original laws. The normalized
reference's fixed-route values and the conditional latent selected subsets
have their exact product distributions. Joining these two statements and
connecting the resulting conditional expectation to chronological physical
kernels is the next routing step.

The physical Cauchy proof retains the feedback event jointly with the first
and second output moments. It supports both pointwise class caps and a
uniform-live mean cap. The complete-row product theorem propagates the same
state through all steps and applies fractional powers only after each
chronological pair. The new coarse birth witness keeps the zero mass and
bounds all other birth mass in the arbitrary-live coordinate. None of these
lemmas authenticates a proposed numerical matrix.

Numerical diagnostic decision: preserve the rounded paper margins, 58.24
bits for K18 and 51.82 bits for K20. With the current K18 one/two-group
endpoints, the remaining larger-occupation union budget is approximately
2^(-66.24324), compared with the retained tail 2^(-75.61228).
At K18's q=51 witness (tilt 0.04, exponent 1), replacing only capped birth
rows costs 0.029 bits. Replacing all nonzero birth rows fails; retaining the
exact one/two-active-packet class rows and using coarse birth rows for j>=3
gives an outward diagnostic of 70.7132 bits for this one occupation. The
same replacement of capped rows costs 1.46 bits at K20's q=615 witness.
These are sampled diagnostics, not complete coverage or Lean certificates.

The compact candidate data now contain 25 packet-weight profiles with
multiplicities and exact birth polynomials from 240 one-packet and 27,000
two-packet inputs. The retained physical row masks agree with `Map22`.
The candidate is only about 6 KB and needs no syndrome-wise Walsh arrays.
Next: test the coarse witness across every K18 occupation with its correct
witness, authenticate the profile classification and small birth census,
then check exact numerical row and coefficient bounds. A file hash or a
successful outward Python replay is not that authentication.

The coefficient lane audits `ConcreteFixedCoefficient` and the retained
`packet_regional_power.py`. The old proof expands ordered products over
finite inputs. Turbo additionally needs the exact selected-count constraint
and paired factors with binomial weights outside the fractional power.
The new recurrence will preserve these features and the single regional
normalizer; binary powering must be checked as ordered convolution, not
as a commutative matrix product.

### Compact witness full numerical replay (2026-11-08)

The coarse-birth alternative now passes complete outward numerical replays
at both target sizes, retaining exact class births only for one or two
active packets and using arbitrary-live birth bounds for three or more.
K18 recomputes all q=3..2048 occupations and retains the existing current-RS8
q=1,2 endpoints. Its new whole bound is 58.2408311085 bits and its tail is
66.4739187234 bits. This preserves the paper's rounded 58.24 claim, not the
old diagnostic's additional decimal digits.

For K20, naively reusing all old witness boundaries fails near q=615.
This is a proof-witness issue: the already available dense witness at
tilt .15, fractional exponent .6, and marker 25/128 applies earlier.
Its normalization was independently compared against the historical
producer and matches the retained physical result at q=616. Using that
same witness for q=507..615 gives a complete tail bound of 100.5938735669
bits and a whole bound of 51.8224664170 bits, including the retained rare,
one-group, and regular two-group components. The exact dyadic union is
below a conservative rational lower bound for the paper's 51.82-bit target.
No encoder or setup distribution changes. No large Walsh witness is needed
by these successful numerical candidates.

These are diagnostics, not Lean certificates. The compact packet-profile
data still require their semantic census theorem; one/two-packet birth
counts require checked evaluation; matrix entries and regional/dense
arithmetic require sound checked upper endpoints. Existing q=1,2 numerical
receipts are not axioms. The sparse lane must connect the exact shuffled
row atoms to one/two-distinct-slot regional placement and the support-class
coefficient calculation. K20's rare outer event must be charged once, with
the regular two-group estimate bounded unconditionally after truncating
the expected measure. It must not condition the setup and then assert
unchanged independence.

The actual routing and fractional tail endpoint lane now also includes
checked `RoutingConditionalSteps22`, `RoutingConditionalKernel22`, and
`RoutingPairedProfiles22`. Root's `MatrixCoefficient22` and
`MatrixPolynomial22` identify the exact constrained ordered sum with its
recurrence and with polynomial coefficients over the noncommutative matrix
semiring. Their binary convolution identity preserves factor order.

Profile-census strategy: the checked production map spans RM(2,6).
Splitting two low variables gives one bit, two affine four-variable forms,
and one quadratic four-variable form. The high-variable restrictions reduce
the census to weight histograms over 2048, 128, and 16 quadratics. These
three histograms are now kernel checked via packed coefficient equivalences.
The production reindexing, uniform restriction fibers, and a small bank of
basis-preserving permutations are being joined to the final 25-profile
census. Do not call the census checked until that join is complete.

### Actual routing fractional endpoint audit (2026-11-08)

Historical proof audit before the new component: `ConcreteFixedCoefficient`
expands ordered finite matrix products and moves only scalar placement weights;
`ScalableEnvelopeIterate` retains the continuous state until the terminal sum.
The Turbo adaptation uses actual potential-slot count fibers, then caps each
fixed-active-set event before its fractional power. Chronological pairs are
formed before entrywise powers, and exact regional binomial profile weights
remain outside those powers.

`RoutingPairedPhysical22` is kernel checked. Its complete physical rows transfer
the full chronological conditional kernel to the eleven-coordinate paired
product, without resetting the state at regional boundaries. Its axiom audit
contains only propext, Classical.choice, and Quot.sound.

`RoutingFractionalEndpoint22` is being checked: it joins the original-law
count disintegration, the support-only conditional first-moment theorem,
physical pairing, and the exact regional DP average. Row soundness and matrix
nonnegativity remain explicit numerical premises. Neither this endpoint nor
the pairing theorem authenticates numerical witnesses or the whole certificate.

`RoutingFractionalEndpoint22` is now kernel checked, including its axiom audit
(only propext, Classical.choice, Quot.sound). For each actual fixed active set S,
its original setup probability is bounded by
`(beta^|S| / z^d)^alpha * sum_k ((pairedRegionDP22 B alpha n |S|)^64 0 k)`.
The denominator and outer domination factor are inside the power; the count
profile average is outside. `RoutingFractionalTail22` is being compiled to sum
these events by occupation, retaining `choose (32*n) q` outside each power and
allowing independent z, alpha, and B per occupation. Its actual-message event
bridge permits the q ≥ 3 tail via `Ico 3 (32*n+1)`.

`RoutingFractionalTail22` is now kernel checked. The occupation union has
binomial placement cardinality outside the fractional power and permits
q-dependent witnesses. The actual-message tail event is dominated under the
same original shared setup; no resampled encoder is introduced. The printed
union theorem has only propext, Classical.choice, Quot.sound. The checked
endpoint and tail still require numerical physical-row premises, coefficient
endpoints, and the separate sparse q=1,2 contributions before a whole bound.
Next routing integration: instantiate q ≥ 3 and feed the checked q-dependent
row/DP witnesses, then add sparse and rare-event budgets exactly once.

### Packet-profile semantic join (2026-11-08)

Historical audit revisited `productionPacked_coordinate`, the small block
checks in `Map22Packed`, and `Map22SpectrumData.Block0.block0_checked`:
checked candidate tables are only used after a semantic evaluation/counting
bridge. The new RM(2,6) route retains that discipline while replacing the
full production-state enumeration with small quadratic histograms and
basis-preserver witnesses. There is no historical affine-zero-set packet
profile counterpart; its new obligation is a checked complete 1024-pair
classification followed by uniform quadratic restriction fibers.

Checked objects now include `Map22PacketProfileProduct` (affine-product
closure and both nibble-weight cases), `Map22PacketProfileAffine` (packed
coordinates and affine-pair sampling equivalence), `Map22PacketProfileBank`
(171 independently checked basis-preserver witnesses), `Map22PacketProfileZero`
(witness reindexing, zero-set cardinality, canonical restrictions and 128/16
fiber-weighted sums), `Map22PacketProfileHistogram` (arbitrary additive
observables of the three checked weight histograms), and
`Map22PacketProfileJoin` (odd profile sum and the three nonempty even
profile sums). These compiled to local oleans with standard axioms.

A direct 1024-state classifier/count expression exceeded the 4GB compiler
limit. Its replacement uses 32 blocks of 32 affine pairs, chunked lookup,
and separately checked kind-list receipts; the class cardinality certificate
is a linear list count. No new production-state experiment is introduced.
`Map22PacketProfileClassifier`, `Classes`, and `Census` are the remaining
join under compilation, so the actual 25-profile census is still open.
Next close the affine-pair class aggregation, then use the production split
equivalence to state the complete actual-map profile observable identity.

### Exact one/two-packet physical birth census lane (2026-11-08)

Historical audit before construction: `Map22PackedBridge` identifies every
packed production-expansion coordinate and weight; `MapSpectrum.block_split_weights`
separates reusable packed images from checked histogram arithmetic. The current
`PhysicalStep22` law is uniform on the actual exact-nonzero packet layer, with
support/label equivalence and cardinality `choose(16,j)*15^j`. The new lane uses
these interfaces without a full-state census or a change to the encoder.

`PhysicalSmallBirth22Packed`, `PhysicalSmallBirth22Inputs`,
`PhysicalSmallBirth22Encoding`, `PhysicalSmallBirth22Data`, and
`PhysicalSmallBirth22CensusDefs` are kernel checked. They establish packed
`A^T x`, zero-state emitted weight, explicit 240/27000 input equivalences,
and agreement of packed class/weight codes with actual `physicalBirthIndex22`.
The data module authenticates 240 single-packet feedback-expansion images.
Two-packet images are XORs of those authenticated images by checked linearity.
The input-coordinate bridge uses a symbolic bit-shift proof, not a brute
function-equality decision over ZMod coordinates.

`PhysicalSmallBirth22Checks` is being compiled to authenticate the exact
class/weight histograms in sixteen small first-support-slot shards. Candidate
integers come from `emit_small_birth22.py`, which visits only 240+27000 inputs
and cross-checks the retained compact candidate; that Python cross-check is
not a proof. Census and exact-law polynomial wrappers are drafted, still
awaiting compilation. Numerical row authentication remains separate.

### Actual sparse rare-event audit (2026-11-08)

Before this component, the sparse lane audited
`ConcreteFixedOneMoment.actual_one_mark_moment`,
`ConcreteFixedRoutedTwo.actual_two_fair_probability`, and
`ScalableFixedSelected.scalableSelectedLayer_eq_route`. The old routed
bounds retain a shared regional permutation for distinct selected marks;
those continuum one/two-bit bounds do not certify Turbo's four-bit shell
endpoints. The Scalable selected-law interface conditions the outer law;
Turbo K20 instead needs an unconditional truncated expected measure.

`Turbo/SparseRare22.lean` is kernel checked. Its rare event is the actual
outer-setup event that some group contains some nonzero message with
sampled packet weight 5..8. The checked bound is
`Pr[rare] <= (16*R) * sum_{v=5}^8 C_v`; at K20 this is 8192 times the
single-group sum. Its complement simultaneously forces packet weight at
least nine for every nonzero group message under the same shared setup.
The coordinate marginal, group bound, global bound, and complement lemma
print only propext, Classical.choice, and Quot.sound.

Elaboration lesson: an ambient `Classical.propDecidable` changes hidden
finite subtype/equality instances and can force normalization of gigantic
finite populations when comparing with `rs8ExpectedPacketCount`. Raising
heartbeats or instance priority did not solve this, and higher priority
also changed the `Pi.instFintype` population. The successful proof retains
the original instances, parameterizes generic union bounds by their event
decidability, and introduces classical decidability only for the precise
bounded existential rare predicates. No finite population is enumerated.

`SparseTruncated22` is being checked to connect this one-time budget to
the actual bad-distance event and its unconditional truncated first moment.
The one/two-group complete-word shell endpoint and checked numerical
coverage remain open; the rare-event result alone closes neither sparse
layer nor a complete K18/K20 certificate.

The semantic packet census gate is now closed at
`Map22PacketProfileCensus.productionPacketProfile_census` (checked 2026-11-08).
For every additive profile observable, it equates the sum over actual
22-bit production states to the complete affine-class/weight-histogram
formula. `Classifier` and `Classes` authenticate all 1024 affine pairs,
the 93/840/90/1 class multiplicities, the nonempty witness alignment and
the empty-zero-set case. `Census` joins these to `productionProfileSplitEquiv`;
its printed axioms are exactly propext, Classical.choice, and Quot.sound.
The explicit 25-entry table from the retained numerical candidate is being
checked against this identity for convenient numerical consumption. The
profile census alone does not authenticate physical rows or final tails.

`Turbo/SparseTruncated22.lean` is now kernel checked. The actual bad-distance
probability is at most the one-time rare budget plus the original joint
expectation of the rare-complement count. Its Fubini identity retains the
outer indicator before summing conditional route/inner failure probabilities.
The pointwise occupation truncation is bounded by the unconditioned sum over
messages whose every active group has at least nine packets. In particular,
this justifies the K20 q=2 heavy restriction without asserting independence
under a conditioned outer setup. All three printed declarations use only
propext, Classical.choice, and Quot.sound. Next sparse step: exact nonzero
regional support/value laws under one common permutation, then chronological
q=1/q=2 shell coefficients and their complete checked witness coverage.

The explicit table gate is also closed at
`Map22PacketProfileTable.productionPacketProfile_table_sum` and
`productionPacketProfile_table_count`. The 25 retained candidate profiles
and multiplicities are now connected to the actual production-state
cardinalities, and their multiplicities sum to 2^22. Both endpoint theorems
compiled with exactly propext, Classical.choice, and Quot.sound. Import
this table for coarse-row numerical authentication; the arbitrary additive
observable identity avoids re-enumerating the 22-bit state space.

### Turbo physical-moment checkpoint (2026-11-08)

The original `ScalableFeedbackFourier` and `ScalableFeedbackMomentExact`
were inspected before extending the feedback proof. Their actual-map adjoint
identity and 22-bit character orthogonality are reused unchanged. The new
obligation is the exact nonzero-packet character factor and occupation
coefficient; the Bernoulli bit factor from Scalable cannot be substituted.

The packet-weight generating factor and the sixteen-slot occupation
polynomial factorization have checked with standard axioms. The final
exact-layer expectation/normalized-coefficient join still has a Lean kernel
resource problem; it is not claimed checked. The complete production
25-profile census is now checked in `Map22PacketProfileTable`, replacing
that formerly explicit hypothesis. The small-birth, sparse-shell and numerical
row connections remain separate obligations. No final K18/K20 endpoint is
claimed by the symbolic tail theorem alone.

### Exact sparse support and shared regional route (2026-11-08)

The sparse lane audited `PacketShuffleMixture.uniformLabeledLaw_apply`,
`RS8SelectedValues22`, `RS8ReferenceRegion22`'s transpose/product proofs,
and `ConcreteFixedCoefficient.fugacity_product_expansion` before adding
this component. `SparseExactRegion22` now checks complete nonzero label
atoms with exact support (15 choices per active packet), their transport
by one common regional permutation, and invariance of the averaged law
under equal support cardinalities. `SparseExactCell22` checks its product
of complete nonzero/zero cell laws. `SparseSupportArray22` checks exact
support-array transposition and the product of independent regional laws,
with one shared permutation acting on every active slot within a region.
All printed declarations use only the standard three axioms.

`SparseShellCoefficient22` also checks the no-carry encoding of the two
support cardinalities into the existing univariate matrix coefficient
checker. A regional Bool-pair has degree `b + 65*c`; across 64 regions its
first coordinate sum is at most 64, so coefficient `u+65*v` (u<=64) is
exactly the separate (u,v) support constraint. This preserves chronological
matrix order and reuses `MatrixPolynomial22`'s checked binary convolution.
It is a semantic encoding; the numerical producer may retain efficient
2D storage. The actual Cv support-measure decomposition and continuous
regional-kernel bridge are the next joins; checked 65/2145 shell witness
coverage is still open.

### Physical profile arithmetic checker design (2026-11-08)

Before this component, audited `ScalableLiveMaxSpectrum` (support and
positive counts turn state maxima into small spectrum maxima),
`ScalableTransferSpectrumModel` (full sum minus the unique zero state,
divided by the live count), and `ScalableSpectrumPolynomial` (checked
multiplicities precede arithmetic expansion). `PhysicalProfile22` reuses
that route with the authenticated 25 packet profiles and with the actual
11 physical weight-class predicates.

The numerical checker audit read `Numeric.FixedDefs`, `Numeric.Fixed`,
and `ConcreteFixedNumeric.rowCheck_sound`. Historical evidence explicitly
rejects reducing rational field arithmetic through kernel `decide` at
scale; integer fixed-point endpoints are the established solution.
`FixedDefs` remains import-light and `Fixed` proves every outward operation
against reals. New `PhysicalProfileNumericDefs22` follows that split. Its
stored coefficient list is updated once per linear factor, avoiding an
exponential recursive coefficient evaluator. The once-proved semantic
bridge uses list coefficient membership and a Horner polynomial.

The retained coarse producer's moment factors are z^r and
(((1+z)^4-z^r)/15), normalized by choose(16,j). The checker moves all 15^j
to the final positive integer divisor, retaining the same polynomial and
reducing intermediate rounding. Fourier factors can be signed; the existing
four-corner interval multiplication handles their cancellations. Integer
cross multiplication checks rational upper/lower targets. No floating logs,
external interval assertions, or standalone candidate equalities become
axioms. Fixed precision must preserve both existing rounded paper margins;
long-tail rescaling/normalization is still a separate future audit.

Small-birth checker resource correction: retain independent first-support-slot
modules `PhysicalSmallBirth22PairChecks0` through `PairChecks14`, with fifteen
down to one 225-input pair gates. Arithmetic candidate definitions now live in
`PhysicalSmallBirth22NumericData` and `NumericCensusDefs`, importing only the
packed map; semantic `Data` and `CensusDefs` retain the actual input/feedback
bridges. This follows the historical Basis/block separation and avoids loading
the unrelated old Majorant graph for each numeric shard. All public declaration
names and input aliases are preserved. The aggregate, census, polynomial, and
marker sources are drafts until their subsequent compile gates pass.

After the import split, compile in this order through the shared gated helper:
NumericData, NumericCensusDefs, PairChecks0..14, PairChecks, Checks, Encoding,
Data, CensusDefs, Census, Polynomial, then MatrixMarker22. Each abbreviated
SmallBirth name has the `PhysicalSmallBirth22` prefix. Rebuild PairChecks0
once against the lightweight imports even though its prior heavy-import proof
artifact exists. Subsequent successful independent shards persist across
continuations. The physical packed/input modules do not need this rebuild.

`PhysicalProfile22` is now fully checked, including generic commutative-monoid
packet-product grouping, actual coordinate weight from the profile, complete
table support, full/live sums and means, all eleven physical class sums and
means, the authenticated class cardinalities, and finite-table class maximum
bounds. The final class-count bridge needed explicit cast-of-sum and
cast-of-list-sum lemmas: broad predicate rewriting across State22 hit a
whnf timeout, and broad norm_cast reached recursion depth. The final proof
reuses the already checked real class sum at the constant-one observable
and casts only the two sums. No full State22 enumeration is performed.

### Compiler resource and checkpoint discipline (2026-11-08)

The source audit found that numeric-only small-birth checks imported the
full semantic setup. The Sol lane is separating packed integer data/checks
from the actual-feedback bridge, following the old checked-data/block design.
Each first-slot pair census now has its own object checkpoint; a stopped
compile no longer discards the whole 120-pair census.

Host measurements showed about 12 GB private commitment per heavy Lean
process, despite a roughly 4 GB working set. Concurrent imports intermittently
failed to read existing object files when global commitment was nearly full.
The local ignored compile helper now serializes these jobs through a named
mutex that is released or abandoned by the operating system on termination.
Source work remains parallel. Do not restart a queued compile as if it failed;
inspect its handle and the actual Lean process. Do not touch unrelated user
processes. A numeric-only pair shard measured about 4.4 GB private commitment,
so it is still using the shared gate.

`PhysicalProfileNumericDefs22` and `PhysicalProfileNumeric22` are checked.
`profileFixCoefficients22_sound` proves each stored integer-interval
coefficient encloses the corresponding coefficient of the real product of
linear factors. Rational upper/lower endpoint checks are proved by exact
integer cross multiplication against `Fix.Mem`. The core endpoint has only
standard axioms. Factor-bank normalization, profile-table adapters and
aggregate/maximum checkers remain under compilation; no numerical row or
final margin is certified by this core theorem alone.

### Exact layer normalization and sparse marker simplification

`PhysicalLayerPolynomial22.packetStepLaw22_expect_layer` and
`packetStepLaw22_expect_coefficient` compiled with only the standard axioms.
A direct application of a generic expectation-to-coefficient identity caused
deep definitional comparison on the 2^64-element packet input type. The checked
replacement first derives the conditional finite sum from the actual pointwise
law, separately rewrites a polynomial coefficient to the same conditional sum,
and joins the equalities. No resource-limit increase or weakened law is used.
The common `degreePolynomial22` definition is used by the output-moment and
signed-character polynomials.

The authorized bounded sparse diagnostics found a much smaller prospective
q=2 witness: common positive marker 24 in both row supports and a single
regional matrix bound. K18 uses tilt 1/625 and gives 105.3712 bits for q=2;
K20 uses tilt 1/2500 and gives 106.4712 bits for the heavy q=2 term. Their
whole numerical sums remain 58.2408311085 and 51.8224664170 bits. These are
candidate calculations, not checked endpoints. The formal route discards
separate support-count constraints using nonnegativity, bounds the total-degree
coefficient by the existing positive-marker theorem, and retains every shell.
K20 still pays the rare low-weight event once and uses the unconditional
truncated heavy expectation. Exact arithmetic and the final joins remain open.

`MatrixMarker22` is kernel checked with only propext, Classical.choice,
Quot.sound. Its generic ordered marker expansion and coefficient inequality
apply to any finite alphabet, including the Bool×Bool sparse common marker.
The paired-region wrapper retains the exact binomial normalization and the
complete matrix state through all64 regions. Matrix power monotonicity and
scalar-power identities are public checked helpers. This authenticates the
positive-marker argument, not any numerical witness or whole probability bound.

Further numeric import reduction: the small-birth NumericData now imports only
core MapSpectrumDefs/Init.List.OfFn, and pair gates use Lean's core kernel decide
tactic. Numeric imports contain neither actual setup nor Mathlib's full tactic
collection; actual Data/CensusDefs still perform the mandatory semantic join.

### Sparse q2 common-marker numerical diagnostic (2026-11-08)

Two authorized bounded proof-numeric diagnostics (not benchmarks) tested
one common positive fugacity for both selected group supports. They used
the current compact physical candidate, exact current RS8 Cv values, and
the original chronological regional placement. No FWHT or new census was
run, and these remain Python/Arb diagnostics rather than Lean certificates.
The first check is `.lake/turbo/q2_common_marker_diagnostic.json`; the
geometry-scaled K20 check is `q2_common_marker_k20_scaled_diagnostic.json`.

At marker x=y=24, K18 tilt1/625 gives a 105.3712048-bit q2 contribution
and whole58.2408311085 bits, retaining every shell including5..8. K20
tilt1/2500 gives a 106.4712223-bit heavy-q2 contribution and whole
51.8224664170 bits, retaining the rare outer budget exactly once and
restricting the unconditioned regular expected measure to weights>=9.
The whole calculations retain the q1 endpoints and the newer coarse-tail
diagnostics. Both exact rounded target comparisons passed in Arb.

The initial K20 use of the unscaled K18 tilt grid failed, as expected from
the fourfold regional length; the separately authorized scaled grid passed.
No numerical experiment is left running. These margins motivate replacing
large q2 coefficient-table certificates by one positive common-marker and
Collatz witness per geometry. The formal path is to discard the separate
(u,v) constraints nonnegatively into total degree u+v, then reuse
`MatrixMarker22.orderedMatrixCoefficient22_marker_le`. This uniformly
covers all2145 shell pairs; it does not assert that any omitted numerical
shell is zero. Checked rational regional rows, marker/Collatz arithmetic,
current-Cv coefficient values and the actual sparse-shell normalization
join are still required before claiming a q2 certificate or whole result.

`SparseSupportMeasure22` is kernel checked: the exact actual Cv/choose(64,v)
support-class weights times uniform15-label laws reproduce every original
row-shuffled word atom and the full active-set array measure.
`SparseRegionKernel22` then checks the exact continuous regional moment and
its lower-tail bound. `SparseActualMoment22` joins these to the original
actual active-set first moment. Their printed declarations use only the
standard three axioms. `SparseSupportCoefficient22` checks exact one/two
row-support sums against the ordered coefficient interfaces.

`SparseExactStep22` is also checked: selected j-subsets with uniform
nonzero labels are exactly `packetStepLaw22 j`; the sparse regional law is
the existing ordered hypergeometric count-profile law mixed with exact
physical-step input laws. These results preserve both distinct slots under
one regional permutation. `SparseSupportPrune22` removes inactive supports
using their exact empty-support atom. Its concrete statement fixes the
canonical subtype Fintype explicitly: `Subtype.fintype` and
`Finset.Subtype.fintype` are extensionally equal populations but are not
definitionally equal, and comparing their large function populations caused
kernel normalization costs. The checked proof does not enumerate them.

Remaining current joins: the generic11/full-state matrix averaging lemma,
regional profile→physical coefficient DP identity and whole-stream envelope
comparison; one/two active-set reindexing and weighted65-shell folds; the
nonnegative separate-constraint→total-degree marker bound; exact current Cv
rational coefficients; and checked rational physical/regional/Collatz
arithmetic for the saved markers. The successful Python markers remain
candidate witnesses until these obligations are discharged.

The pure numeric pair gates now use the helper's whitelisted `-Lightweight`
mode (single light compiler, max2500MB) alongside the default heavy gate.
The largest core-only slot0 shard passed with observed peak private memory
1,164,455,936 bytes and RSS662,974,464 bytes. This is a resource diagnostic,
not a numerical probability certificate. Remaining pair slots are independent
checkpoint modules and do not reload the old Majorant graph.

The polynomial wrapper also drafts the exact j=0 source row. Its audit uses
`RegionalEncoder.stepKernel_empty`, `Kernel.emptyStep_zero`, and
`PhysicalCone22.physicalInitialCone22`. Zero support forces the entire packet
word to be zero, so its original uniform layer law evaluates every observable
at that zero word. Its full kernel is the existing empty-step kernel; from an
entering zero state the birth is exactly delta-zero, with envelope mass1 only
at coordinate0. This avoids the coarse zero-plus-total overcount for j=0.

### Complete q1 compact-family diagnostic (2026-11-08)

The authorized bounded `diagnose_turbo_q1_complete65.py` replay recomputed
all65 shell endpoints from the new compact CoarsePhysical family (exact
j=1 birth and generic Cauchy live rows), rather than assuming that the
retained q1 receipts remain valid after the physical-witness change. It
used all five tilts from the actual K18 receipt and all four tilts from the
actual K20 receipt; no new grid/search was introduced. The ignored checkpoint
is `.lake/turbo/q1_coarse_complete65_diagnostic.json`.

The recomputed q1 margins are58.2456338686490 bits (K18) and
56.2570901999489 bits (K20), matching the retained endpoints to displayed
precision. Combined with the new common-marker q2 candidates, one K20 rare
budget and the newer coarse tail diagnostics, the whole margins are
58.24083110851727 and51.82246641700886 bits; both exact rounded target
comparisons passed. This closes the diagnostic risk that the proposed
uncapped sparse physical witness family loses the requested paper margins.
It remains numerical evidence pending Lean arithmetic authentication.

An authorized exact recombination also replaced K20's heavy-only q2 sum by
the full current-Cv sum over0..64, using the already saved marker witness.
Its q2 margin is102.6824180241 bits and the whole margin remains
51.82246641700886 bits. The unconditional full q2 mean dominates the
rare-complement/truncated q2 mean by nonnegativity. The chosen formal path
can therefore retain the rare budget once while using the same complete
65-shell/2145-pair q2 endpoint family for both geometries; no separate
heavy-row-shuffle transfer is required. The checkpoint is
`q2_common_marker_k20_full_recombination.json`.

`SparseMarker22` now kernel checks the separate-constraint→total-degree
comparison, reuse of the existing positive marker theorem, and the explicit
joint marker matrixR0+2tR1+t^2R2. Its printed declarations use only the
standard three axioms. No numerical experiment remains running in this lane.

`PhysicalProfileNumericFactors22`, `PhysicalProfileNumericAggregate22`,
and `PhysicalProfileNumericTable22` are now checked with standard axioms.
The factors cover first and second moments and signed packet Fourier
coefficients, normalized by choose(16,j)*15^j. The table adapter authenticates
the import-light 25-entry list, its class predicates/counts, the polynomial
profile product, actual State22 class caps/means, and the complete-state
Fourier coefficient average. Raw signed coefficients are aggregated before
the common divisor to retain exact zero cancellations at z=1. Numerical
candidate endpoints have not yet been submitted to these checkers.

These 10^30-scale intervals are for local physical rows. Future regional
DP/marker certificates must retain common choose/t scalars outside the
64th power or use an audited rescaling witness. Rounding normalized tiny
regional entries to absolute 10^-30 and later multiplying a huge z^-d/outer
factor is not a valid route to the final rounded margins.

All15 core pair checkpoint modules0..14 are now kernel checked and persisted:
the240 one-packet inputs and all120 two-slot supports times225 label pairs
(27,000 two-packet inputs). The separate actual-law join is still pending.
The current heavy compile sequence is PairChecks (thin imports), Checks
(finite arithmetic assembly), Encoding, Data (240 actual-image identities),
CensusDefs (actual code bridges), Census (actual-layer counts), Polynomial
(exact normalized class/zero formulas and j=0 source cone). MatrixMarker22
has already passed its complete standard-axiom audit and need not wait for
this remaining birth-law assembly.


The small-birth numeric lane then removed the remaining Mathlib/tactic imports.
All fifteen pure arithmetic shards are now checked. A representative shard
uses about 1.1 GB private memory; the helper permits one such whitelisted
`-Lightweight` job alongside the single heavy job. Semantic wrappers still use
the heavy gate. This is a measured exception to the earlier serialization,
not permission to run arbitrary heavy Lean processes concurrently.

The complete q=1 numerical replay also confirmed the simplified physical-row
family at every one of the 65 shells, using the existing five K18 / four K20
tilts. It reproduces q=1 bounds of 58.245633868649 and 56.257090199949 bits.
The K20 full unconditional q=2 marker bound is 102.6824180241 bits, so the
truncated q=2 contribution can simply be bounded by this full nonnegative
expectation. The whole candidate sums remain 58.2408311085 / 51.8224664170.
No old q=1 receipt is being assumed unchanged without this replay. All these
numerical replays still require the final exact-arithmetic certificates.

The thin `PhysicalSmallBirth22PairChecks` imports and finite arithmetic
`PhysicalSmallBirth22Checks` assembly are now GREEN after the core-only split.
The current restart point is the semantic `Encoding` rebuild, then `Data`,
`CensusDefs`, `Census`, and `Polynomial`. All 120 pair gates, each checking
225 label pairs, remain opaque compiled lemmas; the assembly also checks
the 99-bin sums.
No numeric shard needs to be regenerated or rerun for a semantic-wrapper edit.

### Input coefficient wrapper: finite-instance conversion obstacle

The generic exact-layer coefficient identity is checked in
`PhysicalLayerPolynomial22`. Its concrete `physicalInputMoment22` wrapper
exposed a separate conversion problem. Full elaborated goal comparison
confirmed canonical Pi Fintype instances on both sides, with one difference:
the RHS polynomial definition contains a separately hoisted NeZero(2)
proof inside the nested ZMod Fintype. Broad matching/conversion unfolds
monomials and large function-space enumeration while trying to compare
these representations. A scoped recursion depth 4096, borrowed from the
historical RS8 numeric pattern, reached the 6GB kernel memory limit; it
was removed. This is not a reason to increase resources further.

The next checked interface transports the polynomial's Fintype argument
explicitly using Subsingleton.elim, keeping the actual input law and all
observable functions unchanged. Its helper is public in
`PhysicalInputPolynomial22` so the signed-character polynomial can reuse
it. This helper/wrapper replay is pending. The polynomial product
factorization itself remains checked.

The exact packet-character factor is now kernel checked in
`PhysicalPacketFourier22`: summing all sixteen four-bit labels gives the
signed weight polynomial, and removing zero gives the exact fifteen-label
factor. Printed dependencies are only propext, Classical.choice, and
Quot.sound. This reuses the historical adjoint/character argument rather
than an independent-bit approximation.

A common formalization bottleneck is now diagnosed precisely. Separately
elaborated Fintype instances differ only in a hoisted NeZero proof for
ZMod 2, but definitional comparison can unfold the enormous function-space
enumerator. Scoped depth 4096 still exhausted memory in two joins, so that
workaround was removed. The chosen repair is explicit transport using
Subsingleton.elim on Fintype instances, following the existing proof
irrelevance model; no finite-domain enumeration or trusted axiom is added.
The actual-input polynomial and sparse product joins are being replayed
with this repair. They are not yet counted as closed.

The root's next row adapter preserves the historical split: pointwise
class maxima for arbitrary supported mass, and the actual live-state mean
for the uniformly capped coordinate. `PhysicalProfileWeighted22` and
`PhysicalProfileRows22` are in verification, not certified endpoints.

The semantic Encoding, Data, and CensusDefs rebuilds are now GREEN on the
core-only data graph. Census is being rechecked after replacing a broad
closed-histogram getD rewrite with a generic opaque-function lemma; its
remaining support projections are resolved by scoped cases. Numeric shards
and arithmetic assembly are unchanged and remain GREEN.

The final Polynomial source now drafts hbirth interfaces for actual rows:
`physicalSmallBirth22_zero_exact_cone`, `physicalSmallBirth22_one_cone`,
`physicalSmallBirth22_two_cone`, and their one/two `cone_of_le` wrappers.
The latter accept outward componentwise bounds on the exact normalized
polynomial budget vectors. They partition the same actual birth measure,
following the audited selected-law birth-parts proof. These interfaces are
not GREEN until the pending Polynomial compile and axiom audit complete.

`PhysicalProfileWeighted22` and `PhysicalFeedbackFourier22` are now GREEN
with standard axioms. The latter derives the zero-feedback generating
function and the feedback-active probability from the actual input law.
A scalar Input64 kernel identity is proved before lifting through an
expectation; this avoids deep conversion inside the complete packet-word
sum. The weighted adapter similarly uses explicit congrArg rather than
an unrestricted congr tactic over the full state-space sum.

Endpoint audit: the diagnostic Chernoff cutoffs are integers, not rounded
floating thresholds: 52428 for K18 (N=524288) and 209715 for K20
(N=2097152). These are floor(N/10), and all q1, q2, and tail replays use
these same cutoffs. The exact-arithmetic endpoint must retain them; changing
to ceil(N/10) would spend margin that the rounded claims do not promise.

`SparsePhysicalKernel22` is now GREEN at the default recursion depth.
Its exported scalar, physical-product, exact-region, and normalized-DP
identities use only propext, Classical.choice, and Quot.sound. The region
law averages one shared permutation and exact nonzero packet labels.
The DP divides by choose(16R,j) once, after the physical coefficient fold.

The failed sparse join had a precise cause. The public expectation used
Nat.instNeZeroSucc, whereas the fixed packet law stored its own NeZero
proof inside the PacketStep16 Fintype. Depth 4096 exhausted kernel memory;
that probe was removed. A regular implicit Fintype parameter now extracts
the enumeration from the fixed law before constructing the expectation.
An instance hole would allow eager typeclass search to restore the mismatch.
A checked adapter then exposes the product identity in the public law's
syntax. No enumerator normalization or additional axiom is required.
The eleven-coordinate regional envelope and exact current-Cv coefficient
interface are the next verification gates.

The final q1 composition must retain a separate fugacity for each shell.
The complete65-shell diagnostic chooses among retained tilts before adding
shell contributions. Choosing one global tilt, or minimizing whole q1
sums, does not justify that bound and can lose the K18 margin.
`turboSetup22_activeSet_expect_eq_support_measure` already gives an exact
sum over complete support profiles, with their unconditional expected
weights and complete-valued support laws. The next bound fixes z and its
physical row family separately inside each profile's expectation, before
summing those weights. For one active group, z and the row family depend
only on its support cardinality; grouping profiles then produces the
retained shell-dependent coefficient bound. This composition is a required
verification gate. The common-marker q2 bound uses one fugacity and is
unaffected. No conditioning on the K20 rare-event complement is introduced.

`SparseRegionalEnvelope22` is now GREEN with standard axioms. It applies
physical_list_product_mass_le once to the entire ordered physical stream,
then averages the joint regional count profiles. The continuous state is
preserved across every region. The exact current-Cv coefficient interface
is the pending gate before the remaining sparse shell assembly checks.

A source audit confirmed that the completed q1 replay summed all 65 current
Cv values for both geometries. K20's 56.25709-bit q1 diagnostic therefore
supports the full unconditional q1 path. Together with the approved full
q2 bound, this removes heavy-specific q1 and q2 transfers from the planned
endpoint. The single rare-event charge is retained for the high-q tail.
These numerical statements remain diagnostics until exact arithmetic is
checked against the assembled actual-setup theorems.


The clean `PhysicalInputPolynomial22` replay is now GREEN, as is the official
`PhysicalCharacterPolynomial22` module. The successful coefficient join
first rewrites an opaque equality at the polynomial level, then applies
the already checked exact-layer coefficient theorem. It does not ask Lean
to establish polynomial-enumerator definitional equality under `coeff`.
The character module also factors scalar monomial identities before lifting
them through products. No recursion or memory limits were increased.

`PhysicalSmallBirth22Polynomial` is GREEN, including actual j=0,1,2 cone
interfaces. `SparsePhysicalKernel22` and `SparseRegionalEnvelope22` are
GREEN: the original shared-permutation regions are connected to ordered
coefficient DPs, with cone propagation applied once across all 64 regions.
The profile-to-numeric-bank joins and endpoint arithmetic are still open.

`PhysicalSmallBirth22Polynomial` is GREEN, including its actual normalized
class/zero formulas, exact j=0 kernel and raw-birth cone, and j=1/2 hbirth
interfaces. All printed endpoints use only propext, Classical.choice, Quot.sound.
The actual-row APIs are `physicalSmallBirth22_zero_exact_cone`,
`physicalSmallBirth22_one_cone`, `physicalSmallBirth22_two_cone`, and the
one/two `cone_of_le` adapters for outward polynomial budget vectors.

The next numerical adapter audits `Numeric.Fixed` add/mul/divInt Mem soundness
and reuses the historical Horner pattern and profileFixUpperCheck22. Core
`PhysicalSmallBirth22BudgetDefs`, the shared `PhysicalProfileFamilySmoke22Data`,
and `PhysicalSmallBirth22BudgetCheck` are GREEN. All three retained birth
vectors pass kernel-only integer checks at the existing authorized tiny
family (z=8607079764250578/10^16, endpoint denominator 2^80). Values and exact
zero entries are unchanged. Semantic `Budget` and `BudgetSmoke` remain drafts
until their heavy-gate compile and final actual-cone axiom audits pass.

The sparse absolute512-bit lower-power trace is scoped strictly to the
low-tilt q1/q2 witnesses. Across all retained tilts, lambda*d <168, hence
the diagnostic Chernoff exponent is below243 bits. The sufficient first3
K18 and first2 K20 tilts need fewer than122 bits. The eventual rational
witness must pass the exact trace and its strictly positive final-value
check. High-tilt tail witnesses, including z approximately0.8607, require
explicit exponent/rescaling or a separate checked-log representation.
They must not use the sparse absolute512-bit trace. Root will audit the
historical tail representation before adding that numerical component.

`SparseSpectrum22` is now GREEN with standard axioms. The historical MDS
Ah and actual expected-count generating-function theorem imply the real
polynomial identity. Coefficient extensionality authenticates every Cv,
including shells 0..64 and zero outside 5..64. A checked polynomial map
then identifies Cv with the real cast of its rational coefficient. This
interface uses the current RS8 measure directly, without an RS16 ratio.
The final arithmetic layer will use a common integer denominator rather
than ask the kernel to decide rational arithmetic.

### High-tilt arithmetic audit (after the sparse power-trace scope check)

The absolute 512-bit sparse trace is only for retained q1/q2 parameters;
its positive-final-lower-endpoint guard rejects underflow. It cannot be
reused for the local profile pilot z=8607079764250578/10^16 at a K20 cutoff.

Historical sources inspected: `Numeric/FixedLog.lean` (`flog_mem`,
`flogI_mem`, range reduction and Horner remainder), `Numeric/FixedDefs.lean`,
and `Structured/DenseFourierExactCheck.lean` (`checked_collatz`). The old
proof checks rational transfer inequalities by integer cross products and
uses verified logarithms for exponential conclusions, avoiding direct
small-real powers. That separation should be retained.

For the Turbo tail, an untested substitution of a Collatz bound for the
retained finite matrix power could lose the tight K18 tail budget. Preserve
the ordered coefficient/power calculation. The candidate arithmetic design
is integer mantissas with an explicit shared binary exponent for each
matrix/vector stage; check each rescaling by integer cross multiplication.
Then use the old checked logarithm on bounded scalar mantissas/bases, adding
exponent*log(2), for the large cutoff/exponent comparison. Do not evaluate a
huge rational power or force a tiny tail value onto an absolute fixed grid.
This is an audited next-component plan, not an implemented or checked tail
certificate. Finish the local row family and sparse semantic gates first.

`SparseActiveProfiles22` and `SparseActiveKernel22` are now GREEN with
standard axioms. The actual support sum prunes inactive coordinates and
reindexes the selected groups through an injective labelled embedding.
Each region therefore uses occupancy 0, 1, or 2 of the same shared
permutation. The resulting kernels are the normalized physical DPs,
and their entries are nonnegative when the physical envelopes are.
Verification now proceeds through the support-dependent Chernoff family
and the original setup's occupation-1/occupation-2 assembly.

### Local family authentication checkpoint

`PhysicalProfileRows22` is GREEN, including the actual weighted first and
second moments, coarse birth cone, and feedback-probability complement.
The last conversion compared two decision procedures for the same event.
Allowing definitional reduction expanded the finite probability filter;
an explicit generic DecidablePred equality and an opaque probability
constant avoided enumerating the packet space. No resource cap was raised.

`PhysicalProfileFamilySmoke22Check` is GREEN for all 17 occupations of the
retained rational family. `PhysicalProfileFamilyExtract22` is GREEN and
extracts the individual checks without rerunning numeric arithmetic.
`PhysicalSmallBirth22Budget` and `PhysicalSmallBirth22BudgetSmoke` are GREEN:
the exact retained empty/one/two-packet birth vectors imply cones for the
actual input laws. The whole actual-row family assembly remains in progress.
Return caps already include division by 4194303; dividing twice would be
unsound. The new assembly records that convention explicitly.

### Handoff checkpoint — automated goal cleared

The user requested a good stopping point and a handoff for another session.
`PhysicalProfileAssembly22`, `PhysicalProfileCompleteRows22`,
`PhysicalProfileBirthFamilySmoke22`, and the final
`PhysicalProfileFamilySmoke22` now all pass complete-module checks. The
last wrapper's four declarations use only propext, Classical.choice,
Quot.sound. No task-owned Lean process remains after that check, and all
three subagents stopped. No K18 numerical batch was launched.

Read `TURBO_K18_K20_HANDOFF.md` before resuming. The next semantic blocker
is `SparseShellFold22`, followed by its actual-shell/occupation/q1-choice
dependents. Subsequent finite arithmetic is still open. The scaled-matrix
sources and three-family K18 producer are preserved drafts, not checked
endpoint evidence. The work after commit89b5418a remains uncommitted.

### 2026-11-08 local profile family checkpoint
- `PhysicalInputPolynomial22` now closes at default depth with standard axioms. The expensive boundary was avoided by proving the polynomial equality to its canonical typed `degreePolynomial22` first, then rewriting under the coefficient. No resource-cap increase is retained.
- Integer coefficient caches and integer rational-power upper-bound soundness are kernel checked. The one authorized exact-rational local family at z=8607079764250578/10000000000000000 was generated once; its retained source is `PhysicalProfileFamilySmoke22Data` with common endpoint denominator 2^80.
- `PhysicalProfileFamilySmoke22Check.familyChecked` passed all 17 rows by kernel reduction (axioms: propext only). This authenticates the Bool/Fix checks; the semantic family adapters remain under replay. Return caps are already divided by the live-state count 4194303.
- Routing independently checked the retained exact birth j=0/1/2 vectors against authenticated census coefficients. No witness values changed and no numerical experiment was repeated.
- This local family is not a completed K18/K20 failure bound. High-tilt tail arithmetic must preserve exponent scaling, rather than applying absolute fixed-point rounding to tiny normalized DP entries.

A further q1 source audit found that the complete replay starts each shell
at the trivial probability bound 1. Its saved choices retain that branch
for shells 0..7 in both geometries. In particular, Cv for shells 5..7 are
nonzero and must be charged. The final q1 assembly therefore selects a
trivial or Chernoff bound inside each fixed support expectation, before
summing. `SparseQ1Choice22` drafts this original-setup path; its scalar
support count uses exactly choose(64,v). The final bound retains a separate
current-Cv charge when the trivial branch is selected. The earlier pure
Chernoff-family theorem remains valid but is not by itself the numerical
q1 endpoint. The diagnosed probability budgets have not changed.
- `PhysicalProfileFamilyExtract22` is kernel checked (core-only): the whole-family Bool result exposes each of the 17 canonical cached row checks, and each row exposes its sign, first/second moment, Fourier lower/upper, and return comparison premises. The extraction normalizes left-associated Boolean conjunctions explicitly.

### Routing lane: authenticated small birth endpoint budgets (2026-11-08)
Audited Numeric/FixedDefs and Numeric/Fixed outward arithmetic before adding Horner budget evaluation. PhysicalSmallBirth22BudgetDefs and BudgetCheck are pure checked-data modules; Budget proves Horner interval membership and Boolean-check-to-actual-birth-cone adapters. BudgetSmoke is GREEN for retained rational z=8607079764250578/10^16 and denominator 2^80: exact j0 initial envelope plus authenticated j1/j2 eleven-coordinate cones. All printed endpoints use only propext, Classical.choice, Quot.sound. No new numeric experiment or witness adjustment.
Before whole-family birth assembly, audited PhysicalCoarseBirth22_parts and the existing PhysicalProfileRows22.turboExactBirthMass22_coarse_cone: they partition one actual measure into its zero atom and arbitrary-live remainder; the total bound is the actual q0 input moment. New PhysicalProfileBirthFamily22 only combines that audited constructor with checked first/Fourier family bounds and replaces j0/1/2 by supplied exact cones. Compilation pending dependencies.
- `PhysicalProfileFamily22` now replays with exit 0 and standard axioms: cached first/second cap checks imply actual class maxima or class-2 uniform mean bounds; the integer return comparison implies min(M,sqrt(S F))/4194303 bounded by the supplied final return entry. Only ordinary left-associated Bool projection fixes were needed.

### Next numerical transfer audit: finite scaled powers
Read Numeric/FixedLog (Horner log-series, exact power-of-two reduction, positive endpoint checks) and Structured/DenseFourierExactCheck (exact Int cross-products mapped through positive denominators). Root clarified dense replay uses the finite 11-by-11 power M^8192, not a limiting spectral bound. Planned evaluator preserves that exact power via thirteen normalized squarings, keeping normalized interval matrices materialized and tracking rational scale logarithms separately. Normalization must support multiplication as well as division to prevent small positive powers underflowing the fixed denominator. No generator or numerical run started.
- `PhysicalProfileFamilyFourier22` now exits 0 with standard axioms. Its cached signed coefficient sum equals the checked uncached evaluator and gives actual physical zeroHold upper/lower bounds. Thus every retained family Bool field now has a semantic adapter; parent Assembly and routing birth-family closure can consume these without numerical regeneration.

### 2026-11-08 sparse fold checkpoint — SparseShellFold22 and dependents GREEN
Audit before repair: the failing proofs only rearranged finite sums; no
strategy change was needed. All failures were elaboration engineering:
- Any term-mode unification against a membership in the concrete product
  `range 65 ×ˢ range 65` makes the unifier evaluate `List.range 65` and hit
  the default recursion limit (even `h x` for an identical type). The
  two-coordinate layering is now proved over abstract `range n`
  (`two_weighted_layers`) and instantiated at n=65.
- Terminal-mass identities are proved once over abstract index sets
  (`mass_sum`, `mass_smul`, `mass_weighted_eq`, `mass_weighted_eq₂`); no
  `simp` with `Finset.sum_comm` runs over `Finset (Fin 64)` sums.
- The per-shell choice uses explicit weight definitions and `calc`, not
  `simp_rw` over `List.prod` matrix applications.
- `SparseActualShell22`: `simp_rw [← mul_div_assoc]` had reassociated the
  two-group weight (a/C₁)*(b/C₂); only the outer Chernoff division is
  reassociated now.
- `SparseMarkerFold22`: `rw [Finset.sum_mul_sum]` rewrote the wrong side
  through the let-bound total; the reverse rewrite closes it.
- `SparseQ1Choice22`: `?w #T` is not a higher-order pattern, so the choice
  fold is instantiated explicitly.
Complete-module exit 0 with only propext, Classical.choice, Quot.sound:
`SparseShellFold22`, `SparseActualShell22`, `SparseMarkerFold22`,
`SparseOccupation22`, `SparseQ1Choice22`. No resource limit was raised. No
statement changed; shells 0..7 keep their trivial-branch charge
rs8ExpectedPacketCount(v) = (Cv/choose(64,v))·choose(64,v).
Still unchecked drafts: `SparseFixedDefs22`, `SparseFixed22`,
`SparseSpectrumInteger22`, `MatrixScaledDefs22`, `MatrixScaled22`.

### 2026-11-08 loop tick 1 — remaining sparse drafts GREEN; symbolic q1/q2 closed
- `SparseSpectrumInteger22` GREEN (renamed `eq_div_iff₀` → `eq_div_iff`):
  rs8ExpectedPacketCount v = coeff_v(integer polynomial)/(2^32-1)^8.
- `SparseFixedDefs22` GREEN (core-only). `SparseFixed22` GREEN after adding
  `Mathlib.Algebra.BigOperators.Field` and `open Matrix`: integer matrix
  product upper check, positive Collatz check, and exact finite-power
  Collatz total bound are sound.
- Audit: the symbolic q2 assembly already exists in `SparseOccupation22`
  (`sparseTwoShellBound22_le_common_marker`,
  `turboSetup22_occupation_two_le_common_marker`), matching the diagnostic
  (K18 tilt 1/625, K20 tilt 1/2500, marker t=24, Collatz on the 11x11
  regional marker, power 64). The top-level identity
  `turboSetup22_expected_bad_count_by_occupation` plus
  `turboSetup22_bad_prob_le_mean` provide the occupation split.
- Soundness lemmas still missing in `SparseFixed22` for checks already
  defined in `SparseFixedDefs22`: polynomial product, region normalization,
  q1 shell-step/initial, and the 512-bit lower power trace. These are the
  next symbolic items, then the K18 numerical families.

### 2026-11-08 loop tick 2 — SparseFixedSound22 GREEN
Audit: same split as Numeric.Fixed / DenseFourierExactCheck (Mathlib-free
integer checks, positive-denominator soundness). New `SparseFixedSound22`
(exit 0, standard axioms) proves soundness of every remaining
`SparseFixedDefs22` check:
- `polyMul_sound`: degree-≤2 product check bounds Σ_{antidiagonal d}
  F a * G b, the exact `orderedMatrixDP_add` convolution.
- `normalize_sound`: `sparseRegionDenominator22 R d = choose(R*16,d)`
  (d ≤ 2, R > 0), so the check bounds `(1/choose) • X d`, i.e. the shape of
  `sparseRegionalDP22`.
- `q1Rows_sound`/`q1Rows_terminal`: right-append recurrence
  DP(n+1,v+1) = DP(n,v+1)W₀ + DP(n,v)W₁ derived from `orderedMatrixDP_add`
  (no step reordering); bounds Σ_j orderedMatrixCoefficient … N v 0 j by the
  materialized rows for all v < 65.
- `powerCertificate_sound`: a passing trace gives
  0 < last/scale ≤ (base/scale)^exponent (lower bound for z^d).
Next symbolic item: a binary-powering chain theorem that composes
`polyMul_sound` stages (materialized per stage, one-stage check plus
induction, not a recursive kernel reduction) into a bound on
orderedMatrixDP(Fin.val, C(16,j)•B j, R, d) for d ≤ 2, then
`normalize_sound` into `sparseRegionalDP22` entry bounds. Then the K18
numerical producer review/run.

### 2026-11-08 loop tick 3 — regional power chain GREEN; K18 families generated
- New core-only `SparseRegionalPowerDefs22` (GREEN) and
  `SparseRegionalPower22` (GREEN, standard axioms). MSB-first trace with an
  implicit leading one; each step stores its squared stage and, when the
  bit is set, the stage times the base. One `polyMul` check per stored
  stage; `chain_sound` is an induction over the trace (never a kernel
  reduction of a recursive power). `certificate_sound` adds the exact
  exponent check foldl = R and `normalize_sound`, giving entry bounds on
  `sparseRegionalDP22 B R d` for d ≤ 2. `base_of_rows` reduces the base
  stage to bounds on choose(16,d)·B_d.
- Reviewed `emit_turbo_profile_k18_families22.py`: untrusted producer,
  fail-closed precision gates (relative endpoint rounding < 2^-48), exact
  0/1 boundaries, class-2 uniform mean. Ran it once (single process):
  all three tilts (λ = 1/1250, 3/2500, 1/625; z = round(exp(-λ)·2^80)/2^80)
  pass every gate. z numerators: 1207959065712058707140553,
  1207475978709615532308063, 1206993084903331526962826.
- Generated per-tilt Lean modules by namespace/module substitution from the
  verified Smoke chain (Data, Check, BudgetCheck, BudgetSmoke,
  BirthFamily, Family). Kernel checks queued serially; results pending.
- RESULT: all 18 per-tilt modules exit 0 (each ≤ 45 s). For k = 1,2,3,
  `Spin.Turbo.ProfileFamilyK18Tilt{k}22.physicalRows_sound`,
  `selectedRows_sound`, both nonneg theorems, and `birthSmokeFamily22_sound`
  use only propext/Classical.choice/Quot.sound; `familyChecked` uses only
  propext. The three retained K18 local families are authenticated.
- Next: integer adapters for the rows actually used by the sparse kernel.
  Only physical rows j = 0,1,2 enter the occupancy-≤2 regional coefficient,
  so `base_of_rows` needs integer upper bounds on choose(16,d)·physicalRows d
  for d ≤ 2 (from the family data and birth cones), then a producer for the
  regional trace, q1 rows (65 shells × 64 regions), q2 marker Collatz data,
  and the power-trace endpoints; each certificate gets its own Data/Check
  module and the semantic composition uses SparseRegionalPower22,
  SparseFixedSound22, SparseFixed22, SparseQ1Choice22, SparseOccupation22.

### 2026-11-08 loop tick 4 — row adapter GREEN; Array checks infeasible, switch to Lists
- `SparseRowBaseDefs22` / `SparseRowBase22` GREEN: physical rows j ≤ 2 are
  exact integers over endpointDen (row 0 = exact small-birth list; rows
  i ≥ 1 carry return cap at k = 0 and first cap at k = 2), and
  `base_check_sound` gives choose(16,d)·B_d ≤ base/(den·mult).
- Diagnostic audit (q1_coarse_complete65): K18 L = 2048 groups (R = 128),
  d = 52428; three tilts reproduce the five-tilt q1 (58.2456 bits); chosen
  per shell: trivial 0..7, tilt1 8..13, tilt2 14, tilt3 ≥ 15. The whole
  K18 bound has slack factor only 1.000576 against 2^-58.24, so q1 must
  match the replay closely; shells ≥ 19 contribute < 1e-10 of q1.
- Precision: q1/regional scale 2^256 suffices (tilt3 amplification
  λd ≈ 121 bits; error budget ≪ 2^-71). Base mult = 2^176.
- Kernel feasibility probe (single job): one 65-shell q1 step with the
  flat-`Array` checker exhausts kernel memory (both `#[..]` and `Array.mk`
  literals) — indexed access is not "direct" in the kernel. The same step
  with list traversal (rows as List (List Int), vector×matrix by zipping
  rows, lockstep (next, prev, shifted prev)) passes in ≈ 12 s. Audit: the
  previously verified familyChecked data also used Lists.
  Decision: all large-data checks are re-expressed over Lists
  (new SparseListDefs22 / SparseListSound22). The Array checks stay sound
  but are not used for endpoint data.
- GREEN (standard axioms): `SparseListDefs22` (core-only),
  `SparseListSound22` (vecMat/matMul/matAdd specs, list polyMul and
  normalization soundness), `SparseListQ1Regional22`
  (`certificate_sound`: base check from exact row integers + regional
  chain + normalization ⇒ bounds on `sparseRegionalDP22 B R d`, d ≤ 2;
  `q1Rows_sound`: list q1 rows bound orderedMatrixDP for all 65 shells).
- Next tick: K18 certificate producer (per tilt: base at scale 2^256 with
  mult 2^176, MSB-first trace for R = 128, normalized result, 64 q1 row
  stages; one Data module and one or more `decide +kernel` Check modules
  per tilt, splitting q1 steps across theorems), then the q1 assembly
  theorem: occupation-1 ≤ 2048·Σ_v [trivial: rs8 count | Chernoff:
  (count/choose)/z^d·rows], with z^d lower bounds from the power trace
  and rs8 counts from `SparseSpectrumInteger22`. q2: list Collatz check
  on the regional marker matrix (t = 24, tilt 3).

### 2026-11-08 audit (user request): how the old proof handled kernel memory
Read `scripts/README.md`, `emit_scalable_dense_packed_cover.py`,
`emit_scalable_dense_packed_chunk.py`, `ScalableDensePackedLeaf.lean`.
The old Scalable dense cover is a hybrid: an untrusted Python candidate
(JSON, ~22 MB) is partitioned by an untrusted emitter into many small
generated Lean modules (611 parts, 39 segments, 46,317 leaves). Each part
proves one small `decide +kernel` theorem per item (`leaf{i}_checked`,
`witness{i}_checked`), closes the part's `List.all` by
`simp [refs, List.all, leaf_names]` from those per-item theorems (the kernel
never evaluates one giant Bool), and exposes a sound claim. Segment and
assembly modules combine parts with exact split lemmas. Sources/objects
live in ignored `.lake/fullcover/sources` with a manifest; a sequential,
resumable PowerShell driver compiles one module at a time. No Nat-packing
trick: the memory fix is bounded per-decision size plus modular assembly.
Reused for Turbo: untrusted producer → ignored generated sources + manifest
→ per-item small kernel decisions (one per regional stage, one per q1 row
step, one per power trace) → exact assembly theorem (`hstep` by
`interval_cases` from per-step theorems) → sequential resumable builder.
Different: our items are integer matrix/vector inequalities checked by
list traversal (tick-4 probe), not rational interval leaves.

### 2026-11-08 loop tick 5 — K18 certificate producer; generated modules building
- `scripts/emit_turbo_k18_sparse_cert22.py` (untrusted) emits, per tilt,
  ignored generated modules under `SpinCodes/Turbo/Generated/K18/`
  (local `.gitignore`): `Tilt{k}Regional` (row function from the authenticated
  family data, base, 7 squaring stages for R = 128, normalized result,
  z^d power trace; theorems `regionalChecked`, `powerChecked`),
  `Tilt{k}Rows{c}` (q1 rows, 8 stages per module), `Tilt{k}Steps{c}`
  (one `decide +kernel` theorem per q1 step). Sequential resumable builder
  `.lake/turbo/build_k18_cert.ps1`. Manifest in the same directory.
- The producer asserts the integer RS8 counts equal the receipt's exact
  expected shells for all 65 v.
- Precision correction: scale 2^256 was insufficient — absolute row
  rounding, amplified by 1/z^d ≈ 2^121 and counts up to C_64 ≈ 5·10^36,
  inflated shells ≥ 22 (shell 64 contributed ≈ 2). Error ≈ 2^(264-S); now
  S = 2^384 (mult 2^304). Producer's exact rational q1 = 2^-58.24563386865,
  matching the diagnostic (58.245633868649) to 13 digits; all chosen shells
  8..21 match the diagnostic per shell.
- Tilt 1: `regionalChecked` (propext only, ~8 min kernel) and `powerChecked`
  (no axioms) GREEN; first step module (8 steps) checked in 10 s. Full
  three-tilt build running in background.
- Next: Mathlib assembly per tilt (`rows` by match, `hstep` by
  `interval_cases` from step theorems, `certificate_sound` with
  `profileFamilyRows22_small_eq`, `q1Rows_sound`, `powerCertificate_sound`),
  then the q1 sum theorem with exact rational comparison.
- RESULT: all 54 generated certificate modules exit 0
  (`regionalChecked` propext only; `powerChecked` no axioms; 192 q1 step
  theorems). `scripts/emit_turbo_k18_q1_assembly22.py` emits
  `Generated/K18/Tilt{k}Q1.lean` (k = 1,2,3), all exit 0 with standard
  axioms (~2 min each): `rows` by match, `hstep` by `interval_cases`,
  `regional` (certificate_sound + exact row integers), `q1`, and
  `shell_le v`: Σ_j orderedMatrixCoefficient … (sparseOneRegionChoice22
  physicalRows 128) 64 v 0 j ≤ (Σ_j rows64[v][j])/s for all v < 65, plus
  `zpow_lower`: 0 < last/s ≤ z^52428 with the family's exact rational z.
- Next: K18 q1 endpoint theorem. Instantiate
  `turboSetup22_occupation_one_le_chosen_shells` (R = 128, d = 52428) with
  z v / B v from the chosen tilt (trivial 0..7, tilt1 8..13, tilt2 14,
  tilt3 ≥ 15), bound each Chernoff shell by
  count_v/choose(64,v) · (s/last_k) · (Σ rows64[v])/s, counts from
  `sparseRS8ExpectedPacketCount22_eq_integer_coeff`, then one exact
  rational comparison against the producer's q1 value.

### 2026-11-08 loop tick 7 — exact RS8 counts GREEN
- `SparseRS8Counts22` (standard axioms): computable list mirror of integer
  polynomial arithmetic (`lAdd`, `lMul`, `lPow`) with a proved
  homomorphism `toPoly`; `countList_eq : toPoly countList =
  sparseRS8IntegerPolynomial22`; one kernel decision fixes the 65
  coefficients (`countList_lit`); `count_eq`: rs8ExpectedPacketCount
  points v = countLit[v] / (2^32-1)^8 for injective points.
- Next: `Generated/K18/Q1Endpoint.lean`: instantiate
  `turboSetup22_occupation_one_le_chosen_shells` (R = 128, d = 52428, z and
  B by the chosen tilt, useChernoff v := 8 ≤ v), bound each shell termwise
  with `shell_le`, `zpow_lower`, `count_eq`, and close with integer per-shell
  upper numerators a_v over 2^K (producer-supplied, kernel-checked by
  cross-multiplication) and one integer comparison of 2048·Σ a_v.
- **K18 q1 ENDPOINT CLOSED.** `Generated/K18/Q1Endpoint.lean` (template
  `scripts/templates/k18_q1_endpoint.lean.tmpl`, emitter
  `scripts/emit_turbo_k18_q1_endpoint22.py`): `Spin.Turbo.K18Q1.q1_le`:
  for every injective `points`, the actual sampled-setup occupation-1
  expectation (R = 128, d = 52428) ≤ 2048 · aSum / 2^256 = 2^-58.24563386865;
  axioms propext, Classical.choice, Quot.sound. Per-shell integer checks
  (one `decide +kernel` per tilt), binomials via factorials (kernel
  `Nat.choose` recursion is exponential), counts via `SparseRS8Counts22`.
- **K18 q2 ENDPOINT CLOSED.** `Generated/K18/Q2Endpoint.lean` (template
  `scripts/templates/k18_q2_endpoint.lean.tmpl`, emitter
  `scripts/emit_turbo_k18_q2_endpoint22.py`): `Spin.Turbo.K18Q2.q2_le`:
  occupation-2 expectation ≤ b2/2^256 = 2^-104.18 (diagnostic 105.37; the
  difference is Collatz rounding). Marker t = 24 on the authenticated tilt-3
  regional kernels; Mint dominance, Collatz vector, exact outer weight sum
  W = Wn/(M8·64!·24^64), z^d lower trace, one final integer check. Axioms
  standard. (`positivity` hit recursion limits on 64th powers; explicit
  terms used.)
- K18 budget now: target 2^-58.24; q1 = 2^-58.24563, q2 = 2^-104.18 ⇒
  the q ≥ 3 tail must be ≤ ≈ 2^-66.25 (diagnostic tail 2^-66.474).
  Next: audit the tail replay and the old proof before building it.

### 2026-11-08 K18 q ≥ 3 tail — audit and design
Audited: `scripts/replay_turbo_coarse_k18.py`, `workstreams/k18_margin60/
replay_fractional.py`, Lean `RoutingFractionalTail22`/`Endpoint22` (GREEN,
untracked), `MatrixMarker22.pairedRegionDP22_marker_pow64_le` (GREEN).
- Symbolic interface exists: `turboActiveBadCount22_prob_union_le_paired_dp`
  (setup R = 2n; K18 n = 64) bounds P(∃ q ∈ qs, active set of size q bad) by
  Σ_q C(32n,q)·(β^q / z_q^d)^{α_q}·Σ_k (pairedRegionDP22 B_q α_q n q ^ 64)_0k,
  β = 2^256/65535^8, rows B_q sound for the *selected* kernel (family
  `selectedRows_sound`), and `turboOccupationTail22_prob_le_active_union`
  covers the actual-message tail with qs = Ico 3 (32n+1).
- Replay witnesses: q = 3..384 exact regional (13 families; α = 4/5 at
  tilt 17/200 for q = 62..143), q ≥ 385 "iid conditioning" with marker
  probability p and M^8192 — this dense formula has NO Lean theorem yet.
  Tightest per-q margins: 67.40 (q 144..157), 67.75 (52..61), 70.71 (35..51).
  Tail budget: ≤ ≈ 2^-66.25 against replay 2^-66.474.
- Design (to keep the kernel work bounded and avoid new theory where the
  margin allows):
  (a) q ≥ 225 (and the dense range): use the existing union theorem with
      the GREEN marker bound `pairedRegionDP22_marker_pow64_le` plus a
      Collatz/scaled power for M^(64n); margins there are ≥ 470 bits, so
      the per-region Chernoff loss is affordable. Must be confirmed
      numerically before building (one bounded diagnostic run).
  (b) q = 3..224: exact coefficient q of the paired profile polynomial
      power per family (degree ≤ family q_max), binary powering with
      materialized stages and per-stage binary exponents (entries
      underflow any absolute grid), then the 64-region matrix power with
      exponent scaling, and the scalar prefactor via exact integer checks
      (α = 4/5 handled by fifth-power checks y^5 ≥ x^4).
  (c) z^d for high tilts: exponent-scaled power trace (mantissa, exponent),
      not the absolute 2^384 grid.
  Families needed: selected-row local families at each retained tilt
  (producer gates may fail at high tilts; examine before running).
- Numerical check (single bounded run; per-q optimum over 16 tilts × 999
  values of p, t = p/(1-p)): the replay's dense "iid conditioning" formula
  equals the marker bound `pairedRegionDP22_marker_pow64_le` with
  t = p/(1-p) (verified to 2^-996.6289 at q = 385). So the dense range needs
  NO new theory. Marker loses ≈ 500–700 bits vs exact coefficients for
  q ≤ 384 (≈ 9 bits/region × 64), but still gives ≥ 340 bits at q = 280..384
  and ≥ replay for q ≥ 385 (q → 2048 needs t → ∞; use very large t).
  Decision: exact coefficient families for q = 3..279 (11 tilts: 3/1000,
  3/500, 3/250, 3/125, 1/25, 11/200, 17/200 (α = 4/5, q 62..143), 7/50,
  4/25, 1/5, 6/25); marker + scaled power for q = 280..2048.
- Low-tilt substitution rejected (bounded replay): tilt 3/125 is negative
  from q ≈ 39 (−1677 bits at q = 63); 3/250 negative from q ≈ 27. The tail
  needs the high tilts.
- Precision finding: the family checker evaluates in `Spin.Numeric.Fix`
  (absolute scale 10^30 ≈ 2^-100). High-tilt row quantities reach 2^-207
  (z^128 at tilt 1), so producer gates fail from tilt 1/25 upward with the
  2^-80 grid, and an absolute 2^-100 rounding would inflate them by up to
  2^100 — unaffordable for the dominant tail families (tilts 11/200, 7/50).
  The dependency `Numeric.FixedDefs` is shared with the Scalable proof and
  must not change. Decision: port.
- `scripts/port_fix_hi22.py` copies `Numeric.FixedDefs`/`Fixed` to
  `SpinCodes/Turbo/Hi/` (namespace `Spin.NumericHi`, prec 120, scale 10^120;
  `Fixed.lean` is generic in `scale`, log soundness lives elsewhere) and the
  21-module family chain to `SpinCodes/Turbo/Hi/` (namespace
  `Spin.Turbo.Hi`; enclosing-namespace resolution keeps base references).
  Build order `.lake/turbo/build_hi.ps1`; build running.
- `scripts/emit_turbo_tail_families22.py`: 11 tail tilts (3/1000 … 6/25),
  endpoint grid 2^-384 (same ≈ 2^14 ratio over Fix rounding as before),
  all precision gates pass; six per-tilt modules each in ignored
  `SpinCodes/Turbo/Generated/K18Tail/` (from the verified Smoke chain).
- Also GREEN this tick: `SparseListSquare22` (scaled squaring chain with
  integer exponent per stage, `chain_sound`).
- Subagent result GREEN: `SparseRpowCert22` (core list check
  `sparseRpow45CheckL22 sx sy X Y`: Y ≥ 0 and X^4·sy^5 ≤ Y^5·sx^4 entrywise;
  `sparseRpow45CheckL22_sound`, `..._entrywiseRpow`, `rpow45_scalar_le`,
  `rpow45_mono`, `rpow45_le_of_pow`; standard axioms). Apply with
  M := B_a * B_b, then scale by C(16,a)·C(16,b) for `pairedProfileMatrix22`.
- GREEN: high-precision port, all 23 `SpinCodes/Turbo/Hi/*` modules (no
  errors, no sorryAx). `Hi/RowSelected22` (`rows_eq`: every physical row
  entry, incl. coarse rows j ≥ 3, = integer/den; `selected_eq`: selected row
  d = rowSel/(16^16·den) with a Pascal table pinned by one kernel decision).
- Subagent result GREEN: `SparseScaledPower22` (`scaledPower_sound`:
  m·2^(-E) ≤ (zNum/2^zBits)^d with ~600-bit mantissas; kernel examples incl.
  λ ≈ 1, d = 209715 in 18 steps).
- GREEN: `PairedSquare22.pairedRegionDP22_one_eq`: pairedRegionDP22 B 1 n q =
  (1/C(32n,q)) • regionalCountDP22 B (2n) q (α = 1 paired polynomial is the
  square of the single-step polynomial). Tail-family modules (66) building:
  `.lake/turbo/build_k18_tail_families.ps1`.
- Subagent result GREEN: `SparseListPolyD22` (degree-D list polynomial
  product, `polyMulD_sound`, `PolyBoundD`, `chainD_sound`, `baseD_sound`).
- `TailTerm22.tail_term_le` (hand-written generic per-q bound, α = 1):
  chain bound → `pairedRegionDP22_one_eq` → U0 domination → 6 scaled
  squarings → scaled z^d → one integer check (all powers of two folded into
  T = 256q + Ez + E + Kt − S). Compile queued.
- Exact tail producer `scripts/emit_turbo_k18_tail_exact22.py` (S = 704,
  Kt = 256, flint products) and semantic emitter
  `scripts/emit_turbo_k18_tail_term22.py`. Families 1..6 (q 3..61) generated
  (42 MB); every reported q matches the replay to 3 decimals (e.g. q 52:
  67.755, q 51: 70.713; values printed as 256.000 are the 2^-256 floor).
  Builder `.lake/turbo/build_k18_tail_exact.ps1` running; family 8
  (7/50, q 144..157) generating; families 9..11 next. Subagents: α = 4/5
  paired chain/term (family 7, q 62..143) and marker-range term (q ≥ 280).
- Families 8, 9 match the replay (q 144: 67.397; q 158: 90.482). Family 10
  (q ≥ ~210) and 11 collapsed (q 225: −3136 bits): U0 = ceil(last_q/C(2048,q))
  on the absolute 2^-704 grid rounds tiny normalized coefficients up to 1.
  Fix: domination with its own binary shift, 2^down0·last_q ≤ C·2^up0·U0,
  chain starting at E0 = up0 − down0 (`TailTerm22` patched; producer picks
  the shift so max U0 ≈ 2^(S+8); stage cache `.lake/turbo/tail_stages_k.json`).
  All α = 1 families regenerating; exact builder restarted afterwards.
- Subagent result GREEN: `TailMarker22.tail_marker_le` (q ≥ 280 marker term,
  α = 1; marker identity A_t·A_t, 13 scaled squarings, t = (tn0·2^kp)/(td0·2^kn)
  with tn0, td0 ≤ 2^12; final check stress-tested at q = 2048 (~1.6·10^6-bit
  integers)). Tilt set for q ≥ 280 being chosen by a bounded numerical run
  (`.lake/turbo/marker_cover.json`).
- Marker range tilt set (bounded run): segments 6/25 (q 280–281, family 11,
  81/76 bits), 8/25 (282–344), 2/5 (345–427), 12/25 (428–546), 16/25
  (547–769), 1 (770–1499), 7/4 (1500–2048, ≥ 8584 bits even at q = 2048;
  the replay used 3/2 and 11/5 there, but z^128 at 11/5 ≈ 2^-406 is below
  the 10^-120 Fix). New Hi families 12–17 generated (7/4 trips only the
  diagnostic gate; marker-only, huge margin) and building.
- `scripts/emit_turbo_k18_tail_marker22.py`: shared (family, t) chains on
  the grid t = 2^(j/8) (base check + 13-stage chain once per t), per-q final
  checks via `tail_marker_le`; modules Tail{k}MkRows / Mk{tag} / MkTerm.
- Subagent result GREEN: `PairedFracChain22` (generic step-family chain
  `chainW_sound`, rpow pair check `pairedRpowCheckL22`, base check,
  `pairedChain_sound`) and `TailTermFrac22.tail_term_frac_le` (α = 4/5 per-q
  bound; prefactor via 65535^(8q) and z^d scaled traces with G/yn split).
  Note: its domination step has no binary shift; the family-7 producer
  (`scripts/emit_turbo_k18_tail_frac22.py`, running) reports per q whether
  that underflows (q ≤ 143 expected fine).
- Marker producer fixed (2^(1/64) grid for families 11–15, wide t scan for
  family 17, exact per-q choice): all q = 280..2048 at the 2^-256 floor.
- Producer-side K18 total from manifests: q1 2^-58.245634, q2 2^-104.180,
  tail (2046 occupations) 2^-66.4739, TOTAL 2^-58.240831 (< 2^-58.24;
  replay 58.2408311085). Lean confirmation pending builds + final assembly.
- Final K18 assembly design: `BadSplit22.turboSetup22_bad_prob_split`
  (written, compile queued): P(Z ≥ 1) ≤ E[ZQ_1] + E[ZQ_2] + P(∃ msg with
  occupation in Ico 3 (16R+1) bad). Then q1_le, q2_le, and for the tail
  `turboOccupationTail22_prob_le_active_union` +
  `turboActiveBadCount22_prob_union_le_paired_dp` (n = 64, qs = Ico 3 2049,
  piecewise z/α/B by family) closed termwise by famTerm/frTerm/mkTerm, and
  one integer check N^25·2^1456 ≤ 2^6400 for 2^-58.24.
- Builders: tail families (`build_k18_tail_families.ps1` ×2), exact
  families (`build_k18_tail_exact.ps1`), family 7 (`build_list.ps1`), marker
  families 11–17 (detached `build_list.ps1`, logs `.lake/turbo/mk*_build.log`).
- 2026-11-09: exact builder bug (PowerShell function output swallowed the
  status, so failures did not stop it) fixed with a script-scoped flag.
  `TailTerm22` had failed (ambiguous `zpow_neg` vs `Matrix.zpow_neg`,
  `positivity` recursion on z^52428, fragile exponent identity): fixed,
  recompiling. Throughput: memory pressure (12.8 GB committed vs ~4 GB
  resident) made `Tail1ExChain` take ~9 h. Generated chain/marker modules
  patched in place to use a literal row table `rowL`/`rowT` (one kernel
  equality `rowT_eq` per family) instead of evaluating the Finset-based
  `rowSel` inside every base check. TODO: fold the rowL patch into
  `emit_turbo_k18_tail_exact22.py` / `emit_turbo_k18_tail_marker22.py`
  emitters before any regeneration.
- GREEN: `TailTerm22.tail_term_le` and `BadSplit22.turboSetup22_bad_prob_split`
  (standard axioms). Exact builder relaunched (now fails fast).
- MILESTONE: `Tail1ExTerm.famTerm` GREEN (standard axioms) — first tail family
  closed end to end (authenticated high-precision family → selected rows →
  rowT base → degree-D chain → paired coefficient → region power → z^d →
  integer check). Generated-file tactic fix: hrow uses `unfold row; norm_num`.
  Exact builder relaunched for families 2–6, 8–11.
- 2026-11-09 Peach backend (user-approved): elan + Lean 4.34.0 on
  peach48 (16 cores, 124 GB), Mathlib/packages via `lake exe cache get` at
  the pinned manifest, dependency + worktree oleans uploaded (oleans are
  platform-independent for one Lean version), Windows symlinks in the build
  tree repointed (`fixlinks.sh`). Driver `~/turbo/pbuild.py` builds
  dependency-ordered, 6 parallel jobs, per-module logs in `~/turbo/logs`.
  Local builders stopped.
- Scoped exception (audited): pure data modules (ExBase/ExSq*/ExQData,
  FrR/FrY*/FrBase/FrSq*/FrQData) carry `set_option maxRecDepth 10000`
  because elaborating nested list literals of degree up to 280 exceeds the
  default depth; no proof search runs in these modules (same precedent as
  the audited BudgetCheck modules). MkTerm statements now cast numerators
  as ((N : ℤ) : ℝ) to match `tail_marker_le`.
- Peach pass findings (fixed in generated files; TODO fold into emitters):
  (1) large single list literals exceed recursion depth even at 10000 —
  data defs split into one def per matrix (`name_partI`) plus the list of
  names (kernel reduction unchanged); (2) the 730-case `interval_cases` in
  `mkTerm` exceeds recursion depth — chunked aggregators `mkTermC*`
  (≤ 100 q each) combined by `by_cases`/`omega`.
- 2026-11-09 K18 assembly checkpoint:
  * Earlier Peach failures on `Tail*Sq*` were a race with the patch upload
    (the pass read pre-split sources); rebuilt cleanly, no limit changes.
  * Marker aggregators `mkTerm` failed on late cases: `aMk` was a 730-arm
    literal `match`, so `exact mterm_q` unfolded deep `Nat.decEq` chains.
    All marker bounds are exactly 1, so `aMk` is now the constant 1
    (emitter `emit_turbo_k18_tail_marker22.py` emits the constant whenever
    every bound is 1).
  * `SpinCodes/Turbo/EndpointBits22.lean` GREEN (standard axioms):
    `int_div_two_pow_le_two_rpow_neg`: T ≥ 0, T^n·2^m ≤ 2^(s·n) ⇒
    T/2^s ≤ 2^(-(m/n)). K18 uses s=256, m/n=1456/25=58.24; K20 will use
    2591/50=51.82.
  * `scripts/emit_turbo_k18_endpoint22.py` → `Generated/K18/Endpoint.lean`:
    piecewise zq/αq/Bq/A over 18 pieces (17 families, family 11 split
    Ex 225–279 / Mk 280–281), selector lemmas by `simp only` on breakpoint
    facts, one dispatcher `piece_cases`, tail union with n=64, and a single
    kernel check `(2048·aSum + b2 + Σ_{q∈[3,2048]} A q)^25·2^1456 ≤ 2^6400`.
    Final statement: `K18Final.k18_endpoint`: for injective points,
    P[(turboSetup22 128 points) has a nonzero message of weight ≤ 52428]
    ≤ 2^(-58.24). Pending: build after all term modules are green.
  * GREEN on Peach (standard axioms): `Tail{11..17}MkTerm.mkTerm` (all
    marker families, q = 280..2048) after the constant-`aMk` fix; also
    `Tail1ExTerm`. Remaining: ExChain/ExTerm families 2–6, 8–11 and
    Tail7FrTerm, then `Generated/K18/Endpoint.lean`.
- 2026-11-09 STRATEGY CHANGE (tail squaring chains → packed Kronecker checks).
  Audit: the interpreted degree-D check `sparseRegionalPowerCheckDL22`
  costs ~2 s kernel time per 11×11 matrix product and needs
  (D+1)(D+2)/2 products per squaring. Family 1 (D=5) passed in 301 s; family 2
  (D=10) hit kernel "excessive memory consumption" after 70 min; families up
  to D=279 (and fractional family 7, D=143) were infeasible. The old
  Structured proof never needed high-degree polynomial products, so there is
  no old component to reuse; what IS reused from the old proof is its pattern
  of small core-only kernel-check modules with Mathlib-level soundness proved
  once.
  New: `SparseListPackedDefs22` (core-only) + `SparseListPacked22` (GREEN,
  standard axioms). Each entry (i,k) of a degree-D polynomial matrix is
  packed as Σ_t P_t[i][k]·2^(w t). With all entries in [0,2^v) and
  11(2D+1)4^v < 2^w, the base-2^w digits of Σ_k pack(A_ik)·pack(B_kj) are
  exactly the convolution coefficients (no carries; proof via ℕ[X] evaluation
  at 2^w and digit extraction). The product is truncated mod 2^(w(D+1)) and its
  digits are walked against s·u_d. Soundness `chainPk_sound` yields the same
  `PolyBoundD` interface, so term modules change only in the chain call.
  Fractional family: `PairedFracPacked22` (`pairedChainPk_sound`) and
  `TailTermFracPk22` (`tail_term_frac_pk_le`), mechanical copies of the green
  proofs with the check swapped (GREEN locally).
  Kernel memory: Lean's `-M` counts the ~8 GB Mathlib import, so step modules
  import only the core defs; each squaring step is one module with a header
  declaration plus 11 per-row declarations. Measured: D=279 step ≈ 83 s,
  peak 8.9 GB; D=10 step ≈ 6 s. Limits unchanged (`-M 14000`).
  Generator: `scripts/emit_turbo_k18_tail_packed22.py` (v,w per family in
  `packed_manifest.json`).
- 2026-11-09 K18 ENDPOINT GREEN on Peach: `K18Final.k18_endpoint`
  (P[bad] ≤ 2^-58.24, standard axioms). Further fixes on the way:
  * One process accumulates kernel memory across declarations, so the
    packed steps are split into per-row modules (11 per-entry declarations
    each); family-7 `rpowChecked`/`baseChecked` split into 17 + 12 modules
    with generic assembly lemmas `rpowCheck_of`, `baseGo_append`
    (`PairedFracPacked22`, GREEN).
  * Endpoint selectors: `simp only`/`split_ifs`+`rfl` timed out (rfl unfolds
    real-valued z); replaced by deterministic `unfold; rw [if_neg …, if_pos …]`.
  * Closure audit: 1450 modules, sources identical to the local worktree
    modulo CRLF; oleans being pulled back for a local `#print axioms` audit.
  Next: K20 (R=512, cutoff 209715, 51.82 = 2591/50 bits) — producers for
  q1/q2/tail at K20 tilts, then the same packed/marker pipeline.
- 2026-11-09 K20 scoping (from `.lake/turbo/k20_coarse_full_diagnostic.json`,
  `k20_coarse_extended_diagnostic.json`, `q2_common_marker_k20_full_recombination.json`).
  Budget: rare ≈ 2^-51.9 (dominant), q1 2^-56.257, q2 2^-102.68, tail
  2^-100.594, whole 2^-51.82247 — only 0.0025 bits above 51.82, so every
  component must be certified close to its diagnostic value.
  Lean shape: `turboSetup22_bad_prob_le_rare_truncated` (rare charged once,
  GREEN) + truncated count ≤ q1 + q2 (full unconditional, as in K18) +
  q≥3 tail; R = 512 (n = 256 pairs), 8192 groups, cutoff 209715.
  New infrastructure required beyond K18:
  1. Exact rare budget `(R·16)·Σ_{v=5..8} rs8ExpectedPacketCount points v`
     as a checked rational (C_v closed form) — audit existing
     `rs8ExpectedPacketCount` lemmas first.
  2. Fractional exponents α ∈ {3/5, 7/10} (K18 had only 4/5): generalize
     `sparseRpow45CheckL22`/`PairedFracChain22`/`TailTermFrac22` to α = a/b.
  3. Exact family q∈[248,615] at α = 7/10, D = 615 (packed check; cost
     estimate ~5× the D=279 steps).
  4. q ≥ 616: `iid_fractional_conditioning` markers (fractional α with
     marker probability) — K18's `TailMarker22` covers α = 1 only; needs a
     fractional marker term theorem.
  5. Check whether K20 tail rows rely on the heavy (weight ≥ 9) truncation
     (`turboSparseHeavyMessage22`); if so, the tail union must run on the
     truncated count, not on `turboOccupationTail22_prob_le_active_union`.
  Generators to parameterize by geometry: emit_turbo_k18_{q1,q2}_endpoint22,
  emit_turbo_tail_families22, emit_turbo_k18_tail_{exact,frac,marker,term,
  packed}22, emit_turbo_k18_endpoint22.
- 2026-11-09 K18 CLOSED: local audit (`K18Audit.lean`) of the pulled oleans
  prints standard axioms only for `K18Final.k18_endpoint`.
- 2026-11-09 K20 DESIGN DECISION (audit): the K20 tail replay
  (`scripts/replay_turbo_coarse_k20.py`, q=3..615) and its extension
  (q=616..8192, `k20_small_outer_screen/METHODS.md` "dense components")
  use the same untruncated `CoarsePhysical` per-step model and the
  routing-conditioned fractional first moment, exactly as K18. K18's Lean
  tail families are sound for every message with no heavy (weight ≥ 9)
  hypothesis, and K18 closed without a rare charge. So the K20 Lean endpoint
  uses BadSplit22 (q1 + q2 + tail) WITHOUT the rare budget: strictly sound
  (fewer terms) and it raises the K20 slack from 0.0025 bits to ≈ 4.4 bits
  (q1 ≈ 2^-56.257 dominates; q2 2^-102.7, tail 2^-100.6). The rare-event
  modules (`SparseRare22`, `SparseTruncated22`) remain green but unused.
  Still required: general fractional exponent α = a/b (3/5, 7/10, 4/5) in
  the rpow certificate, paired chain and tail-term theorem; a fractional
  marker theorem for q ≥ 616 (K18's `TailMarker22` is α = 1); K20 q1 (4
  tilts) and q2 (marker 24, tilt 1/2500) producers; geometry R = 512
  (n = 256 two-step macros per region, 64 regions), 8192 groups,
  cutoff 209715, target 2591/50 bits.
- 2026-11-09 K20 infrastructure, GREEN locally (standard axioms):
  * `SparseRpowABDefs22` (core) / `SparseRpowABCert22`: entrywise
    x^(a/b) certificates `X^a·sy^b ≤ Y^b·sx^a`, any b ≠ 0.
  * `PairedFracABDefs22` (core; per-row `pairedRpowABRowL22` +
    `pairedRpowABCheck_of` for per-process splitting) / `PairedFracABChain22`:
    `pairedChainABPk_sound` (paired fractional DP at α = a/b, packed chain).
  * `TailTermFracAB22`: `tail_term_fracAB_le` — one exact tail term at any
    α = a/b and geometry (n pairs/region, 64 regions, C(32n, q), cutoff d,
    z = zNum/2^zb); α = 1 is a/b = 1/1, so all K20 exact families
    (q = 3..615, α ∈ {1, 3/5, 4/5, 7/10}) use this one theorem.
  Remaining for the K20 tail: fractional marker term (q ≥ 616): generalize
  `TailMarker22.tail_marker_le` (single-step α=1 marker, 13 squarings) to the
  paired marker `orderedMatrixMarker22 pairedProfileDegree22
  (pairedProfileMatrix22 B (a/b)) t` with pair certificates Y_ab, using
  `MatrixMarker22.pairedRegionDP22_marker_pow64_le` (already general in α, n),
  log2(64n) = 14 squarings at K20.
- 2026-11-09 `TailMarkerFracAB22.tail_marker_fracAB_le` GREEN locally
  (standard axioms): dense tail term at α = a/b and general geometry via the
  paired marker M_t = Σ t^(a+b) W(a,b), pair certificates Y_ab
  (`tailMarkerPairBaseCheck22`), `pairedRegionDP22_marker_pow64_le`,
  2^len = 64n squarings (`hlen`), `prefactorAB_le`, and one final check
  `yn·sumU·td0^(64q)·2^T ≤ a·C^63·tn0^(64q)`.
  All Lean theorems the K20 tail needs now exist. Next K20 steps (untrusted
  producers + generated check modules, built on Peach):
  1. K20 tail families: physical rows at the K20 tilts (8 exact pieces
     q=3..615 with α ∈ {1,3/5,4/5,7/10}; 15 dense marker pieces q=616..8192,
     tilts/markers from `k20_coarse_full_diagnostic.json` and the K20
     receipts), pair certificates for every α, packed chains (n = 256 → 8
     squarings to 256 pairs; markers 2^14 = 64·256), per-q term modules.
  2. K20 q1 endpoint (4 tilts) and q2 endpoint (marker 24, tilt 1/2500):
     port `emit_turbo_k18_q{1,2}_endpoint22.py` (R = 512, 8192 groups,
     d = 209715).
  3. K20 endpoint assembly: BadSplit22 at R = 512, tail union with n = 256,
     final check `T^50·2^2591 ≤ 2^(256·50)` (51.82 = 2591/50 bits).
- 2026-11-09 K20 progress:
  * q1/q2 GREEN on Peach (standard axioms): `K20Q1.q1_le` (4 tilts 1/5000,
    1/2500, 3/5000, 1/1250; scale raised to 2^704 because z^209715 ≈ 2^-242
    made 2^384 rows hit the precision floor; total 2^-56.25709 = diagnostic)
    and `K20Q2.q2_le` (marker 24 on tilt 1/2500, 2^-101.59). Producers:
    `emit_turbo_profile_k20_families22.py`, `emit_turbo_k20_sparse_cert22.py`,
    `emit_turbo_k20_q1_assembly22.py`, `emit_turbo_k20_q1_endpoint22.py`,
    `emit_turbo_k20_q2_endpoint22.py` (+ templates k20_q{1,2}_endpoint).
  * Tail producer `emit_turbo_k20_tail22.py` (13 families, 14 pair sets,
    8 exact + 15 marker pieces) and emitter `emit_turbo_k20_tail_lean22.py`.
    Exact pieces 1–7 reproduce the diagnostic per-q values exactly
    (q=3: 2^-100.594); marker pieces have ≥ 6000 bits of margin (a = 1).
    Bug fixed: marker base U0 was scaled to 2^8 instead of 2^(S+8).
    Piece 8 (q=248..615, D=615, α=7/10) still computing; build of the other
    22 pieces (1408 modules) running on Peach.
  * `emit_turbo_k20_endpoint22.py` → `Generated/K20/Endpoint.lean`
    (`K20Final.k20_endpoint`, ≤ 2^-51.82, check T^50·2^2591 ≤ 2^12800).
- 2026-11-09 K20 tail replan (audit of the diagnostic): the full K20
  diagnostic's exact witness for q=248..615 has only 20.36 bits of margin and
  fails the 51.82 target; the extended diagnostic re-covers 507..615 with a
  marker. Our degree-615 exact chain also hit a precision floor. New plan
  (`emit_turbo_k20_tail22.py plan` → `.lake/turbo/k20tail/plan.json`):
  exact 248..315 (tilt 17/200, α 7/10) + 9 greedy marker chunks over
  316..615 (pair sets 17/200@7/10 and 3/20@3/5, per-chunk t). Exact pieces
  now use `TailTermFrac22.tail_term_fracAB_shift_le` (GREEN; domination shift
  2^down0/2^up0 as in K18's `TailTerm22`) because DP_q/C(8192,q) ≈ 2^-650
  underflowed U0. Result: 32 pieces; exact pieces match the diagnostic
  per-q values; every marker term ≤ 2^-256 with ≥ 1883 bits to spare.
  Full K20 build (2714 modules incl. endpoint) running on Peach.
- 2026-11-09 K20 CLOSED: `K20Final.k20_endpoint` (≤ 2^-51.82, certified
  ≈ 2^-56.257) GREEN on Peach and in the local audit `K20Audit.lean`
  (standard axioms; closure sources identical modulo CRLF). Both K18 and
  K20 endpoints are closed and audited.
- 2026-11-09 Tightened endpoints: K18 ≤ 2^-58.24083, K20 ≤ 2^-56.25709
  (GREEN, locally audited). Paper updated to cite Lean as the authority.
