# Archived proof history through 2026-08-13

This file preserves the former append-only status ledger. It contains
superseded frontiers and next-step recommendations. Use `PROOF_STATUS.md` for
the current source of truth.

# Proof Status And High-Risk TODOs

Updated: 2026-08-12

Purpose: keep the post-restructure proof state explicit.  This is not a new
theorem source; it is the working checklist for what is in good narrative
shape, what is checkable, and what remains too risky to hide in prose.

## Current Narrative State

The broad paper restructure is done.  The compiled manuscript now has one
clean spine:

1. introduction and preliminaries;
2. modular first-moment framework;
3. accumulator warmup;
4. random sliding dense construction;
5. structured local-code construction.

The original exploratory branches remain preserved but uncompiled:
`innerDense.tex`, `integration.tex`, `outerExpandAcc.tex`, and
`innerSparse.tex`.  Topic notes under `explorations/` retain the dense finite
diagnostics, block-recursive BCH inner attempts, BCH spectra work, outer
spectrum comparisons, sparse inner notes, and expander/accumulator outer notes.

## Riffle Packet-Permutation Proof Frontier

Status: g=2 has a complete outward certificate; g=4 is an active diagnostic
proof search and is not yet an end-to-end theorem.

The Claude second-opinion audit motivated two changes that are now implemented:

- positive floating-point underflow in the packet-profile DP fails closed;
  structural zero transitions remain exact zeros;
- `explorations/riffle_group_chain_proof.md` now states the packet-profile
  transitivity lemma and the required independent lane-bijection invariant.

For g=4, exact rational per-cell BSP refinement is strongly effective.  The
original worst full-support cell `s1fr022089` fell from branch
`+67198.7277518722` to outward-certified branch upper
`-4283.2348443477`; its overlap-safe cell union term is at most
`-4229.5073186902`.  The verifier checked 23 exact BSP nodes, 12 leaves,
177 leaf-vertex inequalities, and the complete 37,750-cell source geometry.
Artifact:
`out/g4_worst_cell_bsp_depth5_outward.json`, SHA-256
`9e912d87d2614cf376376ed235406bf6a8468e6a4426f11eefbad314f53db759`.
This proves the refinement for one source cell only.

A 128-cell contribution queue plus two rounds of targeted atlas enrichment
reduced the updated diagnostic union from `+67252.4553` to `+32688.4715`;
the largest remaining unprocessed source-cell term is `+32687.2217`.
Most processed cells now close.  Exact BSP replay isolated ten round-three
tuning profiles.  Deep self-tuning closed three, and neighboring-witness
multi-starts closed five more.  Two adjacent profiles remain open under the
current 65-state inner relaxation after 30 broad multi-start runs:

- `(442062,18425,38364,25436,1)`: best combined bound `+451.721903445`,
  hence `563.136968461` bits above the per-profile target;
- `(442348,19721,37834,24384,1)`: best combined bound `+968.431442062`,
  hence `1079.846507078` bits above target.

The broad-reseed artifact is
`out/g4_bsp_round3_last2_broad_reseed.json`, SHA-256
`c92c726e6a555f2f5834163aeb267a2824fb642fef3e3367a64535401955fcda`.
This is evidence about the old upper-bound machinery, not a low-distance
witness.

Subsequent work closed both profiles with a new outer lemma.  The 128 graph
holes occupy distinct band-zero physical groups.  Conditional on a graph bit
`b`, an affected g=4 group has the exact packet polynomial `P^15 P_b`, rather
than paying the adversarial one-bit replacement ratio.  Averaging the 128
hole factors with the exact 24-dimensional graph spectrum gives reusable
multivariate linear-BL and total-spectrum affine witnesses.  Exact-rational
polynomial evaluation plus outward logarithms reports:

- first profile: combined upper `-1052.0540829448`, margin at least
  `940.6390179284` bits;
- second profile: combined upper `-694.1334409284`, margin at least
  `582.7183759120` bits.

Artifacts are `out/g4_hard_profile1_exact_graph_linear_outward.json` and
`out/g4_hard_profile2_exact_graph_spectrum_outward.json`.  The first 128-cell
BSP rerun with 2,114 upgraded graph-aware witnesses closes 125 cells and
leaves three cells with pointwise atlas failures.  Their exact replay yields
nine new integer tuning anchors; graph upgrade closes two immediately, while
seven still need neighboring-witness multi-start tuning (current diagnostic
gaps are 1.4k--11.4k bits).  The next priority is that bounded seven-profile
multi-start round, followed by another 128-cell sweep.  Only after the queue
stabilizes should the BSP batch be integrated into the global outward ledger.

## Theorem-Facing Or Checkable Claims

### Accumulator Warmup

Status: theorem-facing warmup.

The accumulator exact enumerator and contraction are in `innerAcc.tex`.  The
combined warmup theorem says a log-memory random sliding dense outer plus the
accumulator inner gives linear distance for sufficiently small constants.  The
role is explanatory; it is not the intended concrete construction.

### Random Sliding Dense Construction

Status: theorem-facing asymptotic line.

The random sliding dense parent section is `randomDenseConstruction.tex`; it
inputs `outerDense.tex`, `innerDenseScalar.tex`, and `integrationDense.tex`.
The key claims are:

- `thm:dense-dense-integration`: linear distance under explicit parameter
  inequalities.
- `thm:dense-dense-explicit-constant`: concrete `0.109` asymptotic checkpoint.

Current verifier:

```powershell
python scripts\verify_dense_claims.py --delta 0.109
```

Latest checked status: sampled-grid verification passes with worst gap
`-0.003107731647`.

### Structured Local-Code Construction

Status: two theorem-facing finite checkable certificate rows.

The structured parent section is `localCodeStructured.tex`; it inputs the
local-code outer interface, full-split inner interface, two exact certificate
rows, and the BCH projection note.  The frozen full-dimension instantiation is:

- outer: direct sum of `4096` copies of `RM(4,9) [512,256,32]`;
- inner: full-split EBCH `[128,64,22]` with `b=64`;
- length: `N=2^21`;
- target: `delta=.09`, `d=188743`;
- theorem-safe rational/outward first moment: `E[Z_d] <= 2^-37.27`;
- conclusion: some realization is a binary
  `[2^21,2^20,d_min >= 188744]` code.

The promoted field-symbol instantiation is:

- outer: `16383` independent EBCH `[128,64,22]` message blocks plus one
  componentwise-XOR parity block;
- dimension: `2^20-64=1048512`;
- outer minimum weight: at least `44`;
- theorem-safe rational/outward first moment:
  `417/2^50 <= 2^-41.29`;
- conclusion: some binary `[2^21,1048512,d_min >= 188744]` code exists, and
  scalar extension gives the same parameters over `GF(2^128)`.

Current verifier:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
```

Latest checked status: verifier passes and reports
`current_checked_ledger_total_log2,-37.278528`; the underlying stored value is
`-37.278527626...`.  This remains a sharper checked-log diagnostic.  The
theorem-facing rational/outward certificate instead reports
`fullsplit_complete_rational_37_27_bits_status,PASS` and uses `2^-37.27`.

### Smaller Exact Outer And BCH Projection Rows

Status: exact EBCH-128 parity-block alternate; BCH256/BCH512 projections
remain heuristic.

The primary smaller-block alternate uses `16383` independent EBCH
`[128,64,22]` message blocks and one componentwise-XOR parity block.  It is an
explicit codimension-64 subcode of the ambient direct sum of `16384` EBCH
blocks, has outer weight at least 44, and has dimension `2^20-64=1048512`.
The standalone rational/outward verifier
`scripts/certify_ebch128_outer_fullsplit.py` covers every outer weight and
reports the complete diagnostic exponent `-41.306810170135`; exact integer
arithmetic proves `417/2^50 <= 2^-41.29`.  Scalar extension through the same
binary generator gives the identical dimension and distance lower bound for
vectors of `GF(2^128)` symbols.  Encoding the parity requires 1048512
128-bit XORs and no field multiplication.  The ambient full-dimension lane
still has only 20.80 bits; a checked weight-22 suffix family lower-bounds its
first moment by `2^-36.493931`, so it cannot reach 40 bits by upper-bound
cleanup.  A random codimension-20 subcode remains a 40.80-bit theoretical
alternative.  This parity-block lane is now promoted as
`thm:fullsplit-ebch-xor-parity-certificate-009`; see
`scripts/ebch128_outer_fullsplit_certificate.md`.

The BCH256 and BCH512 rows are spectrum-model projections, not theorem claims.
They should not be promoted until an exact spectrum or rigorous low-weight
envelope replaces the modeled spectrum and the affected tails are recomputed.
The BCH512 code construction is no longer heuristic:
`scripts/build_bch512_256_candidate.py` exactly regenerates and authenticates a
complement-symmetric `[512,256,>=62]` subcode.  What remains heuristic is only
its spectrum.  The corrected projection uses even weights and the
`2^(k-(n-1))*binom(n,w)` baseline.

The first exact obstruction audit is recorded by
`scripts/audit_bch512_spectrum_obstruction.py`.  The extended parent has dual
distance at least `16`, but the best strength-15 Christoffel/constant-weight
generic bounds at weights `62,94,118` remain `56.84,59.92,81.89` bits above
the current projection targets.  Therefore a theorem upgrade requires
BCH-specific cancellation, specialized spectrum enumeration, or a different
local code with a certified spectrum; another parameter-only Delsarte estimate
is not expected to close the recorded gap.

## High-Risk Items Already Handled

- The broad section structure is now real in the compiled paper, not only in
  planning notes.
- The structured local-code line is now one parent section with local-code
  independent outer and inner interfaces before the RM/EBCH instantiation.
- The random sliding dense line is now one parent section containing the outer,
  inner, and asymptotic integration.
- The dense `0.109` tiny-window proof now uses an explicit
  `xi_tiny=1/4` handoff instead of hiding the linear-window `xi=8`.
- The dense `0.109` low-weight window now uses the polynomial-free global
  outer envelope `A_h <= 3^h`, avoiding the loose `n^{O(1)}` prefactor.
- The scalar dense inner section now includes a dependency map identifying the
  ON/OFF lemmas as local machinery and pointing to the packaged envelopes
  consumed by the integration theorem.
- The dense integration handoff now uses an explicit checked endpoint
  `eta_1=0.99` in the linear-window criterion and treats the top endpoint
  `h=n` separately instead of applying the outer linear exponent at `eta=1`.
- The full-split finite certificate now states the block-time convention
  (`T=B-i+1` from the first active block), the `late_blocks=5949` split, and
  the unconstrained terminal-state convention in both paper prose and manifest
  metadata.
- The generic first-moment theorem and the block-outer interface now bound the
  joint event that the inner restriction is noninjective or the image has
  minimum distance at most the target.  This repairs the former implicit
  assumption that every square, length-preserving inner map is invertible.
- The full-split one-step turnoff law now states its probability space:
  conditional on block occupancies and preceding transition weights, the
  input and state supports are independent uniform subsets of their specified
  sizes.
- The finite theorem is factored through named complete `h<=500` and `h>=501`
  certificate lemmas; the theorem itself only combines their exact dyadic
  thresholds and concludes the certified dimension and distance.
- The finite certificate states its audit layers explicitly.  The legacy
  manifest total still re-sums floating logarithmic rows within stored
  tolerances, while independent exact/rational and outward-rounded artifacts
  now cover both the complete `h<=500` family and every `h>=501` post-prefix
  family.
- The finite manifest now gives explicit `PASS` status fields for the exact
  RM support check, the prefix interior audit, and the prefix placement-ratio
  certificate, rather than leaving those checks implicit in successful
  execution.
- The finite certificate prose now distinguishes the default manifest verifier
  from full regeneration: the default pass re-sums committed row artifacts and
  manifest interval constants, while the high/complement/early/late interval
  families have opt-in recomputation flags.
- A full opt-in interval recomputation audit is now tracked at
  `scripts/fullsplit_interval_recompute_audit.md`.  It recomputes all 45
  high, complement-high, early-postprefix, early-accelerated, and late
  postprefix rows and records `PASS` under upper-bound semantics; no row
  recomputed above its stored manifest bound, and several high-interval rows
  recomputed to slightly safer values.
- The dominant RM/EBCH prefix mechanism is now isolated in the manifest:
  after exact RM reweighting the \(h=32,r=1,\mathrm{gap}=1..4000\) row is the
  peak, the full \(32\le h\le500,e\le8\) prefix family is only about
  `0.002423` bits above that peak, and the companion ridge-shape check records
  the first-ridge ratio thresholds.
- The same peak-to-family inflation is now a thresholded verifier check:
  total family inflation, total non-peak remainder, \(h=32\)-slice inflation,
  \(h=32\)-slice remainder, and \(h>32\) separation all have explicit gates.
- The dominant peak row arithmetic is now recomputed as its own manifest check:
  exact \(A_{32}=4096\,A^{RM}_{32}\), exact first-gap placement
  numerator/denominator, and the checked `T=5949,H=31` inner knot table value.
- The dominant peak now also has an independent rational upper-bound audit.
  It constructs the EBCH split law with exact fractions, proves the finite
  termination atom is below the rational cap `2^-63` at the monotone endpoint,
  uses the rational Chernoff pole `2333/2373`, rounds the survival probability
  upward to an explicit 96-bit dyadic, bounds the `e=1..8` termination terms by
  `sum_e binom(H+1,e)2^(-63e)`, and checks the complete RM/placement/inner peak against
  `49/2^43` by exact cross multiplication.  The resulting displayed upper
  exponent is about `-37.385767`.  This remains a useful independent peak
  cross-check; the complete low-weight family is now covered by the aggregate
  rational certificate below.
- The dominant `T=5949,H=31` inner knot is now regenerated inside the verifier
  from the EBCH spectrum and the `e<=1` full-split episode formula under the
  stored effective-turnoff convention, then compared against the knot table.
- The ultra-late prefix `T<5949` is now a thresholded verifier check using
  exact RM direct-sum outer coefficients; the total, peak, split, above-split
  mass, and above-split gap all have explicit gates.
- The same ultra-late prefix now has a fully rational aggregate check: the
  verifier sums all 114 supported terms
  `A_h binom(64*5949,h)/binom(2^21,h)` as exact fractions and proves the total
  is below `2^-41` by integer cross multiplication.  Its displayed exact-sum
  exponent is about `-41.113442`; logarithms are used only for reporting.
- The complete `h<=500` first-active family now has an independent
  exact/outward-rounded certificate.  It covers all six gap buckets, keeps
  `e=0..8` in the three prefix buckets and `e=0..16` in the three early
  buckets, and covers the remaining `e>=9`/`e>=17` tails.  Exact RM counts,
  exact placement sums, a 1024-bit outward dyadic inner envelope, and exact
  final cross multiplication give a rational upper bound whose diagnostic
  logarithm is about `-37.276548513006`.  The theorem-safe threshold is the
  exact rational `6793/2^50`; the integer inequality
  `6793^100 <= 2^1273` proves it is at most `2^-37.27` without evaluating a
  transcendental logarithm.
- The complete `h>=501` post-prefix aggregate now has an independent
  rational/outward-rounded certificate.  Exact rational arithmetic handles the
  EBCH split law, MGFs, turnoff envelope, adjacent-ratio peak locations, and 56
  endpoint-in-`T` monotonicity checks.  Decimal logarithms are evaluated with
  directed outward rounding at precision 100.  The diagnostic aggregate upper
  logarithm is `-183.800029543464`, and the theorem-facing gate is `2^-180`.
  The default verifier validates the SHA-256-protected artifact; full
  regeneration is available with
  `python scripts\certify_fullsplit_postprefix_rational.py --recompute-artifact`.
- Both rational certificate paths structurally validate the committed RM and
  EBCH weight tables and record their SHA-256 fingerprints.  The post-prefix
  artifact is rejected if either source table changes.  These checks bind and
  internally validate the inputs; they do not independently derive either
  spectrum from a code definition.  In particular, the EBCH table's
  mathematical provenance remains a declared input.
- The two independent bounds are now combined exactly:
  `(6793*2^130+1)/2^180`.  The default verifier checks
  `(6793*2^130+1)^100 <= 2^14273`, proving that the complete `1<=h<=N`
  rational/outward first moment is at most `2^-37.27` with integer arithmetic.
- The RM/EBCH finite result is stated with the combined rational/outward first
  moment as its theorem-facing bound; the sharper checked-log ledger remains a
  diagnostic, and BCH rows remain separated as projections.
- The current LaTeX log is clean after two passes: no warnings, overfulls,
  undefined references, multiply-defined labels, or rerun requests were found.

## Explicit TODOs Before Stronger Claims

### Riffle eight-bit packet co-design branch

Status: active exploratory proof branch, not a theorem claim.

The corrected construction keeps both binary local codes at the committed
extended BCH `[128,64,22]` parameters.  Only the global permutation atom is
eight bits: eight permuted packets are concatenated into each unchanged
64-bit recursive-inner input.  The abandoned `[16,8,4]` small-inner probe is
not part of the proof target.

Current exact endpoints are recorded in
`explorations/riffle_group_chain_proof.md`:

- `scripts/certify_packet8_permutation.py` proves the former universal
  late-suffix event is at most about `2^-72.965` and checks the exact
  capacity-eight occupancy law;
- `scripts/analyze_packet8_impulse.py` exhaustively checks every nonzero
  packet value and slot, including an inner total impulse of at least 22;
- `scripts/certify_packet8_dense_orbit.py` gives an exact-integer balanced
  packet-profile endpoint with about 130330 bits of margin beyond the
  `2^-40` target, without crediting the graph codimension.

The low-entropy repeated-packet endpoint currently uses an exact local DP and
floating pole/matrix arithmetic.  At conservative nonzero density `21/128`,
the repeated `0xff` family has a diagnostic margin of roughly 65329 bits after
its 64-dimensional family count and the 24-bit graph charge.  The remaining
global gate is a packet-profile overlap ledger: high packet-value entropy must
enter the inverse-orbit/Hamming-ball branch, while every lower-entropy outer
profile must enter a hardened packet transfer/outer-spectrum branch.  No
intra-packet mixer is currently justified by the proof diagnostics.

The repeated-value scan is now exhaustive over all 255 nonzero bytes.  A
cross-checked byte trellis finds `0x18` worst at approximately
`-42116.902488` after the provisional `64-24` family charge; adjacent two-bit
bytes occupy the next six rows.  The generalized bounded-alphabet OFF/LIVE
probe retains exact multinomial conditioning and exact turnoff coefficients.
An equal `0x18/0x30` profile remains near `-42046.72` after a provisional
`128-24` charge, and a related three-value linear alphabet is safer.  These
are diagnostic rows, not outer-family counts.

The packet-weight-profile transfer and the exact high-profile certificate are
also implemented.  The latter checks both inverse-orbit/Hamming-ball and a
full-message bijection moment with integer comparisons.  At the balanced
profile their union exponents are respectively about `-130370.69` and
`-133156.25`.  On the inspected repeated-to-balanced interpolation, the
bijection branch becomes safe near 85.6 percent and the orbit branch near
85.8 percent.  Before adding an outer profile count, the 85-percent profile
was uncovered by 5919--8647 bits because both high branches unioned all `2^K`
messages.  This identified the need for a three-band packet-profile enumerator
using the exact BCH band spectra, within-band coordinate permutations, and
lane permutations; it was a proof-method gap rather than distance evidence.

That outer enumerator now has a global Finner form:
`A_outer(t) <= 2^K S_3(t)^(B/3)`, with the exact random-lane group polynomial
inside `S_3`.  Its binary64 coefficient probe closes the repeated-to-balanced
interpolation from roughly 75 percent onward (`-8426.95` bits at 75 percent
and `-52129.64` at 80 percent), so the intermediate and exact high-profile
branches overlap.  Graph/puncture correction and outward arithmetic remain to
be inserted after full numerical coverage.

The surviving endpoint phase has also been isolated.  Exact two-band binary
ranks for all tested sizes through 64 tiles have nullity one, consisting only
of the zero/all-ones packet phases.  The exact `[85,64]` 42/43 split spectrum
has been constructed from its 21-dimensional dual.  The earlier structured
Holder estimate was corrected: every constant-packet variable has degree
eight, and the resulting degree-8 Finner witness lowers the universal
first-two-band expected solution-count bound from about `2^152832` to about
`2^4096.001`.  A cumulative support fugacity plus the 65-state inner transfer
is diagnostically safe from 208 active packets onward (`-43.93` bits at the
handoff and rapidly improving).  The exact all-ones trajectory itself remains
safe by about 681253 bits.

The former `21..207` exact-constant residue has been removed.  Exact vertex
counting, pair capacity one, and fixed-band outside distances prove that every
nonzero final constant-packet word has at least 40 active packets, including
the 128 punctures and graph replacements.  The exact averaged 41/43 punctured
spectrum has digest
`1e5e7bf5679098703b0d90c41e7718be1f1ceb4c9224ae7dcd7f976ea55d51e9`.
A heterogeneous degree-8 witness uses exactly 16256 ordinary and 128
punctured data factors; for a fixed final packet assignment, graph matching
costs at most `2^-24`.  The OFF-start renewal transfer is already about
`-65.61` bits at support 40 and improves rapidly.

The all-ones complement is closed by an exact integer support-graph potential.
Every local transition obeys
`emitted >= 6 + phi(next)-phi(current)-3*zero_packets`, so `q` zero-packet
defects emit at least `196608-3q`.  This is deterministically above the target
through `q=2621` (188745 at the endpoint), and the support-dependent
probabilistic row at `q=2622` is safe by thousands of bits.  Exact auxiliary
spectra retained for audit are the dimension-21 band-two-zero split (digest
`d9f9a66f0bf8f16a543b8cdeff8ee33f4febeaf6ba733591cdb50dc80379252b`)
and the `[86,64]` band-(1,2) split (digest
`859217d417f433e1361aba51a40e05110b55d25c2a2a6473015754a76ff79b17`).

The former intermediate global-correlation gap is now closed numerically by a
stronger algebraic inequality.  The three fixed-band projection kernels in
the 64-bit EBCH message space have dimensions `22,21,21`, and their joined 64
basis vectors have exact rank 64 (certificate digest
`f2c0a6e227151fac4ecdadd640c81272c2fd6888ffedbf9c3adc7317d3acc1be`).
Thus they form a direct-sum decomposition.  For every global message subspace
`V`, this gives

```
sum_(band,tile) dim L_(band,tile)(V) >= 2 dim(V).
```

The finite-field linear Brascamp--Lieb inequality therefore applies to the
tile factors with coefficient `1/2`, improving the previous exponent-three
Finner outer bound to

```
S_2(t) = 2^-64 sum_w C(64,w) R_w(t)^2,
A_outer(t) <= 2^K S_2(t)^(B/2).
```

`scripts/certify_ebch_three_band_kernel_decomposition.py` verifies the exact
linear algebra.  The binary64 coefficient probe in
`scripts/probe_packet8_three_band_linear_bl2.py` lowers the rounded 70-percent
profile from `+35380.312595` to `-41409.964997` bits after the exact bijection
moment at the symmetric coefficient point.

The direct-sum proof gives the larger linear BL polytope as well.  If the band
coefficients are `p_0,p_1,p_2`, the dimension condition holds whenever every
pair sums to at least one.  Indeed, writing `K_b=ker(P_b)` and
`r_b=dim(V intersect K_b)`, directness gives `sum_b r_b <= dim(V)`, so

```
sum_b p_b dim P_b(V)
 = (sum_b p_b) dim(V) - sum_b p_b r_b
 >= dim(V)
```

whenever `sum_b p_b - 1 >= max_b p_b`.  Optimizing the valid family
`(p_0,p_1,p_2)=(1-p,p,p)`, `1/2 <= p < 1`, is especially effective near the
constant-packet endpoint.  Direct fugacity tuning of the robust 65-state
inner transfer in `scripts/probe_packet8_weight_profile_state_tuned.py` removes
the remaining optimizer mismatch.  On the rounded repeated-to-balanced grid,
fixed witnesses now give

```
balanced fraction     outer linear BL      inner transfer      combined
0.10                    28591.888265       -146606.267293     -118014.379029
0.20                   221262.820468       -231578.802177      -10315.981710
0.30                   352367.269866       -354645.760973       -2278.491106
0.40                   463414.396837       -471045.463052       -7631.066215
0.50                   571886.406714       -634336.009369      -62449.602655
0.60                   678351.863905       -692302.522810      -13950.658905
```

The 5-percent row is safer still: its optimized outer coefficient is
`-90126.453793` before the inner charge.  Thus the entire inspected
repeated-to-balanced interpolation grid overlaps numerically with the exact
constant endpoint and the high-profile branches.  This is not yet the full
nine-dimensional profile ledger: a rigorous reduction from every profile to
a finite outward-certified cover is the active gate, followed by the exact
graph/puncture replacement.  Unequal *degree-three Finner* weights improve the
old 70-percent row by only 0.35 bits, and exact one-block tile-pair
correlations are below binary64 resolution; the successful mechanism is the
linear BL kernel polytope plus direct inner tuning.  Evidence continues to
favor proof success, not construction failure or intra-packet mixing.

The first full-simplex diagnostics are now explicit.  All eight nonzero pure
packet-weight vertices admit positive-floor witnesses with hundreds of
thousands of bits of margin; at inactive fugacity floor `0.1`, the weakest
pure row is still about `-152432.73` bits.  Regularized midpoint and adaptive
witnesses reduce the sampled landscape but show that a flat atlas is
inefficient: low-dimensional faces require support-aware fugacity floors.
The coupled universal choice `C(8,j)t_j f_j=1`, which would cover the whole
simplex with one witness, fails by roughly 960000 bits and is retained only as
a discarded minimax diagnostic.

One zero-heavy profile exposed the importance of retaining both rigorous
outer mechanisms:

```
(150014,48151,48529,15450,0,0,0,0,0).
```

The fully tuned robust inner exponent is about `-334107.67`.  Linear BL alone
gives outer exponent `406756.39` and leaves a `+72648.72` gap, but the exact
16384-block EBCH total-spectrum envelope gives `251821.02`, closing the same
profile by about `-82286.66` bits.  The simplex probe and adaptive atlas now
evaluate both outer bounds.  This was an integration defect in the diagnostic
cover, not evidence against the construction.  The active numerical gate is
therefore a hierarchical support-stratified atlas (sharp face witnesses plus
positive-floor neighborhood witnesses), followed by a genuine convex-cell
certificate; sampled coverage alone is not a proof.

The first paired cutting-plane run confirms that hierarchy quantitatively.  At
profile `(205799,0,0,32424,0,0,0,23921,0)`, the broad inactive-floor `0.2`
witness misses by `+133722.39` bits, while a separately tuned sharp
support-local witness closes by `-15288.89` bits.  Five subsequent sharp/broad
pairs all close their selected rows; their sharp margins range from about
`-21586` to `-251096` bits.  The reusable atlas cache avoids rebuilding every
outer/inner fixed witness on each resume.  These remain sampled binary64
diagnostics, but they identify support stratification as necessary proof
organization rather than optional optimization.

The profile simplex now has an exact finite ordered-chamber decomposition.
For each ordering of the nine class counts, the corresponding chamber is a
simplex whose vertices are the uniform profiles on the nine nested nonempty
subsets.  Across all `9!` chambers there are only `511` distinct uniform-subset
vertices.  Chambers beginning at the forbidden all-zero outer word are
clipped at outer weight 21, adding only `255` distinct clipped vertices.
`scripts/probe_packet8_profile_ordered_chambers.py` exhausts this geometry
without sampling.  The paired broad/sharp atlas now covers all `511/511`
uniform vertices and all `255/255` clipped vertices in binary64.  The final
uniform-vertex worst value was `-252.198771` versus the per-profile target
`-168.700990`.  Vertex closure is diagnostic because the fixed witnesses and
comparisons are not yet outward rounded, but it is exhaustive over this finite
vertex set.

Whole-edge coverage is also stronger than the former shared-endpoint test
showed.  `scripts/probe_packet8_profile_edge_interval_cover.py` uses convexity
of each fixed-witness exponent along an edge, computes its safe interval, and
unions the intervals of all witnesses.  `scripts/probe_packet8_edge_interval_atlas.py`
closed all `2287/2287` relevant Hasse edges in binary64.  An ordered chamber is
an eight-simplex, however, so its one-skeleton also contains every edge between
strictly nested prefix subsets, not only Hasse neighbors.  The complete
comparable-edge family has 18405 members.  The last completed full audit plus
two 20-anchor long-chord passes cover `17077/18405`; 1328 remain.  The monotone
residual list is stored in `out/packet8_comparable_edge_residual.json` and
validates its edge-cover source digest and witness-name prefix before reuse.
Every targeted chord midpoint closed; this is substantial evidence against a
one-dimensional construction obstruction, but full one-skeleton coverage is
not complete.

The first higher-dimensional audit is
`scripts/probe_packet8_nested_face_landscape.py`.  It exhausts all 198580
strict nested-prefix triangle faces.  Before face-specific tuning, 17526
triangles had one fixed witness safe at all three vertices, and
`193121/198580` barycenters were covered.  Five greedy face anchors raised the
barycenter count to `193258/198580`, but the fifth selected face exposed the
first profile that did not close under the current generic transfer witness:

```
face     100 < 1e0 < 1ee
profile  (0,12483,12483,12483,0,34328,34328,34329,121710)
target   -168.700990
broad    +70977.391179
sharp    +23165.930932
```

The sharp support-local row improves the old atlas by roughly 71000 bits but
still misses the per-profile target by about 23335 bits.  This is a binary64,
coarsely optimized upper-bound failure, not a bad codeword or lower bound on
ensemble failure.  It is nevertheless the first material warning that the
generic packet-weight-profile atlas may be an insufficient proof method.
Atlas growth is paused.  Before resuming the full proof, the planned
construction-health diagnostic is a small-length expected input-output weight
enumerator adapted from `enumerator_paper/turbo.tex` and
`enumerator_paper/AccPoly.tex`, using exact packet-orbit denominators rather
than bit-weight binomials.  No C++ construction change or intra-packet mixer
has been made or justified.

An independent diagnostic review changed the first experiment, without
changing the construction.  An exact packet IOWE must track the counts of all
256 concrete byte values; the nine packet-weight counts are not generally an
exact interleaver orbit.  Moreover, the 256-tile graph-hole rule fixes the
full layout size, so smaller end-to-end instances are surrogates rather than
literal scaled copies.  The exact composition remains

```
E[E_(w,h)] = sum_A O_(w,A) I_(A,h) / (M! / product_v A_v!),
```

but this is not the cheapest first construction-health test.

`scripts/probe_packet8_hard_profile_product_kernel.py` now attacks the first
non-closing face directly in a sampled full-size legal outer instance.  It
uses the unique cyclic-basis message `0x8e63cf44efd4fa21` whose EBCH word is
all ones, includes the actual random 24-bit graph coupling and 128 graph-hole
bits, and anchors the target's 121710 weight-eight packets only at non-hole
positions.  Those homogeneous constraints split across the 16384 data EBCH
blocks.  Exact elimination gave product-kernel dimensions 88097 and 87992 on
two independent layout seeds; the median local nullity was five.  Thus the
hard face is not dismissed by a shortage of local EBCH degrees of freedom.

For seed 20260810, one exact graph-coupled basis-coordinate search started at
profile `(0,0,0,0,0,0,0,84,262060)` and reached
`(0,14,274,2076,2577,27649,40893,43808,144853)`, reducing packet-profile L1
error to 83528 before single-direction greedy descent stalled.  This is not
an attained hard profile and not a low final codeword.  It does establish a
concrete exact search path and shows that the product kernel moves materially
toward the face.  The next attack step is local-domain or paired-direction
descent inside each small EBCH kernel, followed by an exact inner-weight test
only if the target profile is attained.  Profile-conditioned thermodynamic
integration remains the best fallback for measuring how much of the roughly
0.712 bits-per-inner-group deficit belongs to the robust inner relaxation.

That fallback has now been prototyped.  Exact whole-local-kernel enumeration
through nullity eight lowered the seed-20260810 profile L1 error to 73966;
100000 random improving direction-pair trials lowered it only to 73840.  A
second full layout/anchor seed stopped at 75106, with the same qualitative
deficit in packet weights one through three.  Nullity-ten enumeration reached
the same basin.  The random-anchor attack therefore did not attain the face,
but its failure is not an infeasibility result because an attainable face may
have a highly correlated weight-eight packet set.

`scripts/probe_packet8_hard_profile_thermodynamic.py` samples the exact frozen
systematic inner conditional on the scaled hard nine-class profile.  Its
ensemble includes packet order, uniform concrete byte values within each
weight class, and every 64-coordinate state permutation.  Parallel tempering
and thermodynamic integration gave the following mixing-sensitive Chernoff
diagnostics at beta 0.6:

```
inner groups   estimated Chernoff log2 exponent per group
16             -19.5102
32             -20.1885
32, seed 2     -20.1518
64             -20.7665
```

The longer-burn-in second 32-group run gave -19.9891 over all measurement
batches and -19.9330 using only the last half.  It recorded ten complete
hot-to-cold-to-hot replica trips; the last four cold-end energy batches were
stable near 861.3.  The outer hard-profile envelope costs 16.4328 bits per
group and the per-profile target requires the inner to supply only -16.4380
bits per group.  Thus the least favorable stabilized finite-size estimate has
about 3.50 bits per group of diagnostic room.  This strongly suggests that
the +23335-bit face miss is slack in the robust inner transfer, not evidence
of poor true distance.  It is not a proof: the MCMC means, quadrature, and
finite-size extrapolation are not certified, and the exact outer codeword
profile was not attained.

The exploration/proof follow-up has now isolated and certified the missing
inner correlation.  For incoming state support `R`, current input `U`, and
accumulator drive `W=U xor R`, `scripts/probe_packet8_drive_stratified_caps.py`
computes the exact local mass `H[q,d,y]` and the maximum point mass `c[q,d]`.
For a fixed emitted weight and next-state weight, all drive strata must share
the same systematic-EBCH split column.  After the substitution
`n[d,r]=x[d,r]/c[q,d]`, the resulting upper problem is a rank-one
transportation problem and is solved by the rearrangement inequality.  The
lemma and its Collatz composition are written in
`explorations/packet8_drive_stratified_transfer.md`.

At the first hard face, the binary64 fixed-witness diagnostic improves the
inner MGF by `81113.802017` bits (`2.475397` bits per group), moving the old
combined value from `+23165.930932` to `-57947.871085`.  This diagnostic is
reproduced by `scripts/probe_packet8_hard_face_drive_witness.py`.

More importantly, `scripts/certify_packet8_hard_face_drive_inner.py` evaluates
the local dynamic programs with one-ulp directed rounding, solves all final
transportation problems with exact rational arithmetic, and checks a frozen
positive Collatz vector with outward Decimal logarithms.  The certified inner
MGF upper is `300323.502865628` bits.  The end-to-end script
`scripts/certify_packet8_hard_face_drive_composition.py` uses explicit integer
outer packet variables, the symmetric `p=1/2` three-band linear
Brascamp--Lieb inequality, and a worst-adjacent-variable charge for each of
the 128 graph-hole bit replacements.  Its complete hard-face union ledger is
at most `-48839.7238144`, leaving `48799.7238144` bits below the `2^-40`
budget.  Thus the first non-closing face is now outward certified without a
mixer or a construction change.

The remaining-face audit is informative but not a complete atlas theorem.
The single hard-face witness covers only 49 of the 5322 previously uncovered
nested-triangle barycenters because its inactive-class fugacities are
support-specific.  Retuning the same shared-drive lemma at the two worst,
structurally different residuals gives binary64 combined values
`-67516.355816` and `-164678.161830`.  This shows that witness reuse, rather
than the local abstraction, is the immediate obstruction.  Closing the whole
simplex still requires generating a retuned shared-drive atlas and outward
certifying its finite face cover; the current hard-face certificate must not
be described as that global cover.

## Permutation-group curve (`g=1,2,4,8,16,32,64`)

The frozen no-mixer inner proof machinery is now parameterized by the
permutation atom width `g`, without changing the systematic `[128,64,22]`
local code, accumulator, graph-hole rule, or construction semantics.  The
new exact-local factorization fixes the drive atom `W=U xor R`.  Accumulator
emission depends on the concrete `W`, while the weighted number of compatible
`(U,R)` pairs depends only on `wt(W)`.  Consequently
`scripts/packet_group_drive_stratified.py` replaces direct `4^g` atom
enumeration by a polynomial parity/weight DP and supports every divisor

```
g in {1,2,4,8,16,32,64}.
```

`scripts/packet_group_outer_profile.py` and
`scripts/packet_group_profile_bound.py` likewise parameterize the exact
profile normalization, symmetric outer moments, linear and total-spectrum
outer bounds, graph-replacement charge, point caps, shared EBCH-column
transport, and Collatz witness.  At `g=8`, the generic histogram, point caps,
shared operator, and composed hard-face value agree with the original
packet-8 implementation to binary64 roundoff.  The outward generic
certificate reproduces the original conclusion.

The following is the independently retuned **mapped-hard-family curve**.  It
is useful for comparing `g`, but it is not a worst-profile curve and must not
be presented as one.

| `g` | binary64 diagnostic margin (bits) | outward-certified margin (bits) |
|---:|---:|---:|
| 1 | 255021.230 | 243698.772 |
| 2 | 312455.002 | 296562.087 |
| 4 | 355032.237 | 342997.457 |
| 8 | 325861.137 | 322812.282 |
| 16 | 235852.182 | 233590.167 |
| 32 | 253968.324 | 247807.296 |
| 64 | 98959.943 | 95923.004 |

All seven rows pass outward for that one mapped family.  The nonmonotone
shape is genuine optimization behavior: changing `g` changes both the number
of atoms and the exact drive/state compatibility table.  It is not evidence
that `g=4` has globally better distance than `g=8`.

The endpoint and residual diagnostics sharpen the status at `g=64`:

- all 64 nonzero pure packet-class endpoints close in binary64; the weakest
  is class 32 with `127701.688` bits of margin;
- a first mixed residual closes with `9401.378` bits and a low-density
  residual closes with `64939.779` bits;
- the first apparent nonpositive mixed profile initially missed by
  `107317.176` bits, but this was fugacity-tuning slack.  Pinning its 48 absent
  classes to the zero boundary and coordinate-refining only its 16 occupied
  non-anchor classes changes the verified inner exponent from
  `-461344.222` to `-601585.824`.  With the exact-total-weight outer value
  `567957.301`, the combined value is `-33628.522` against target
  `-704.096`, leaving `32924.426` bits.

`refine_inner_sparse_fugacities` implements that reusable sparse-profile
step.  It is a sound change of Cauchy parameters, not a construction change
or a new probabilistic assumption.

### Precise global-cover obstruction

The generalized work has **not** produced the requested finite global
profile cover, so there is not yet a worst-profile security-margin curve for
all seven `g`.  This is already unresolved at `g=8`: the outward result above
certifies the first hard face, while the existing nested-face audit still has
5273 uncovered barycenters before retuning.

The new experiments isolate why a naive polynomial full-support cover does
not solve this.  A fixed witness has profile bound

```
F(a) = C - <a,q> - log2 Q(a).
```

It is convex in `a`.  Therefore a fixed (or convexly mixed) witness can be
checked on the vertices of a total-weight slab.  The script
`scripts/probe_packet_group_global_mixture_cover.py` implements this finite
LP reduction, including the exact weight-21 clipping and locally frozen
exact-spectrum outer poles.

However, a middle class exposes an unavoidable tradeoff in the present
point-cap abstraction.  Aggregate class balancing wants
`f_j proportional to 1/binom(g,j)`.  At `g=64`, this makes the representative
fugacity of class 32 about `2^-60`, while rare classes have much larger
representative fugacities.  The exact per-drive point mass then sees those
rare representatives even when their aggregate fugacity is small.  The
opposite normalization (`f_j=1` for the favored middle class) keeps the point
cap sharp but gives the favored class aggregate mass `binom(64,32)`, making a
full-support witness unusable away from that face.  Decreasing a class-mass
floor to `10^-3` still left representative middle-class endpoint probes
roughly 1.35 million bits high, whereas the corresponding zero-support sharp
endpoint closes by 127702 bits.  This is a support-decomposition obstruction,
not evidence of a low-weight codeword.

The next proof step is consequently not more scalar pole tuning.  It is a
support-stratified finite cover (with zero fugacities outside each face), or a
stronger local inequality that retains enough concrete-drive mass
distribution to combine support faces without one global point cap.  Until
that is done and outward checked, only the mapped-family rows and selected
faces above are certified; no global distance theorem for arbitrary `g`
should be claimed.

### `g=2` ladder rung: current end-to-end status

The first ladder target is now the complete `g=2` profile simplex.  Here
`M=N/2=1048576`, a profile is `(a0,a1,a2)`, and there are
`binom(M+2,2)` profiles.  A uniform union allocation for 40 total security
bits therefore requires every profile branch to be at most
`-79.00000412758034` bits.

The diagnostic cover combines three kinds of valid fixed branches: frozen
outer packet-profile witnesses, frozen shared-drive outer/inner witnesses,
and the exact all-message bijection moment.  Exact edge scans and all
pointwise interior profiles optimized so far are safe by tens of thousands to
hundreds of thousands of bits.  This is strong evidence for the `g=2` rung,
but it is not yet a global certificate: the current finite cell atlas still
leaves unresolved integer profiles.

One important implementation defect in the diagnostic cover was found and
fixed.  The edge scan often selected the exact-spectrum total-weight outer
Cauchy bound, but `build_packet_group_curve_witnesses.py` froze only the much
weaker linear-profile outer branch.  A profile `(891220,157354,2)` consequently
appeared to miss the target by `6111.150` bits.  Freezing the winning
total-weight pole as the globally valid affine witness

```
C(r) - (a1 + 2 a2) log2(r)
```

changes that same profile to `-118933.545`, or `118854.545` bits of margin.
The reusable implementation is `fixed_total_weight_outer`; the edge builder
and full-profile atlas now preserve the winning outer branch.

With 63 rebuilt edge witnesses and 19 full-support witnesses, the latest
30000-node diagnostic cell pass still does not close.  Its next exact
unresolved profile is `(958312,60796,29468)`.  The preceding residual
`(914042,31395,103139)` optimizes to `-82986.246`, leaving `82907.246` bits,
so the observed failures remain missing-atlas geometry rather than pointwise
bound failures.  The next step is to automate the residual-to-witness loop,
then outward-round only the selected witnesses and cell duals.  No end-to-end
`g=2` distance certificate should be claimed before both steps finish.

That automation is now implemented.  `probe_packet_group_g2_triangular_cover.py`
uses the exact five-vertex integer hull, one fixed witness (or one fixed convex
mixture) per triangle, and emits residual profiles directly consumable by
`probe_packet_group_fixed_atlas_parallel.py`.  Peach provides 16 physical
cores for one coordinated batch-tuning job.  Three 16-profile rounds reduced
sampled uncovered profiles from `839` to `248`, then `46`, then `0`; every one
of the 48 new pointwise witnesses closes, with the weakest observed margin
still `55742.673` bits.  This is strong diagnostic evidence that no pointwise
`g=2` obstruction has been encountered.

The diagnostic integer cover now **does close**.  Witness-anchor Delaunay
seeding, followed by exact mesh auditing and adaptive refinement, produced
37124 convexly covered triangles, of which 21784 use fixed mixtures.  The
remaining 4517 terminal cells were exhaustively lattice-enumerated: they have
13551 cell-profile occurrences but only 6198 distinct profiles, and all 6198
have explicit safe singleton witnesses.  There are no nonterminal unresolved
cells and no uncovered integer residuals.  Both
`complete_interior_integer_cover` and `complete_global_integer_cover` pass.
The continuous real hull is intentionally not claimed covered; the hybrid
triangle-plus-singleton theorem needs only the finite integer profile set.

`certify_packet_group_triangle_ledger.py` checks exact mesh geometry, terminal
lattice enumeration and ownership, source hashes, support eligibility, exact
rational mixture weights, outward component and weighted vertex bounds, and
the final 40-bit union ledger.  The checkpointed parallel outward run on Peach
has now completed and **the `g=2` end-to-end integer-profile certificate
passes**.  The rigorous global union interval is

```
[-40.09421546136248628735, -40.09413878241724607147].
```

Thus at distance `floor(0.09 N)=188743`, the construction has at least
`40.0941387824` bits of certified first-moment security: it exceeds the
requested 40-bit threshold by `0.0941387824` bits.  The worst owned profile is
`(939836,43595,65145)` and is certified by a fixed two-witness mixture.  The
outward verifier hardened 54 distinct witnesses, checked 89918 mixture-vertex
inequalities, 159834 component-witness vertex inequalities, and 6198 terminal
singleton inequalities.  Its exact geometry audit covers all 41641 cells
(37124 convex cells plus 4517 terminal cells) and all
`binom(1048578,2)=549757386753` packet-weight profiles.

Canonical artifacts:

- diagnostic hybrid ledger: `out/g2_triangle_hybrid_complete.json`, SHA-256
  `d45cf8ffacdabf7d74669db60eea958b01a5ee30fc72b290b4cab8ca986e2333`;
- outward certificate: `out/g2_triangle_hybrid_outward_parallel.json`, SHA-256
  `c03d188708ee5dbc65c586cd7d35eb5cb554ed12af47b89f77d4d032f1580875`;
- outward inequality ledger SHA-256
  `0841a46c7b45766cefaefa0142ab733e7f3ad3cd765e7bbed5660827499d0e75`.

The uniform worst-profile allocation above has now been superseded by a
cell-local weighted union ledger.  For each convex triangle, the verifier
multiplies its outward maximum vertex bound by the exact closed-triangle
lattice count; it then adds each globally deduplicated terminal singleton
once.  Shared closed-cell boundaries are deliberately counted more than once,
which is a safe nonnegative overcount.  Exact Pick counts and an outward
43,322-term log-sum-exp give

```
[-62.71954806419136023328, -62.71947138524163556114].
```

Accordingly the optimized `g=2` certificate has **62.7194713852 bits of
total security margin**, or `22.7194713852` bits beyond the requested 40.
This uses the same construction, witnesses, and exact geometric cover; only
the final union accounting changed.  Its aggregation ledger SHA-256 is
`43ecbb47c2b61e9cb20890053e879a6adf863e14923a16f4a188a9d88e2e4b40`.
The generated outward report is
`out/g2_triangle_hybrid_outward_cell_sum.json`, SHA-256
`e7245efa0c7ff6bf6cf25189f69501bbc26805ee8e7a27e5e7178d18ebebfd1d`.

The next ladder step is `g=4`.  Because the `g=2` certificate has only 0.094
bits under the obsolete global-max allocation but 22.719 bits under the
cell-local ledger, use the latter aggregation as the template for `g=4`.

### `g=4` ladder rung: diagnostic coverage and exact-cover interface

The `g=4` profile has five weight classes, `M=N/4=524288` atoms, and exactly

```
binom(M+4,4) - 717
```

feasible integer profiles after excluding total physical weights below 21
(about `2^71.4151` profiles).  Every inspected support-uniform and support-hull
profile closes pointwise.  A profile-shaped exact total-spectrum outer branch
also repaired the four initial failures of the cheap linear-BL/total-weight
atlas, leaving roughly 18800 bits at those points.

An anchor-aligned 4-dimensional Delaunay diagnostic then improved as follows:

```
atlas anchors       covered feasible volume       unresolved volume
384                    44.0427731913%                 55.9572268087%
510                    80.7020698980%                 19.2979301020%
638                    90.6674805886%                  9.3325194114%
893                    97.8192929330%                  2.1807070670%
```

The first three enrichment batches closed every selected profile.  A
volume-biased fourth batch exposed three pointwise misses among 256 profiles.
Deeper tuning closed two (`+2599` to `-7551` and `+17964` to `-4006` bits),
but the profile

```
(449970,30553,28453,10278,5034)
```

remains around `+14089` bits against the uniform `-111.415` target.  Its
current split is approximately `+184828` outer and `-170337` inner.  Full
class-by-class fugacity refinement does not close it, and the optimized
linear-BL outer coefficient is already at its valid-polytope boundary.  This
is therefore a localized outer/profile-realizability proof bottleneck, not a
numerical-margin issue.  It is not yet evidence of an actual low-weight
codeword.

Exact geometry is now separated from the diagnostic Delaunay mesh.
`probe_packet_group_g4_stellar_cover.py` emits a deterministic four-simplex
rational root triangulation and recursive five-cone stellar subdivisions.
`certify_packet_group_g4_anchor_mesh.py` checks exact determinants, facets,
ownership, rational mixture weights, source hashes, outward witness values,
and the final profile-count union factor.  A barycenter-only five-wave probe
was inefficient (8866 failed leaves), and 128 tuned stellar centroids reduced
that only to 8738.  The active refinement is to insert already-tuned atlas
profiles as arbitrary exact interior stellar points, aligning exact cells
with the witnesses without trusting floating Qhull adjacency.

No end-to-end `g=4` theorem is claimed yet.  The next gate is an exact
atlas-aligned stellar cover, followed by outward hardening; the persistent
profile above may require a sharper outer branch or an outer-feasibility
exclusion.

The geometry gate has since been resolved more cleanly.  A deterministic
rational regular triangulation keeps every exact profile anchor fixed and
uses Qhull only to propose lifted lower facets.  Every proposed retained facet
is reconstructed over `Fraction`, checked to be strictly below every other
anchor, and the projected complex is independently audited for exact
nondegeneracy, paired interior facets, clipped-hull boundary facets, and total
volume.  On the 1405-anchor enriched atlas this produces 27610 exact
4-simplices with zero degeneracies and passes every exact audit.  The parallel
facet checker gives identical deterministic output to the serial checker.

The first exact-mesh mixture pass covers 25364 of 27610 cells.  Adding 263
exact-residual witnesses raises this to 25646 cells.  A one-wave exact stellar
split confirms that simple centroid cones remain too anisotropic, so the
active strategy is exact regular retriangulation plus cell-local aggregation,
not deep stellar recursion.

The formerly persistent point

```
(375821,42364,60218,35981,9904)
```

was isolated to the old adversarial graph-hole replacement tax.  A new exact
graph/puncture total-weight lemma uses distinct-tile Maclaurin averaging,
independent coordinate bijections, exact hypergeometric band occupancy, and
the exact 24-dimensional graph-subcode spectrum.  Its fully outward result is

```
combined interval upper  -111.474649983044596...
uniform target            -111.415065016424250...
certified local margin       0.059584966620346... bits.
```

The witness artifact SHA-256 is
`2ee194712b537324ec2ef0f037d477cc61801ba741aeea265e7645d03e79352e`;
the graph spectrum SHA-256 is
`79a3f8f34280996d46ccd076c515ef1ccdc3b1f80f389ddaf2e7a24286c4f2b3`.
The canonical `exact_graph_puncture_total_weight` branch is supported by the
shared outward hardener.  Thus there is currently no known pointwise `g=4`
obstruction; the remaining gate is finite exact-cell aggregation and outward
hardening of the complete ledger.

The first global exact regular mesh nevertheless exposed a geometry defect:
some simplices joined vertices from different exact supports.  Sparse
zero-fugacity witnesses are infinite off their support, so such cells had
maximum branch bounds near `+178000` bits even though their sampled profiles
were safe.  The proof domain is now partitioned into the 31 exact positive
supports (30 feasible; the support `{0}` is empty).  On each active support
we shift `a_j=1+b_j`, build its exact integer hull, and triangulate it
independently.  These strata are disjoint on integer profiles, so their
cell-local union terms may be added without any cross-support interpolation.

For full support, the exact integer hull has 15 vertices and 9 maximal facets,
including three nontrivial kink inequalities beyond physical weight 21.  The
producer's hull-only regular mesh has 20 simplices.  An independent verifier
reconstructs the finite low-weight candidate set, checks extremality and every
supporting facet, obtains the same volume from a different 18-simplex pulling
triangulation, and then audits the producer mesh's facets, orientations, and
volume.  With 954 atlas anchors added, the first full-support mesh had 969
unique anchors and 24358 exact cells.  Every cell admitted a finite fixed
rational mixture, and 20578 cells passed the conservative uniform per-profile
target.

The full-support cell-local sum does not yet close: its first diagnostic upper
was about `+922322.53` bits, dominated by one long simplex spanning three
pure-residual corners.  This is confirmed geometry slack rather than a
pointwise obstruction.  The 128 highest-contribution cell centroids were
rounded to legal full-support integer profiles and independently tuned; every
one closed pointwise, with the weakest still about 103171 bits below its
target.  Retriangulating with those anchors reduced the worst cell branch from
about `+922259` to `+882949` bits, too slowly for repeated isotropic centroid
refinement.  The next refinement should therefore resolve the low-dimensional
pure-corner faces (or witness-dominance chambers) directly, rather than add
another undirected centroid wave.  No end-to-end `g=4` theorem is claimed yet.

That structured refinement is effective.  Adding a tuned `q=8` barycentric
grid on the recurring pure `(1,2,3)` face and rerunning the *global* exact
regular triangulation reduced the worst branch to about `+497145` bits.  A
coarse feasible `q=4` barycentric grid over all five pure directions reduced
it again to about `+301859` bits; every one of the 69 new grid profiles closed
pointwise, with at least 125199 bits of margin.  Thus the residual is still a
convex-cell localization defect, but structured grids remove it far faster
than isotropic centroid waves.

An attempted local face-star replacement was deliberately rejected by the
exact mesh audit.  Subdividing the selected triangle with edge midpoints or an
`m=3` grid creates hanging facets in neighboring simplices that contain an
edge of the triangle but not the entire triangle (the first rejected facet was
`(0,1,4,15)`).  Consequently face nodes are currently inserted as global
anchors followed by a complete exact regular retriangulation.  Any future
local refinement ledger must propagate the induced edge subdivisions through
the adjacent stars; the verifier does not accept the nonconforming shortcut.

Uniformly doubling the full 4-dimensional grid would eventually localize the
cells but is not the scalable route to `g=8`.  The active next diagnostic uses
the fact that all fixed witness bounds share the same `-log2 Q(a)` term, so
pairwise witness differences are affine.  Pairwise dominance intersections
inside the remaining high-contribution cells can therefore be used to place
targeted rational/integer discovery anchors.  These affine regions are not
the certificate themselves: the resulting global exact mesh must still assign
one fixed outward-certified witness or rational mixture to every cell.

That dominance refinement has now been exercised through several bounded
rounds.  A complete one-dimensional lower-envelope scan produced legal
transition profiles, all of which tuned pointwise; the hardest observed
transition still closed by about 1944 bits.  With 1653 full-support anchors
and 37750 exact cells, dual column generation over the entire finite witness
atlas gives the current diagnostic

```
full-support cell-local log2 upper     +67252.4553
largest cell branch                   +67198.7278
```

The mixture solver uses at most five active witnesses per cell and prices all
eligible atlas columns against the LP dual; the maximum reported final pricing
gap is below `8.5e-5` bits.  Thus finite-atlas mixture truncation is no longer
the main owner of the residual.  The full-support cover still does **not**
close, and no end-to-end `g=4` theorem is claimed.

The exact per-cell BSP route has since replaced global simplex interpolation
on the leading queue.  An exact graph/puncture-averaged outer branch removed
the dominant adversarial 128-hole tax.  Targeted multistart inner tuning then
closed all seven point profiles left by the first top-128 BSP round; the
strongest formerly residual point has about 3270.88 bits of diagnostic
pointwise margin.  Replaying the exact BSP partition initially reduced the
sole remaining positive top-128 cell (`s1fr033825`) from about `+1672.72` to
`+711.95` bits and isolated two nearby integer tuning profiles.  Independent
coordinate descent stalled there, but a joint Powell search over the four
free log-fugacity ratios and log-pole closed both profiles by more than 2270
bits.  The next exact BSP replay puts `s1fr033825` at about `-548.08` bits and
all 128 processed cells pass their uniform target.  This is a complete
diagnostic top-128 BSP frontier, not yet a global or outward theorem.

Expanding the frozen atlas to the top 512 source cells exposes the next broad
profile region.  The first pass has 656 pointwise-failure leaves, 1496 exact
rational failure-vertex occurrences, and 534 distinct rounded integer tuning
profiles.  Its largest discovery gap is about 15288.39 bits, the largest
processed-cell contribution is about `+36182.15`, and the largest unprocessed
source-cell contribution is about `+19153.80`.  Most positive cells have only
two leaves, so the next gate is a batched joint-tuning wave followed by BSP
replay rather than a new geometry construction.

The first bounded top-512 tuning wave selected the 32 largest-gap integer
representatives.  With the accelerated diagnostic tuner it completed in
288.74 seconds on 16 Peach workers: 15 profiles closed under ordinary deep
coordinate tuning and 17 did not.  The failures cluster into a few profile
families; the two hardest are about 20.0--20.5k bits short.  The next run is
therefore joint Powell tuning of those 17 profiles, seeded by the 15 newly
closed witnesses, followed by another exact BSP replay before selecting more
of the 534-profile census.

Subsequent bounded waves have moved that frontier substantially.  Replaying
the first two coordinate/joint batches reduced the top-512 diagnostic union
from about `+31591.56` to `+22142.34` bits and the largest exact pointwise gap
to about 7053.42 bits.  A second 64-profile batch (32 worst plus 32
farthest-first representatives) and a depth-12 paired Powell continuation
closed 8 of its 32 worst profiles.  Replaying all frozen witnesses then put
every processed top-512 contribution below the largest unprocessed source
cell: both the updated queue maximum and the next unprocessed ceiling were
about `+19153.80` bits.  This is a diagnostic breadth-frontier milestone, not
a theorem.

Expanding the same exact BSP discovery to the top 1024 source cells lowered
the unprocessed ceiling to about `+12078.98` bits.  The processed maximum is
initially about `+18774.53` bits, localized to the newly exposed profile
region.  That first top-1024 pass has 1116 failure leaves, 1924 exact rational
failure-vertex occurrences, and 722 distinct rounded integer tuning profiles;
the largest pointwise tuning gap is about 7280.48 bits.  The active next wave
again takes 32 worst and 32 geometrically diverse profiles, followed by paired
joint optimization and BSP replay.

That wave plus a targeted cross-profile seed bank peeled the leading
`s1fr033687` cell from about `+18584.99` to `+13144.95` bits.  Its two rounded
residual profiles now close pointwise: one by about 2046.82 bits and the other,
after three local continuations, by about 38.89 bits.  The full top-1024 replay
then moved the processed maximum to `s1fr032620` at about `+17712.39` bits.
That cell has one rounded residual profile `[444306,21517,37620,17842,3003]`.
Cross-profile Powell reduced its pointwise miss from about 7207.07 to 746.08
bits; two continuations reached 422.45 bits short.  Upgrading the same frozen
inner witness to the exact graph-conditioned total-spectrum outer recovered
122.20 bits (300.25 bits short), while reoptimizing the inner tilt under that
outer regressed to 346.62 bits short.  Further Powell depth is therefore not
the next move.  The next bounded attack should compare/tune alternate outer
families at this exact profile (including a dedicated exact-graph spectrum
multistart) or use a cell-level fixed mixture; only after that should the
top-1024 BSP be replayed again.

The proof-search hot path is now accelerated without changing the verifier.
Profiling assigned about 88% of one inner evaluation to the Python
`SharedDriveStratifiedKernel.apply` loop.  The optional native backend in
`packet_group_shared_drive_native.cpp`, loaded through `packet_group_native.py`,
implements only that fixed 65-state greedy transport recurrence; Python
remains the reference and fallback.  On Peach, 64 randomized vector checks
agreed with the Python implementation to maximum relative error
`4.32e-16`, and a 17-application kernel run improved from about `0.7002 s` to
`0.00820 s` (about 85x).  The complete four-profile deep tuner improved from
about 25 minutes to `152.6 s` (about 10x) with exactly identical reported
combined values.  Build it on Peach with

```
python scripts/build_packet_group_native.py
```

The next performance owner is fugacity-dependent histogram and point-cap
construction, not the 65-state recurrence.  Compile that table builder before
large `g=8` witness campaigns; the exact/outward certificate evaluators should
remain independent of this diagnostic acceleration.

The diagnostic point-cap builder is now native as well.  Its shared-prefix
convolution recursion agrees with the Python reference within 5 ulps for g=4
and 4 ulps for g=8 (the evaluated probability is unchanged at reported bit
precision).  On Peach it reduced the isolated cap phase from about 0.0830 s to
0.00660 s for g=4 (12.6x) and from 0.1214 s to 0.00855 s for g=8 (14.2x); a
complete representative g=4 inner objective improved from 0.280 s to 0.211 s
(1.33x).  The outward proof evaluator still uses its independent
directed-rounding implementation.

The second independent audit found a distinct finite-iteration defect in the
diagnostic search.  The Collatz test-vector bound is valid at every positive
power-iteration vector, but its finite-block objective is not monotone in the
iteration number.  The former tuner optimized with 16 iterations and scored
with the final 40-iteration vector.  For the residual profile
`[444306,21517,37620,17842,3003]`, retaining the best vector on a 320-step
trajectory selects iteration 135 and changes the combined diagnostic from
about `+188.84` to `-134.80314`.  The profile therefore has about 23.38808
bits of diagnostic margin without changing its pole or fugacities.  Freezing
that 65-entry binary64 vector and replaying the independent outward transport
gives a combined interval upper of
`-134.803121210017343850429763...`, hence a certified profile margin lower of
`23.388056193593...` bits.  The artifact is
`out/g4_cell_bsp_s1fr032620_best_collatz_outward.json`.

This closes the rounded profile, not the complete source cell.  A fresh exact
BSP discovery with the frozen vector still has four pointwise-failure leaves
and cell-local contribution about `+13840.32` bits.  Exact residual extraction
finds 18 failing rational-vertex occurrences and nine distinct rounded tuning
profiles; the largest current envelope gap is about 6931.78 bits.  Thus the
under-iterated Collatz vector was a real and substantial proof-search defect,
but it was not the sole owner of the cell-level frontier.  The next bounded
global pass must re-freeze the best trajectory vector for the atlas, rebuild
the failure census, and only then retune the surviving profiles.

The corrected top-1024 campaign now freezes the best vector from 320 Collatz
steps for every atlas witness.  It has also added 288 targeted witnesses from
coordinate and bounded joint searches.  Four greedy BSP replays improve
different source cells.  Selecting the best complete BSP tree independently
for each source cell gives a largest diagnostic contribution of about
`+16662.68` bits.  This cellwise selection is sound as a discovery operation,
but it is not yet an outward global certificate.

A greedy replay can replace a useful old partition with a worse new
partition.  The leaf-reattachment pass avoids this regression.  It preserves
every exact BSP split and minimizes each leaf bound over the expanded frozen
atlas.  Reattaching 5587 witnesses to the best-of-four geometry reduces the
largest contribution to about `+16514.51` bits.  The leading cell is
`s1fr033931`.  Its exact failing vertex rounds to
`[429359,24531,38395,28039,3964]`.

At this vertex, the corrected pointwise gap fell from about 3396.53 bits to
2200.47 bits after the targeted atlas wave.  A four-start joint search that
uses the full 320-step Collatz objective lowers the gap again to about
2017.55 bits.  The corresponding exact-rational envelope value is about
`+1906.13`, against the uniform threshold `-111.415075...`.  Thus Collatz
under-convergence and greedy partition regression both caused measurable
slack.  Neither issue explains the remaining leading gap.  The next proof
search should enrich the witness parameterization or its initialization at
this radial profile before another broad BSP replay.

### P0: Must Be Resolved Before Theorem Upgrades

- Replace BCH/random-like projection spectra with exact spectra or rigorous
  BCH-specific low-weight envelopes; the exact generic obstruction audit is
  still more than 50 bits too weak.
- Recompute the `h>2000` tail under any new BCH-like outer model.
- Keep the BCH projection rows visibly heuristic until the two items above are
  done.

### P1: Proof-Hardening And Audit Items

- For a fully formal audit packet, store raw regenerated interval artifacts or
  make the canonical verifier require the full opt-in recomputation pass under
  outward-rounded interval arithmetic.
- Replace the sampled-grid dense `0.109` exponent-gap check by an interval
  certificate if a fully formal computer-assisted proof is needed.
- Further compress the scalar dense inner ON/OFF proof chain if the main paper
  remains too long; the current narrative now marks which lemmas are local
  scaffolding versus consumed interfaces.
- Keep the bare `verify_dense_claims.py` default from being confused with the
  manuscript theorem: the current theorem check is explicitly
  `--delta 0.109`.

### P2: Presentation And Future-Construction Items

- Shorten dense outer proof machinery only if the main paper remains too long.
- Build a cleaner outer-spectrum comparison note for sliding banded, random
  block, RM, BCH-like, and possible sliding-structured outers.
- Preserve `outerExpandAcc.tex` and `innerSparse.tex` as separate construction
  lines until they are either revived or intentionally retired.

## Verification Gates

Use these as the basic post-edit checks:

```powershell
python scripts\verify_dense_claims.py --delta 0.109
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
python scripts\build_bch512_256_candidate.py
python scripts\audit_bch512_spectrum_obstruction.py
python scripts\compare_outer_modes_fullsplit.py --delta 0.09
pdflatex -interaction=nonstopmode main_permConv.tex
pdflatex -interaction=nonstopmode main_permConv.tex
```

Do not run multiple long verification or benchmark commands at the same time.
