# Merge summary: finite and asymptotic theory

Current BCH addition: `landscape_db/BCH_GROWTH_ANALYSIS.md` explains the
continuous growth trends from 212 new Q1 tuples. New sources are
`activation_q1_refresh.py`, `bch_growth_model.py`, `study_bch_growth.py`,
`verify_bch_refresh.py`, `plot_bch_growth.py`, and `test_bch_growth.py`.
The transfer retains exact uniformity after zero-input epochs. Exact
rational tests dominate an actual GF(4) encoder; a separate 90-digit replay
agrees within 1.081e-12 log units on 236 coefficients. The full local suite
passes 77 tests. No frozen producer or existing receipt changes.

At fixed t64/s20, each doubling of k costs nearly one Q1 margin bit.
Removing the row-count factor reveals distinct state-controlled small-BCH
curves and BCH-64/128 plateaus explained by first-activation geometry.
The prior four-size forecast should not be treated as reliable. The next
small task is a higher-occupation extension that handles zero syndromes,
followed by matched t/s comparisons. BCH-256 remains secondary bounded
evidence. No full-distance claim or change outside this workstream is
proposed, and higher-occupation proofs are still a dependency for parameter
selection. The figures and all generated data remain local and ignored.

GitHub integration is source-only. The finite analysis code, documentation,
and catalog configuration are integrated on the curated main history without
the local experiment commits. Generated outputs stay unchanged on disk and
are ignored. The results and validation counts below describe that local
research snapshot; a fresh checkout must generate data before running its
data-dependent reports and audits. No PR or data release is created.

The completed diagnostic grid has the 3,108 native tuples described in
`landscape_db/COMPLETE_GRID_PLAN.md`. It adds a resumable Q1 producer and
native coefficient recurrence, a kernel-aware general occupation engine,
simultaneous spectrum caps for one reused random outer, and explicit
coverage reports. Q1..64 and composition-preserving Q2/Q3/Q4 now cover
every tuple, and the typed Q65..L cover completes the occupation range.
Exact-region composition boxes give a full RM(4,9) diagnostic
of 42.969 bits at k=2^16, t64/s18. Balanced probabilities and variance spectrum caps
give three full-range random length-512 diagnostics near 60 bits; the
shared setup-event budget limits those margins. Schema version 4 keeps
conditional setup events explicit, and the union reporter checks event
containment before combining different receipts. The stronger Q2..64
whole-grid sweep is complete. The typed dense cover accounts for every
remaining occupation. Its coarse bounds and counting fallbacks
remain explicit; they do not establish distance failure. Selected outward
replay remains a separate task. See `landscape_db/GRID_FINDINGS.md` for
the final t/s comparison and `landscape_db/complete_landscape_audit.json`
for the authenticated final coverage snapshot.

The final snapshot has 263,262 observations and passes 72 tests. Its
strict completion audit verifies 13,352 dependencies and 203,739 selected
union components; all 3,108 native tuples have full evaluated coverage,
four have positive full diagnostics, and none is a new outward certificate.
The repository now includes tested
t/s operation frontiers and 473 exact-family Q1 engineering projections,
with largest-size holdout errors and window sensitivity. The bounded
BCH-256 result appears only in a secondary comparison. Verified gzip
snapshots replace tracked raw SQLite/CSV convenience files; local query
and restore commands retain the ordinary file interfaces.

Historical addition: `landscape_db/SMALL_STATE_AND_Q2_STUDY.md` records 297
smaller-state Q1 screens and 36 activation-aware Q2 screens. The changes
add an explicit three-state native pair kernel, a Windows build script,
exact-arithmetic tests, source-bound receipts, and Q1/Q2 comparison exports.
The database preserves pair weights and their witnesses. The original
tranche is committed as 88aa11a, with checkout-hash repair 9d3b7be.
The continued work remains a diagnostic study. A full test run passed 58
tests, followed by focused validation of the batched producer and shared
region ladder. All 9,541 registered receipt dependencies passed the content
and checkout-line-ending audit. No benchmark or change to another worktree
was made.

Current status: the historical RM certificate discussed below requires an
activation-state re-audit. The current parameter work is documented in
`landscape_db/ACTIVATION_PARAMETER_STUDY.md` and the final section of this
summary. Its 715 preferred observations are Q1 diagnostics, not new full
certificates. The database exposes this distinction explicitly.

## Historical fixed-RM certificate, under activation re-audit

The fixed-RM small-k proof template is now a complete outward theorem. At
message dimension (k=2^{16}) and output length (N=2^{17}), one fixed
RM(4,9) ([512,256,32]) constituent is reused in all 256 outer rows. The
structured route uses independent uniform row-coordinate and region
permutations. The inner is the fixed audited RM2Sub (t=64,s=14) map, with
one independent nonzero field scalar per epoch.

The authenticated Q=1--256 union proves

\[
  d_{\min}\ge13{,}108=\lfloor0.10N\rfloor+1
\]

except with probability below (2^{-42.5779817562755}). Q=1 is the
bottleneck. The complete outward Q=31--256 union has 1293.444358402287 bits,
so dense arithmetic is not close to the acceptance boundary.

The integration owner should treat this as a proved finite result for a
separately parameterized structured-route SPIN instance. It is not a result
for the frozen Structured SPIN (B=256, t=128, s=19) implementation, and it has
no benchmark or decoder claim. The controlled theorem statement is
`small_k_replay/RM2SUB_RM49_FINITE_CERTIFICATE.md`. The final manifest,
checker, and receipt are
`small_k_replay/RM2SUB_RM49_FULL_CERTIFICATE_MANIFEST.json`,
`small_k_replay/certify_rm2sub_rm49_full_distance.py`, and
`small_k_replay/rm2sub_rm49_t64_s14_full_distance_outward.json`.

The recommended next task is to reuse the certified engine for the planned
10% parameter frontier. Report, for each message length, both the fastest
proved exact-spectrum instance and the fastest conditional extrapolation.
Keep sub-40 margins visible, and attach XOR/runtime metadata rather than
optimizing margin alone.

## Objective and scope

Final finite update: the one-stage sparse-mixer--accumulator route is now an
unconditional finite certificate for the stated randomized variant.  At
message dimension \(K=2^{20}\) and output length \(N=2^{21}\), it proves
minimum distance at least \(228{,}590\), hence relative distance at least
10.9 percent, except with probability below \(2^{-41.3621}\).  Setup samples
one rank-tested \([512,256]\) degree-33 sparse-mixer--accumulator constituent
and reuses it at all 4,096 outer positions.  The transfer uses uniform
row-coordinate permutations, bit transpose, uniform region permutations,
and RandomStepConv with memory 22.  This is a close SPIN variant, not the
frozen Structured SPIN/RM2Sub construction.

This workstream formulated finite and asymptotic theorem targets for SPIN
codes.  It made the probability spaces, admissible-length rules, wrapper
losses, parameter schedules, and unresolved constituent obligations explicit.
It also proves a complete Random SPIN proof-model theorem at rate \(1/2\) and
relative distance \(0.11002\), using logarithmic outer blocks and logarithmic
convolution memory.  For random-block Accumulator SPIN, it proves an exact
rate--block--distance tradeoff with logarithmic blocks and one accumulator.
For a random rate-half outer paired with the fixed Structured SPIN
RM2Sub-S19 inner, it proves a complete asymptotic distance theorem at
relative distance \(0.11\). The schedule is
\(B=9\log_2N+O(1)\). Exact fixed-occupation bounds, an exact four-state
small-density certificate, and 31 outward compact-density boxes cover every
occupation from one through \(L\).

The companion linear-time audit supplies a separate structured-outer
fallback at distance \(0.101\). It uses a selected Golay--BA outer,
\(B=\Theta((\ln N)^2)\), and linear ordinary and transposed encoding work.
Its manifest, receipts, verifiers, and frozen dependencies were hash-audited
successfully from the companion worktree.

The follow-on BA--RM2Sub analysis now proves a complete scalable structured
theorem. It samples one Golay--BA-3 code and reuses it at every outer-block
position. The existing BA selection event holds with probability \(1-o(1)\),
so its failure is paid once. A 101,910-box outward verifier proves a concave
majorant of the complete BA variational spectrum. A second outward verifier
uses 621 rational-witness boxes to close every occupation density from
\(10^{-4}\) through \(1\) at relative distance \(0.11\). Exact transfer of the
four-state sparse certificate and the \(Q=1,2\) certificates closes those
regimes. A weight-coupled transfer closes every fixed \(Q\ge3\). The result
has rate \(1/2\), logarithmic outer blocks with
\(B=(39/4)\log_2N+O(1)\), and linear ordinary and transposed encoding work.
The theorem manifest binds the companion linear-time certificate and outer
certificate, the local distance receipts, and the frozen RM2Sub dependencies.

A finite pass at message dimension \(2^{20}\) compares the next three
admissible Golay--BA block sizes. The best current design point is \(B=240\),
\(L=8832\), and \(N_+=2{,}119{,}680\). After shortening to the requested
dimension, the rate is \(0.494686\). The implementation comparison selected
independent BA rows for the first proof interface: their transposed online map
was 2.699% faster than the same-binary B=256 BCH baseline. One-sided
binary64 verifiers now prove more than 47 bits of margin for occupation one
and more than 84 bits for the union of occupations 2 through 64 under the
revised \([23,217]\) conditioning window. The denser
occupations and efficient conditional setup remain open.

A separate global Expand--Convolute construction now gives the strongest
unconditional finite result. Its parent is a binary two-sided regular
expander with degrees \((10,5)\), followed by a memory-15 wrapped random
convolution. An outward first-moment certificate proves parent distance at
least 228,608 except with probability below \(2^{-50.2054082202}\). Shortening
nine input coordinates and puncturing eighteen output coordinates gives an
exact binary \([2^{21},2^{20}]\) code with distance at least 228,590, hence
relative distance at least 10.9%, with the same failure bound. Ordinary and
transposed encoding have linear work for fixed degrees and memory. This
result is not a structured
SPIN or repeated-BA theorem, has no decoder claim, and has not been
benchmarked.

The primary expander follow-up instead replaces BA locally. An exact-size
\([512,256]\) candidate uses fourteen left-regular regions and a memory-15
wrapped convolution. Its binary64 first moment gives 57.7008 bits against
nonzero kernel words. After removing zero outputs, its \(z=0.1\) weighted
moment is \(2^{-63.7288}\). At \(z=0.5\), its dense moment differs from the
random-linear value by only 0.0016 bits. The missing finite statement is a
second-factorial-moment bound for the one sampled constituent. The input pair
types are explicit, so this route avoids the unknown BCH genus-two
enumerator. The preferred iterated candidate is a degree-14 expander followed
by independently interleaved accumulators. At two accumulators, its positive
\(z=0.1\) moment is \(2^{-64.8489}\). At three accumulators, its dense
\(z=0.5\) moment is only 0.0189 bits above random. Its pair proof uses
the existing accumulator kernel and a new explicit regional-expander pair
enumerator; it has no memory-15 convolution state. The regional kernel is now
an exact coefficient identity and has been checked against direct enumeration
in 176 small cases. Composing it with the accumulator kernel completely
defines every finite second factorial moment. What remains open is an outward,
non-materialized length-512 evaluation; the full pair-type state space has
22,632,705 entries.

The finite cap gate now selects the useful layer range. Block Expand--2 has
only 32.43 bits against positive output mass outside the existing cap support
and is rejected. Block Expand--3 has 51.69 bits; stages four and five improve
that to 52.41 and 52.67 bits. Under the provisional variance target
\(\operatorname{Var}(A_w)\le2\mathbb E[A_w]\), a budget-optimized Block
Expand--5 cap event produces low/central/high band majorants below the
certified random-code values. New outward \(Q=1,2\) evaluations and inherited
\(Q\ge3\) bounds give conditional distance-failure margin 51.6440 bits.
An exact-integer regional calculation followed by five outward accumulator
transitions certifies every shell mean. Under the factor-two variance premise,
the outer cap event has 42.6283 bits and the combined result has 42.6254 bits.
Exact models at \((K,B)=(4,8),(6,12),(8,16)\) give five-stage maximum variance
factors 1.3505, 1.2397, and 1.1659. The recommended target is therefore Block
Expand--5 with variance factor two. Only the length-512 factor-two lemma
remains open, so this is not yet an unconditional theorem.

All changes are confined to `workstreams/finite_asymptotic_theory/`.  The
frozen Structured SPIN sources, receipts, naming file, roadmap, and manuscript
sources were read but not modified. Workstream-local certificate scripts were
run sequentially.

## Files created

- `NEXT_CHAT_HANDOFF.md`: self-contained continuation guide for the concrete
  Structured SPIN objective, the current sparse-mixer concentration proof,
  the certified evidence, reproduction commands, and remaining obligations.
- `FINITE_K20_ATTEMPT_LOG.md`: consolidated decision log for every finite
  `k=2^20` route tried so far, including its probability space, strongest
  result, failed proof step, implementation evidence, and present disposition.
- `EXPANDER_CODE_ALTERNATIVE_AUDIT.md`: construction-level comparison of the
  certified two-sided regular Expand--Convolute alternative with the repeated
  constituent SPIN routes.
- `BLOCK_EXPAND_CONSTITUENT_ROUTE.md`: exact local block-EC ensemble,
  first-moment diagnostics, exact regional pair kernel, finite pair-moment
  pipeline, and iterated Block Expand--\(t\) constructions.
- `verify_block_expand_pair_region_small.py` and
  `block_expand_pair_region_small_exact.json`: exact rational verification of
  the regional four-symbol coefficient formula against exhaustive edge
  assignments on small instances.
- `evaluate_block_expand_cap_budget.py`,
  `block_expand_cap_budget_diagnostic.json`,
  `evaluate_block_expand_q1_cap_transfer.py`, and
  `block_expand3_q1_cap_transfer_diagnostic.json`: finite support gate,
  recentered-cap costs, and the occupation-one check.
- `analyze_block_expand_pair_moments_small.py` and
  `block_expand_pair_moments_K4_B8_exact.json`,
  `block_expand_pair_moments_K6_B12_exact.json`, and
  `block_expand_pair_moments_K8_B16_exact.json`: exact small-model shell means,
  second factorial moments, and variances.
- `certify_block_expand5_F2_conditional_transfer.py` and
  `block_expand5_F2_conditional_transfer_outward.json`: all-occupation outward
  transfer for the frozen Block Expand caps.
- `certify_block_expand5_mean_cap_outward.py` and
  `block_expand5_mean_cap_outward.json`: outward shell means, conditional cap
  event, and combined 42.6254-bit accounting.
- `build_block_expand5_conditional_manifest.py` and
  `FINITE_K20_BLOCK_EXPAND5_F2_CONDITIONAL_MANIFEST.json`: hash-bound artifact
  inventory with the remaining variance hypothesis stated explicitly.
- `FINITE_K20_EXPANDER_EC_D10_M15_CERTIFICATE.md`: finite theorem, exact
  probability space, first-moment partition, power-of-two wrapper, and scope.
- `generate_expander_ec_k20_candidate.py`,
  `expander_ec_d10_m15_k20_parent_candidate.json`,
  `certify_expander_ec_k20_wrapper.py`, and
  `expander_ec_d10_m15_k20_wrapper_outward.json`: frozen parent candidate,
  independent structural check, outward verification, and wrapper receipt.
- `build_expander_ec_k20_manifest.py` and
  `FINITE_K20_EXPANDER_EC_D10_M15_MANIFEST.json`: hash-bound artifact and
  external-verifier dependency manifest for the finite expander theorem.
- `FINITE_LENGTH_FRAMEWORK.md`: exact weight-slice and type-orbit first-moment
  theorems, existence logic, probability spaces, and certificate interface.
- `ARBITRARY_LENGTH_WRAPPER.md`: admissibility predicate, padding,
  puncturing, shortening, adjacent-block option, and a single exact-length
  theorem.
- `ASYMPTOTIC_SCALING.md`: schedules for Accumulator SPIN, Random SPIN, and
  Structured SPIN, with sparse and bulk margin scales.
- `ACCUMULATOR_SPIN_ASYMPTOTIC.md`: sharp hypergeometric contraction, exact
  random-block theorem, admissible-length rule, and explicit rate-half
  constants for one accumulator.
- `BIT_TRANSPOSE_ACCUMULATOR_OBSTRUCTION.md`: exact adjacent-pair obstruction
  for replacing the uniform interleaver by a pure bit transpose.
- `REGION_SHUFFLED_TRANSPOSE_ACCUMULATOR.md`: exact two-state transfer for a
  transpose followed by independent region permutations, plus fixed-support
  asymptotic saddles.
- `STRUCTURED_TRANSLATION_FROM_RANDOM_BASELINE.md`: exact symmetrized
  spectrum-to-region-profile identity, a nonnegative-transfer comparison with
  the random outer, and the resulting logarithmic-block criterion.
- `RANDOM_OUTER_RM2SUB_RAMP.md`: minimal scalable random-outer ensemble with
  the fixed RM2Sub-S19 inner, inner bijectivity, and the exact conditional
  one-active transfer theorem.
- `RM2SUB_ONE_ACTIVE_CONTINUUM.md`: exact full-state epoch kernel,
  long-region two-state limit, occupation-one spectral exponent, optimized
  logarithmic-block threshold, and an exact rational \(c=18/5\) certificate.
- `RM2SUB_TWO_ACTIVE_CONTINUUM.md`: exact two-position region law,
  same-epoch collision term, two-impulse limiting matrix, occupation-two
  exponent, and a second exact \(c=18/5\) certificate.
- `RM2SUB_FIXED_OCCUPATION_CONTINUUM.md`: ordered reset integrals,
  Dirichlet/beta coefficient formula, and the exact common-tilt theorem for
  every fixed outer occupation.
- `RM2SUB_UNIFORM_FIXED_OCCUPATION.md`: beta-integral domination, a uniform
  Perron bound for all fixed \(Q\ge3\), and the exact rational \(c=9\)
  certificate at relative distance \(0.11\).
- `RM2SUB_GROWING_SPARSE_LIFT.md`: exact nonzero-row relaxation, placement
  collision estimates, a uniform zero-epoch moment bound, and the precise
  finite-kernel lemma needed for \(Q=Q_N\to\infty\).
- `RM2SUB_DENSE_OCCUPATION.md`: finite Bernoulli conditioning, exact
  constituent moments, rigorous three- and four-state transfers, the exact
  small-density certificate, and the outward compact-density certificate.
- `RANDOM_OUTER_RM2SUB_CERTIFICATE.md`: complete rate-half distance-\(0.11\)
  theorem for the random outer with \(B=9\log_2N+O(1)\).
- `STRUCTURED_CERTIFICATE_CROSSCHECK.md`: hash audit and proof-scope comparison
  for the companion structured \(0.101\), linear-time certificate.
- `SINGLE_SAMPLED_BA_RM2SUB_D11.md`: complete one-sampled
  Golay--BA-3/RM2Sub-S19 theorem at rate \(1/2\), relative distance \(0.11\),
  \(B=(39/4)\log_2N+O(1)\), and linear ordinary and transposed work.
- `certify_golay_ba_concave_majorant.py` and
  `golay_ba3_concave_majorant.json`: outward affine-support certificate for a
  concave majorant of the BA variational spectrum.
- `certify_golay_ba_rm2sub_joint_interval.py` and
  `golay_ba3_rm2sub_joint_interval_d11.json`: outward positive-occupation
  certificate with rational change-of-measure and Collatz witnesses.
- `certify_golay_ba_rm2sub_sparse.py` and
  `golay_ba3_rm2sub_sparse_d11.json`: exact transfer of the fixed-occupation
  and four-state sparse certificates, including the block constant \(39/4\).
- `WEIGHT_COUPLED_FIXED_OCCUPATION.md`,
  `certify_golay_ba_rm2sub_weight_coupled_fixed.py`, and
  `golay_ba3_rm2sub_weight_coupled_fixed_d11.json`: row-weight fugacity
  transfer and 62 outward endpoint checks for every fixed \(Q\ge3\).
- `FINITE_K20_COMPARISON.md`: finite \(k=2^{20}\) comparison of admissible
  \(B=216,240,264\) instances, explicit conditioned probability spaces,
  online work and setup storage, and the remaining 40-bit obligations.
- `FINITE_K20_INDEPENDENT_SETUP.md`: frozen mathematical setup law for the
  independent-row B=240 target, shortening convention, recurrence, proved
  partial bounds, performance evidence, and open setup interface.
- `FINITE_K20_PROOF_OBLIGATIONS.md`: active ledger for the exact 40-bit
  theorem target, occupation coverage, RM2Sub audit, outward arithmetic,
  certificate packaging, and performance limits.
- `FINITE_K20_DENSE_TRANSFER_TARGET.md`: exact weight-band mixture reduction
  for occupations 65 through 8,832, the multitype finite RM2Sub interface, and
  the two acceptable positive-compression lemmas still to prove.
- `certify_golay_ba_rm2sub_finite_one_active.py` and
  `golay_ba3_rm2sub_finite_B240_q1_outward_k20_d11.json`: one-sided
  occupation-one certificate covering every permitted BA weight and proving
  a contribution below \(2^{-47}\).
- `certify_golay_ba_rm2sub_finite_q2_64.py` and
  `golay_ba3_rm2sub_finite_B240_q2_64_outward_k20_d11.json`: one-sided
  Bernoulli-envelope certificate for every occupation 2 through 64, with an
  aggregate contribution below \(2^{-89}\).
- `build_finite_k20_partial_manifest.py` and
  `FINITE_K20_PARTIAL_CERTIFICATE_MANIFEST.json`: reproducibility commands,
  hashes, proved occupation ranges, and an explicit list of claims that remain
  open. The manifest is intentionally partial and does not claim the final
  40-bit theorem.
- `evaluate_golay_ba_rm2sub_finite.py` and the
  `golay_ba3_rm2sub_finite_*_k20_d11.json` receipts: exact-in-form binary64
  BA spectra and one-active RM2Sub evaluations at the three block sizes.
- `golay_ba3_B*_conditioned_spectrum_upper.json`: pointwise conditional
  expected-spectrum bounds obtained from the finite Markov event.
- `evaluate_golay_ba_rm2sub_finite_holder.py` and its corrected Holder
  receipt: full-spectrum dense-occupation diagnostic. The corrected bound is
  still vacuous at some large occupations and is retained as an open-boundary
  record, not as evidence against the construction.
- `build_single_sampled_ba_manifest.py` and
  `SINGLE_SAMPLED_BA_RM2SUB_CERTIFICATE_MANIFEST.json`: SHA-256 binding for
  the new theorem, verifiers, receipts, companion outer proof, and frozen
  dependencies.
- `analyze_golay_ba_rm2sub_joint.py`: exact finite BA-spectrum and
  variational-asymptotic binary64 diagnostic for the joint exponent.
- `golay_ba3_rm2sub_joint_B480_d11.json`,
  `golay_ba3_rm2sub_joint_B960_d11.json`, and
  `golay_ba3_rm2sub_joint_asymptotic_d11.json`: diagnostic receipts for the
  corrected value-tilted saddle.
- `golay_ba3_rm2sub_joint_expurgated_hull_d11.json`: mixed-type diagnostic
  using the proposed expurgated-random concave hull.
- `golay_ba3_rm2sub_joint_expurgated_hull_sparse_d11.json`: small-occupation
  diagnostic showing where the three-state puncturing loss requires the
  four-state transfer.
- `golay_ba3_rm2sub_joint_expurgated_hull_sparse_four_state_d11.json`:
  four-state small-density diagnostic supporting \(c=6\) as the next
  logarithmic-block certificate target.
- `RANDOM_OUTER_RM2SUB_CERTIFICATE_MANIFEST.json`: SHA-256 binding for the
  random-outer theorem, fixed- and dense-occupation verifiers, receipts, and
  frozen RM2Sub selection.
- `STRUCTURED_SPIN_THEOREM_TARGET.md`: finite and asymptotic conditional
  theorems for structured orbits, plus obligations for the frozen
  constituents.
- `BA_ASYMPTOTIC_EXPLORATION.md`: exact repeated-constituent BA ensemble,
  a proved sparse-path obstruction, polynomial block schedules, EBCH128
  diagnostics, and the next theorem target.
- `RANDOM_SPIN_PROOF_AUDIT.md`: exact random-block outer law,
   logarithmic schedule, audit of the old random-convolution proof, an exact
   transfer matrix, parameterized sparse sufficient conditions, and a complete
   proof-model asymptotic theorem at relative distance \(0.11002\).
- `analyze_ba_asymptotics.py`: sequential binary64 diagnostic for comparing
  exact BA weighted spectra with their sparse boundary paths.
- `evaluate_random_spin_finite.py`: exact-in-form binary64 first-moment
  evaluator with an independent exact integer validation.
- `certify_random_spin_linear_exponent.py`: interval check for the endpoint
  inequalities and schedule margins in the Random SPIN proof.
- `certify_accumulator_spin_asymptotic.py`: interval check for a strict
  Accumulator SPIN rate, distance, and block-constant triple.
- `analyze_bit_transpose_accumulator_obstruction.py`: exact random-subspace
  intersection calculation for a pure transpose.
- `analyze_region_shuffled_transpose_accumulator.py`: binary64 evaluation of
  the fixed-active-block continuum saddle for region shuffles.
- `compare_structured_outer_to_random.py`: binary64 diagnostic for the
  pointwise excess of a local expected spectrum over the uniform-random
  injection spectrum.
- `evaluate_random_outer_rm2sub_one_active.py`: binary64 one-active evaluator
  for a uniform random outer paired with the frozen RM2Sub-S19 transfer.
- `derive_rm2sub_one_active_continuum.py`: binary64 optimizer for the proved
  continuum formula and exact rational checker for the \(c=18/5\) choice.
- `derive_rm2sub_two_active_continuum.py`: binary64 optimizer for the
  occupation-two formula and exact rational checker for the same schedule.
- `analyze_rm2sub_fixed_occupation_continuum.py`: sequential binary64
  evaluator for the exact fixed-occupation continuum formula.
- `certify_rm2sub_uniform_fixed_occupation.py`: exact rational checker for
  the all-fixed-\(Q\ge3\), \(c=9\) certificate.
- `analyze_rm2sub_dense_occupation.py`: exact integer association audit and
  binary64 witness search for the dense transfer.
- `certify_rm2sub_dense_small.py`: exact rational Bernstein certificate for
  every \(0<Q/L\le10^{-4}\).
- `certify_rm2sub_dense_compact_interval.py`: 100-digit outward interval
  verifier for every \(10^{-4}\le Q/L\le1\).
- `random_outer_rm2sub_one_active_B256_d11.json`: diagnostic output at the
  frozen dimensions, including a same-grid structured-spectrum comparison.
- `rm2sub_one_active_continuum.json`: continuum constants and the exact
  rational inequality receipt.
- `rm2sub_two_active_continuum.json`: occupation-two constants, limiting
  matrices, and the exact rational inequality receipt.
- `rm2sub_fixed_occupation_continuum_d11.json`: optimized common-tilt
  diagnostics through sampled occupation \(Q=512\).
- `rm2sub_uniform_fixed_occupation_d11.json`: exact rational bounds for the
  uniform fixed-occupation certificate.
- `rm2sub_dense_small_exact_d11.json`: exact small-density Bernstein receipt.
- `rm2sub_dense_compact_interval_d11.json`: 31-box outward compact-density
  receipt.
- `optimize_random_spin_asymptotic_constants.py`: binary64 diagnostic search
  for simple sparse-saddle, memory, and block constants; selected candidates
  still require the interval checker.
- `MERGE_SUMMARY.md`: this handoff.

`BRIEF.md` was not changed.

## Decisions recommended to the integration owner

1. Use the exact event
   \[
     Z_D=|\{x\ne0:\operatorname{wt}(E(x))<D\}|
   \]
   so that \(Z_D=0\) means injectivity and minimum distance at least \(D\).
2. State the finite theorem first for a fixed injective outer map and a
   uniform global interleaver.  State random-outer factorization only under
   explicit independence.
3. Use the type-orbit theorem for Structured SPIN.  Do not substitute a
   Hamming-slice transfer probability unless slice uniformity is proved.
4. Make each construction manifest publish its mathematical admissible-length
   set.  Keep implementation tiling constraints separate.
5. Use next-admissible padding plus fixed-coordinate puncturing as the default
   arbitrary-length wrapper.  Treat adjacent block sizes as an optional
   construction theorem, not as free rounding.
6. Report finite results as \((K,N,D,\mu_D)\), where \(\mu_D\) is an outward
   upper bound on the expected number of bad nonzero codewords.
7. Use \(B_N=\Theta(\log N)\) when the one-block bad mass is exponentially
   small in \(B_N\).  If it is only \(B_N^{-p+o(1)}\), use a sublinear
   schedule \(B_N=N^{\beta+o(1)}\) with \(\beta>1/(p+1)\), unless a finer
   structured-inner argument proves stronger suppression.
8. State the expected global margin as \(\Theta(\log N)\) and the many-block
   margin as \(\Theta(N)\), under the corresponding uniform hypotheses.
9. Keep the frozen Structured SPIN (B=256, t=128, s=19) statement finite.
   Any scalable constituent family is a family-level extension, not a change
   to the frozen instance.
10. Resolve the naming/construction mismatch before paper integration:
    `innerAcc.tex` and the current dense integration use a random sliding dense
    outer, while `SPIN_NAMING.md` defines Accumulator SPIN and Random SPIN with
    random block outers.
11. For Structured SPIN, benchmark a candidate outer against a random outer
    paired with the same structured inner.  Do not import the
    random-convolution exponent into RM2Sub without proving an RM2Sub transfer
    comparison.
12. Use the fixed \(t=128,s=19\) RM2Sub constituent in the first scalable
    intermediate model.  Let the even outer length \(B_N\) grow while
    requiring \(128\mid L_N\).  Scaling \(t\) with \(B_N\) would require a
    new constituent family and should be deferred.

The integration owner should also normalize the notation in `framework.tex`.
That file first defines `A_h` as an expected spectrum, then later takes an
expectation of `A_h` again in the random-outer corollary.  The finite framework
deliverable instead uses a realized count for a fixed outer map and introduces
the expectation exactly once.  Its good-event corollary also retains the
indicator of the outer good event, avoiding an implicit conditioning step.

## Claims established and evidence

### Proved within this workstream

The following claims follow from elementary linearity, counting, or standard
code operations.  Complete arguments appear in the deliverables.

- For a fixed injective outer code, uniform interleaver, and independent inner
  map,
  \[
    \mathbb E[Z_D]=\sum_h A_h^{\mathrm{out}}p_h^{\mathrm{in}}(D).
  \]
- The same expectation upper-bounds setup failure by Markov's inequality.
- The exact type-orbit identity remains valid for an arbitrary joint setup
  distribution.
- Input padding preserves the certified distance.
- Puncturing \(q\) fixed coordinates loses at most \(q\) distance and remains
  injective when \(q<D\).
- A bounded-gap admissible family extends to all lengths with an
  \(o(1)\) relative rate and distance loss.
- Under the stated sparse product bound, \(B_N=c\log_2N\) gives a vanishing
  one-block union bound when \(c>(a+1)/\lambda\).
- For a uniform weight-\(h\) input to one accumulator, the exact low-output
  probability is a hypergeometric majority tail.  Consequently,
  \[
    p_{N,h}(\lfloor\delta N\rfloor)
    \le[2\sqrt{\delta(1-\delta)}]^h
  \]
  for every \(N\), \(h\), and \(0<\delta<1/2\).
- For random rate-\(R\) block injections followed by one global permutation
  and one accumulator, define
  \[
    \lambda_A(R,\delta)
    :=1-R-\log_2(1+2\sqrt{\delta(1-\delta)}).
  \]
  The ensemble has relative distance greater than \(\delta\) with setup
  failure probability \(o(1)\) whenever \(\lambda_A>0\) and
  \(B_N\ge c\log_2N\) with \(c\lambda_A\ge1\).
- At rate one half, every \(c>2\) yields a positive distance.  The certified
  ceilings are \(0.0037634\) for \(c=3\), \(0.0225363\) for \(c=8\), and
  \(0.0330838\) for \(c=17\).  The ceiling tends to \(0.0449101\) as
  \(c\to\infty\).
- A pure bit transpose fails at rate one half.  If adjacent outer block
  subspaces share a nonzero word, the two equal rows produce accumulator
  output weight at most \(B\).  Two independent random half-dimensional
  subspaces intersect nontrivially with probability at least \(1/3\).  Hence
  \[
    \Pr[d_{\min}>B]\le(2/3)^{\lfloor L/2\rfloor}.
  \]
  Every pure-transpose family with \(B_N=o(N)\) has vanishing relative
  distance in probability.
- For transpose plus independent region permutations, conditioning on \(q\)
  active blocks gives an exact two-state kernel per region.  A binomial mixture
  over the number of region ones, followed by a \(B\)-th matrix power, gives a
  finite first-moment upper bound.  The derivation includes the explicit cost
  of relaxing active rows from uniform nonzero words to uniform binary words.
  Per-block coordinate permutations are omitted because the random outer law
  already makes every active row uniform nonzero.
- For a structured local outer with expected spectrum
  \(\overline A_{B,h}\), independent uniform coordinate permutations give the
  symmetrized profile polynomial
  \[
    S_B(\boldsymbol y)
    =
    \sum_h\frac{\overline A_{B,h}}{\binom Bh}
    e_h(\boldsymbol y).
  \]
  For \(q\) independently sampled active blocks, the coefficient of
  \(\boldsymbol y^{\boldsymbol a}\) in \(S_B^q\) is the expected number of
  messages with transposed region-count profile \(\boldsymbol a\).  Combining
  these coefficients with the conditional region kernels gives an exact
  finite first-moment bound for a region-shuffled transpose followed by one
  accumulator.
- Let \(g_B\) be the maximum log-ratio of the structured expected local
  spectrum to the uniform-random injection spectrum.  Every nonnegative
  region-profile transfer at occupation \(q\) is at most \(2^{g_Bq}\) times
  its random-outer counterpart.  Hence \(g_B=o(B)\) preserves the limiting
  sparse exponent and logarithmic block constant, provided the reference
  bound uses the same inner kernel.
- For every fixed multiplier schedule, the RM2Sub recurrence
  \[
    Y_r=X_r+A(Q_r),
    \qquad
    Q_{r+1}=\alpha_rQ_r+C(X_r)
  \]
  is bijective as a map from the complete input sequence to the complete
  output sequence.  Sequential inversion recovers \(X_r\) from \(Y_r\) and
  the already recovered state.
- For the all-active random-outer class, relaxing every nonzero outer row to a
  uniform binary row makes the complete inner input uniform.  Inner
  bijectivity preserves that law.  Consequently
  \[
    \mathbb E[Z_{D,L}]
    \le
    \left(\frac{2^{K_B}-1}{1-2^{-B}}\right)^L
    2^{-N}\sum_{w=0}^{D}\binom Nw.
  \]
  At rate one half, this class has exponent
  \(H_2(\delta)-1/2+o(1)\).
- For one active random outer block, the region permutations reduce each
  active region to one uniformly placed impulse.  If \(W_0,W_1\) are valid
  two-state epoch envelopes, the inactive and active region matrices are
  \[
    R_0=W_0^{L/t},
    \qquad
    R_1=\frac{t}{L}\sum_r W_0^rW_1W_0^{L/t-1-r}.
  \]
  Coefficient extraction from \((R_0+uR_1)^B\), followed by the exact random
  nonzero-row density, gives the one-active first-moment bound in
  RANDOM_OUTER_RM2SUB_RAMP.md.
- The finite two-state RM2Sub envelope is not asymptotically iterable: it
  pays the punctured-live factor once per zero epoch. Using the exact
  \(2^{19}\)-state epoch kernel first and then taking \(L\to\infty\) gives
  \[
    K_0=\begin{pmatrix}1&0\\0&a\end{pmatrix},
    \qquad
    K_1=\begin{pmatrix}0&f\\f/M&(1-1/M)a\end{pmatrix},
  \]
  where \(M=2^{19}-1\), \(a=e^{-p\theta}\),
  \(f=(1-e^{-p\theta})/(p\theta)\), and \(p=2^{18}/M\). The Perron root of
  \(K_0+K_1\) gives
  \[
    \log_2\mathbb E Z_{D,1}
    \le \log_2L-\eta_{19}(\delta)B+o(B).
  \]
  At \(\delta=0.11002\), binary64 optimization gives
  \(\eta_{19}=0.2787705115\ldots\) and threshold
  \(c>3.587179987\ldots\). An exact rational inequality certifies
  \(c=18/5\) and polynomial decay exponent greater than
  \(0.0026858604\). This proves occupation one only.
- With two active outer blocks, a region receives zero, one, or two impulses.
  Two impulses occupy distinct bit positions, and their exact same-epoch
  collision probability is \(127/(L-1)\). The collision remains in the
  finite transfer and vanishes in the continuum. The two-impulse limit,
  combined with \(K_0+2K_1+K_2\), gives the occupation-two threshold
  \(c>3.587809926\ldots\) at \(\delta=0.11002\). An exact rational inequality
  certifies \(c=18/5\) and decay exponent greater than \(0.0049419781\).
  The optimized occupation-two threshold is slightly larger than the
  occupation-one threshold.
- For every fixed occupation \(Q\), the region matrix has the ordered-integral
  representation
  \[
    T_Q(\theta)=\sum_{a=0}^Q\binom Qa K_a(\theta).
  \]
  A Dirichlet-spacing decomposition converts this matrix into nonnegative
  path coefficients multiplied by beta Laplace transforms. It proves
  \[
    \log_2\mathbb E Z_{\lfloor\delta N\rfloor,Q}
    \le Q\log_2L-\eta_Q(\delta)B+o_Q(B)
  \]
  for every fixed \(Q\). At \(\delta=0.11\), the binary64 common-tilt
  threshold first exceeds \(c=18/5\) at sampled \(Q=16\). It reaches
  \(4.10383\) at \(Q=128\) and \(4.55763\) at \(Q=512\). This diagnostic
  does not prove a distance obstruction; it shows that the existing small-\(Q\)
  constant is not a uniform certificate.
- The beta-transform admits a nonasymptotic coefficient domination. For
  \(y\in(0,1]\), \(\tau>0\), and
  \[
    G_Q(y,\tau):=
    \sup_{0\le x\le1}
    \left\{-\tau x+
    \left(1+\frac1Q\right)\ln(1-x+x/y)\right\},
  \]
  it gives
  \[
    \lambda_Q(Q\tau/p)
    \le
    \left(e^{G_Q(y,\tau)}\rho(H D_y)\right)^Q.
  \]
  Choosing \(y=1/2,\tau=11/10\) bounds all fixed \(Q\ge3\) by the
  \(Q=3\) case. Exact rational arithmetic certifies
  \[
    \eta_Q(0.11)/Q>0.11604381273505275.
  \]
  Hence \(B=9\log_2N+O(1)\) gives decay exponent greater than
  \(0.04439431461547Q\) for every fixed \(Q\ge3\). Together with the sharper
  occupation-one and occupation-two receipts, this closes every fixed
  occupation at relative distance \(0.11\).
- For growing occupation, conditioning \(Q\) uniform binary rows to be
  nonzero costs exactly at most \((1-2^{-B})^{-Q}\). A fixed pair of active
  rows lands in the same or adjacent RM2Sub epoch with probability at most
  \(383/(L-1)\). Hoeffding's lemma bounds the cumulative zero-epoch
  discretization cost by
  \[
    \exp\!\left(\frac{Q^2\tau^2t}{8p^2L}\right)
  \]
  per region. These proved estimates all have logarithmic scale
  \(O(Q^2/L)\). They remain a valid but superseded cluster-comparison route.
- The completed finite lift instead conditions an iid Bernoulli candidate
  process on exactly \(Q\) marks per region. An exact four-state transfer
  removes the repeated puncturing loss at vanishing density. A three-state
  transfer is sharper on compact positive densities. These two transfers
  close every growing occupation without assuming \(Q^2B/L=o(1)\).
- For a repeated fixed constituent of minimum distance \(d_0\), followed by
  \(\ell\) independent uniform-interleaver accumulator stages, the minimum
  sparse path \(h_0=d_0\), \(h_j=\lceil h_{j-1}/2\rceil\) contributes
  \(\Omega(B^{-p_\ell(d_0)})\) to every fixed-\(\rho\) weighted spectrum,
  where \(p_\ell(d_0)=\sum_{j=1}^{\ell}h_j-1\).  Thus this weight-only
  interface cannot satisfy an exponentially small one-block bound.
- For independent uniform rate-half block injections,
  \[
    G_B(z)
    =
    (2^{B/2}-1)\frac{(1+z)^B-1}{2^B-1},
  \]
  and the complete expected generating function is
  \((1+G_B(z))^{N/B}\).  For every \(z<\sqrt2-1\), a sufficiently large
  logarithmic block constant makes its nonzero part \(o(1)\).
- For the non-wrapping random convolution and a uniform weight-\(h\)
  interleaved input, the pair consisting of remaining input weight and clipped
  trailing-zero count is an exact time-inhomogeneous Markov chain.  Its
  transition probabilities are given in RANDOM_SPIN_PROOF_AUDIT.md.
- The same convolution has an exact \((m+1)\)-state bivariate transfer matrix
  that marks input and output weight.  Coefficient extraction gives the exact
  output-weight transform on every Hamming slice, and a nonnegative-
  coefficient Chernoff bound reduces the distance tail to a sparse matrix
  power.
- In the sparse saddle regime, the transfer matrix Perron root obeys
  \(\lambda_m-1=O(b^m)\), where \(b=(1+a)/2<1\).  Consequently,
  \(m=\lceil\gamma\log_2N\rceil\) suppresses its exponential perturbation when
  \(\gamma>1/\log_2(1/b)\).
- The limiting Perron root is
  \[
    \max\{1,(1+a)(1+s)/2\}.
  \]
  Its optimized coefficient exponent closes every fixed linear input-weight
  class at relative distance \(0.11002\).  Interval arithmetic certifies a
  maximum bulk exponent below \(-2.3718\times10^{-5}\).
- With \(m_N=\lceil2\log_2N\rceil\), the convolution satisfies
  \[
    p_{N,h}(0.09)\le C N^4\sqrt h\,(0.285750)^h
    \qquad (1\le h\le N/25).
  \]
  Because \(0.285750<\sqrt2-1\), the exact outer generating function closes
  this range when \(B_N\ge42\log_2N\).
  This remains a proved conservative baseline.
- The sparse argument is parameterized by cutoff \(\eta_*\), interior saddle
  \(\chi\), memory constant \(\gamma\), and block constant \(c_R\).  The
  exact sufficient inequalities appear as (36)--(38) in
  RANDOM_SPIN_PROOF_AUDIT.md.
- For
  \((\eta_*,\chi,\gamma,c_R)=(1/5000,101/100,51/50,17)\), interval arithmetic
  certifies a memory margin greater than \(0.0186697\), a block-schedule margin
  greater than \(0.0576243\), and sparse contraction
  \(\rho_*<0.221264<\sqrt2-1\).
- On the optimized admissible family, the proof-model Random SPIN ensemble has
  rate \(1/2\), relative distance greater than \(0.11002\), and setup failure
  probability \(o(1)\).  Both the outer block and convolution memory are
  logarithmic.  The admissible gaps are \(o(N)\), so the arbitrary-length
  wrapper preserves the asymptotic rate and distance.
- For the random-outer RM2Sub family, the exact four-state Collatz witness
  proves a natural-log exponent at most
  \(-780897(QB)/6665600\) when \(0<Q/L\le10^{-4}\), before the support
  choice. The schedule \(B=9\log_2N+O(1)\) pays that support cost with a
  strict margin greater than \(0.04013QB\).
- Thirty-one rational-witness boxes, evaluated with 100-digit outward
  intervals, cover every \(10^{-4}\le Q/L\le1\). The largest accepted
  exponent upper bound is \(-4.3689\times10^{-7}\) natural units per output
  bit.
- Combining these bounds with the fixed-occupation continuum proves rate one
  half and asymptotic relative distance \(0.11\) for the complete
  random-outer RM2Sub family.
- One sampled Golay--BA-3 outer satisfies the existing pointwise good-spectrum
  event with probability \(1-o(1)\) and may be reused across all outer-block
  positions.
- After row-specific Bernoulli conditioning, a uniform region permutation
  dominates the iid Bernoulli law with the same mean by at most \(L+1\).
  Across \(B\) regions this costs \(\exp(o(N))\). This lemma and the proved
  concave spectrum majorant reduce arbitrary row-weight mixtures to their
  mean without an exponential heterogeneity penalty.
- The outward BA majorant and joint RM2Sub interval receipts, together with
  the exact sparse transfer, prove rate \(1/2\), distance \(0.11\),
  \(B=(39/4)\log_2N+O(1)\), and linear ordinary and transposed work for the
  scalable one-sampled Golay--BA-3/RM2Sub-S19 family.

### Proved facts imported from the source material

- A uniform permutation makes a fixed weight word uniform on its Hamming
  slice (`framework.tex`).
- The accumulator input-output enumerator and its simple low-output
  contraction are derived exactly (`innerAcc.tex`).
- The random dense inner material derives explicit finite transfer envelopes
  and analytic exponent functions (`innerDenseScalar.tex`).

These imported facts were used to formulate targets.  This workstream did not
independently replay their algebra or numerical checks.

### Conditional statements

- The Accumulator SPIN existence theorem is unconditional for its stated
  random-block proof-model ensemble.  A finite outward-rounded certificate
  remains open.
- For the region-shuffled transpose, every fixed number \(q\) of active blocks
  closes under its stated saddle condition.  A complete theorem requires a
  bound uniform for \(q\to\infty\) and a separate linear-\(q\) argument.
- The Random SPIN existence theorem is unconditional for its stated
  proof-model ensemble.  Smaller constants, larger distance, and a finite
  outward-rounded certificate remain optimization tasks.
- The Structured SPIN asymptotic theorem follows from the sparse-profile
  product bound (S) and many-block exponent bound (B).
- Arbitrary-length asymptotics require bounded gaps between admissible lengths.

### Diagnostics only

At \(k=2^{20}\), the nearest-binary64 finite evaluator gives conditioned
one-active margins of 41.780, 48.019, and 54.795 bits for
\(B=216,240,264\), respectively. The corresponding output-length overheads
over \(2k\) are 0.195, 1.074, and 1.514 percent. The B=240 values for
occupations 1 through 64 are no longer merely diagnostic: outward verifiers
prove a combined margin above 47 bits. The dense range is not closed by the
current worst-shell or Holder comparisons.

The frozen Structured SPIN receipts report nearest-binary64 margins of about
55.8646 bits at relative distance 0.11, 61.0696 bits at 0.10, and 66.2701 bits
at 0.09.  The one-active class is limiting.  These values remain diagnostics
because the outer spectrum is modeled and the arithmetic is not outward
rounded.

The modeled ParityFanout-31x33 expected spectrum at \((B,K_B)=(256,128)\)
has pointwise excess
\[
  g_{256}=0.0065094697569\text{ bits}
\]
over the exact uniform-random injection spectrum.  The maximum is at weight
241.  The excess is only \(2.83681\times10^{-6}\) bits on weights 37 through
219 and is at binary64 roundoff scale on weights 64 through 192.  This
explains the modeled Bernoulli-envelope cost of \(128.0065094697569\) bits per
active block, but it is not an authenticated spectrum statement.

With the same frozen RM2Sub-S19 entrywise transfer and the same \(0.01\)
log-surprisal grid, a uniform random \([256,128]\) outer gives a one-active
margin of 54.879111670 bits at relative distance \(0.11\).  The modeled
ParityFanout-31x33 spectrum gives 55.864654977 bits.  The modeled structured
outer is therefore 0.985543307 bits better in this class.  The random-outer
common-tilt bound gives 51.876707694 bits.  These are binary64 diagnostics,
not outward-rounded certificates.
The exact-in-form all-active relaxation gives 188.379385174 bits of margin
at the same frozen length and distance.  Its only numerical operation is the
binary64 evaluation of the Hamming-ball sum.

The Random SPIN `0.109` checkpoint includes a sampled-grid exponent check.  A
sampled grid does not prove the required uniform inequality.

For the region-shuffled transpose at rate one half and distance \(0.02\), the
fixed-support diagnostic requires \(c=3.2835\) at \(q=1\), \(6.2196\) at
\(q=8\), and \(6.9265\) at \(q=256\).  The values approach the
uniform-interleaver constant \(6.9514\).  This convergence is numerical and
does not close the growing-\(q\) proof obligation.

The new exact-in-form binary64 evaluator reports
\(\log_2\mu=-4.30265\) at \((N,B,m,D)=(512,32,24,46)\) and
\(\log_2\mu=-3.28969\) at \((1024,32,28,92)\).  The control point with
\((N,B,m,D)=(512,16,24,46)\) gives \(\log_2\mu=2.02265\).  These are
diagnostics, not outward-rounded finite certificates.

For the authenticated EBCH \([128,64,22]\) spectrum followed by two
accumulators, the binary64 weighted-spectrum diagnostic at
\(\rho\in\{0.05,0.1,0.2\}\) has tail log-log slopes close to \(-16\), matching
the proved sparse-path lower-bound exponent.  This supports, but does not
prove, a matching \(O(B^{-16})\) upper bound.  The EBCH32 and shortened-XBCH64
spectrum filenames referenced by older receipts are absent, so they were not
silently reconstructed.

A legacy packet-transpose/random-step construction reports a complete
outward certificate with about 195.2839 bits of margin at distance 188744 for
its declared finite ensemble.  Its factored packet interleaver and relaxed
outer model differ from the canonical Random SPIN definition.  The separate
BCH/bit-shuffle RandomStepConv endpoint remains a modeled-spectrum,
nearest-binary64 diagnostic.

## Assumptions and open issues

### Framework

- The random-outer weight factorization needs an ensemble supported on
  injective outer maps and independence from the interleaver and inner map.
- If setup components share randomness, the theorem must retain their joint
  expectation or condition on the shared setup.
- A structured interleaver needs a complete orbit type and a proof that the
  chosen type is sufficient for the transfer law.

### Frozen Structured SPIN

- Prove a spectrum or sufficient weighted envelope for the actual
  extended-BCH-based outer after ParityFanout-31x33.
- Define the setup distributions behind the source seeds and prove the needed
  independence and uniformity statements.
- Prove the exact orbit law for the coordinate permutation, transpose, and
  region permutations.
- Prove that the RM2Sub-S19 transfer ledger matches the frozen recurrence and
  covers every boundary case.
- Replace the modeled spectrum and nearest-binary64 ledger with authenticated
  inputs and outward-rounded verification.

### Asymptotic families

- For every proposed outer family, prove a one-block weighted bound \(M_B\)
  and choose a schedule satisfying \((N/B_N)M_{B_N}=o(1)\).
- For a candidate random-like structured outer, prove \(g_B=o(B)\), or prove
  the weaker transfer-weighted analogue, uniformly over the scalable family.
  A fixed \(B=256\) spectrum comparison has no asymptotic force by itself.
- For a structured outer at distance \(0.11\), prove a uniform spectrum or
  transfer-weighted comparison with the random rate-half outer. The selected
  Golay--BA family has a complete \(0.101\) certificate but its current
  half-Bernoulli spectrum cost is too large for the near-GV endpoint.
- Determine whether a second accumulator admits a comparably sharp closed
  contraction and materially improves the rate-half distance ceiling beyond
  \(0.0449101\).
- For a transpose followed by independent region shuffles, derive the joint
  column-profile law and multiply the conditional two-state accumulator
  kernels.  The pure-transpose obstruction does not settle this richer case.
- For the repeated EBCH128 plus two-accumulator BA family, prove or refute the
  matching upper bound \(M_{B,2}(\rho)=O(B^{-16})\).  The present lower bound
  and diagnostics do not establish it.
- Determine whether a sharper sparse prefactor can reduce the certified Random
  SPIN block constant below 17.  The optimized distance \(0.11002\) is within
  \(7.9\times10^{-6}\) of the rate-half Gilbert--Varshamov point; closing that
  final gap is also open.  The proof-model existence theorem itself is closed.
- The random-outer and selected Golay--BA certificates both keep
  \((t,s)=(128,19)\) constant. Other structured outer families still require
  their own uniform comparison.
- The one-sampled Golay--BA theorem is closed with block constant \(39/4\).
  The same construction and proof architecture give linear ordinary and
  transposed work. The fixed-occupation proof is now weight-coupled. Reducing
  the constant further requires a stronger norm witness or a weight-coupled
  four-state sparse transfer; it is not a gap in the \(c=39/4\) theorem.
- The finite Golay \(B=240\), \(k=2^{20}\) calculation with RM2Sub-S19 does
  not yet yield a 40-bit complete certificate. Replacing the inner by one
  shared random Toeplitz convolution gives a separate complete 180-bit
  certificate for a fixed repeated outer and route. The power-of-two
  EBCH32--ParityFanout \(B=256\) route now has a complete 10% mathematical
  certificate, including bounded setup failure, with 45.323934 bits of
  combined margin. Its exact setup test is exhaustive and impractical;
  efficient setup and the RM2Sub implementation-equivalence audit remain open.

## Conflicts and dependencies

- **Paper architecture:** the paper should use the type-orbit theorem for the
  structured section and reserve the weight-slice theorem for a uniform global
  interleaver.
- **Linear-time audit:** the companion certificate proves that the selected
  Golay--BA family keeps \((t,s)=(128,19)\) and has linear ordinary and
  transposed work at distance \(0.101\).
- **Implementation cleanup:** mathematical admissibility must not be inferred
  from frozen tiling constants.  Cleanup should preserve the frozen instance,
  while a future scalable implementation uses a separate manifest.
- **Naming:** current analytic source files use a sliding dense outer, whereas
  the canonical named Accumulator SPIN and Random SPIN variants use random
  block outers.

No merge conflict was created because this workstream changed no shared file.

## Earlier RM2Sub-S19 finite task

For the RM2Sub-S19 route, the independent-row B=240
transposed map beat the same-binary BCH baseline by 2.699%. A conditioning
window sweep replaced \([25,215]\) by \([23,217]\): the new outward receipts
give 46.479 bits for occupation one and 84.617 bits for occupations 2 through
64. The old window's available conditional-spectrum bound loses 598.911 bits
at \(Q=L\), while the revised window loses 81.029 bits and leaves a 152.637-bit
binary64 central-shell margin. The next smallest useful task is the
finite weight-coupled transfer for occupations 65 through 8,832. This task is
not needed by the separate Toeplitz certificate. In parallel,
setup needs an efficient exact or authenticated \(\mathcal G_{240}\) test.
`FINITE_K20_PROOF_OBLIGATIONS.md` gives the complete ordered ledger.

The dense work now follows the two-track plan in
`FINITE_K20_TWO_TRACK_DENSE_PLAN.md`. At 11%, a common-norm column-Hölder
relaxation fails by more than 1.7 million bits and is rejected. At 10.9%, all
116 sampled five-band compositions at \(Q=L\) close; the smallest sampled
margin is 6,016 bits. A complete barycentric composition cover and all
occupations below \(L\) remain open.

## Finite \(k=2^{20}\) EBCH32--ParityFanout--BA checkpoint

The same candidate now has a complete outward first-moment certificate at
relative distance 10%.  With \(N=2^{21}\) and \(d=209715\), the certificate
proves
\[
  \Pr[d_{\min}\le209715]<2^{-51}.
\]
The dense occupations \(65\) through \(8192\) collapse to a one-band
19-interval cover and have 468.225918 outward bits of margin.  Reusing the
stronger 11% receipts for occupations 1 through 64 gives a combined margin of
51.845632 bits.  `FINITE_K20_D10_CERTIFICATE.md` states the probability space,
convex interval argument, artifacts, and remaining implementation obligations.

A bounded setup wrapper now makes this an unconditional finite algorithm. It
tries at most six independent BA candidates per row and accepts a candidate
only after exact enumeration of all \(2^{128}\) row messages verifies the
weight window \([24,232]\). Exact rational combination of setup abort and bad
distance proves

\[
 \Pr[\text{setup aborts or }d_{\min}\le209715]<2^{-45},
\]

with 45.323934 display bits of margin. This closes the mathematical 40-bit
target while leaving practical setup open: the exhaustive test can require
\(6\cdot8192\cdot2^{128}\) row-message evaluations. The frozen implementation
baseline is under
`constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/`;
the exact EBCH32--BA row law still requires an implementation-equivalence
audit there.

The active power-of-two candidate uses \(B=256\), \(L=8192\),
\(N=2^{21}\), \(k=2^{20}\), and \(d=230686\). Each independently sampled row
uses eight genuine \([32,16,8]\) extended-BCH constituents, a disjoint
ParityFanout-31x33 map, and two independently interleaved terminated
accumulators. The exact setup law and spectrum identity, including the
deterministic all-ones transition \(256\mapsto223\), are in
`FINITE_K20_EBCH32_PARITYFANOUT_SETUP.md`.

New outward results are:

- setup acceptance probability at least 0.998816961305865 and expected row
  trials at most 1.0011844399324061;
- \(Q=1\) expected bad-word margin 51 certified integer bits (51.845632 bits
  for display); and
- aggregate \(Q=2,\ldots,64\) margin 95 certified integer bits (95.572310
  bits for display).

The full 11% theorem is still open. The dense \(Q=65,\ldots,8192\)
diagnostic cover is resumable after 8,060 processed tetrahedra: 3,957
accepted, 149 pending, zero rejected. An exact audit reconstructs all 4,106
leaves and verifies the explicit root-volume identity. A 256-bit Arb verifier
covers all 3,716 nondelegated nonempty accepted cells and gives 43.525570
outward bits for their aggregate; 234 accepted cells have empty lattice
boxes.

Seven geometric leaves contain only the external three-band type
\((3344,1,21)\), where the three-band relaxation fails. A local split of
weights \([24,100]\) into \([24,54]\) and \([55,100]\) yields a gap-free
17-interval cover of all 3,345 internal compositions. Its 256-bit Arb
aggregate margin is 198.029788 bits. The main partition delegates those seven
leaves to this local certificate. The partial dense sums still omit 149
pending cells and do not form a full dense certificate.

Before merging a theorem claim, complete and outward-verify the dense cover,
add the delegated split-low bound exactly once,
audit the RM2Sub and region-permutation interfaces, provide an efficient
\(G_{256}\) setup test, and benchmark an implementation of this exact outer.
The existing frozen ParityFanout implementation benchmark does not measure
the EBCH32--BA pipeline.

## Finite repeated-EBCH128 11% closure audit

The fixed repeated \([128,64,22]\) extended-BCH candidate does not yet have
a complete finite 11% certificate. The prior 26.714-bit number covers only
occupation one and uses nearest-binary64 arithmetic. An exact uniform-shell
diagnostic covers occupations 2 through 100 with 28.450 bits. Their combined
diagnostic margin is 26.335 bits.

A parity-aware body comparison removes the old one-bit-per-active-row charge:
it represents a uniform even row by 127 fair coordinates and omits the final
parity-determined region. This closes sampled occupations through \(Q=16100\),
but it fails by about 14,545 bits at \(Q=L=16384\). The failure is structural.
At this length, an ideal random rate-one-half code has only about 178 bits of
first-moment margin at 11%, so discarding one complete region is unaffordable.

The exact one-coordinate marginal does not supply a shortcut. A common-norm
column-Hölder calculation at \(Q=L\) fails by about 1.77 million bits and is
rejected. Closing the candidate requires a 128-region transfer that preserves
the exact BCH row dependence, or a proved RM2Sub parity-region decoupling
lemma. It must also cover mixed configurations containing the unique all-one
BCH word and must be verified with outward arithmetic.

`FINITE_K20_REPEATED_EBCH128_D11_STATUS.md` states the exact probability
space, receipts, correction, and remaining lemma. The integration owner
should not describe this candidate as certified at 11%. The next smallest
useful task is the parity-region decoupling lemma; it is worth resuming only
if the repeated-EBCH128 implementation remains competitive.

### RandomStepConv comparison for the same EBCH128 outer

Replacing RM2Sub by RandomStepConv-M30 produces a complete outward
certificate. At

\[
 (L,L_0,N,D)=(16560,16384,2119680,233165),
\]

the verifier proves \(d_{\min}\ge233165\) except with probability below
\(2^{-26.1921836158}\). Thus the smaller fixed constituent closes every
occupation without any random-outer spectrum test.

The initial all-occupation even-row majorant fails by 11,624 bits at
\(Q=3787\). Its Bernoulli coefficient bound overweights completely inactive
coordinate regions. It also aligns every row's parity pivot into one final
region and deletes that region's inputs. At the dense endpoint this produces
an artificial 16,560-position silent interval. Raising the random state from
30 to 45 bits does not repair the intermediate-occupation loss. A separate
support-coverage bound was also rejected because keeping only one difference
per covered region creates artificial post-collision gaps.

These are proof-relaxation failures, not bad-code witnesses. The closing
proof assigns the parity pivot independently in each uniform-even reference
row. It then uses a positive multinomial-pivot transfer with two coefficient
fugacities. Occupation one uses the exact BCH spectrum, occupations 2 through
99 use an exact sparse recurrence, and the dispersed transfer covers the
remaining occupations. See
`FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_CERTIFICATE.md` and
`FINITE_K20_REPEATED_EBCH128_RANDOMSTEP_CONV_STATUS.md`.

`FINITE_K20_EBCH128_OUTER_DEFINITION.md` now fixes the constituent generator
polynomial and coordinate convention, the repeated-row map, and both setup
permutation layers. It also records the remaining provenance
obligation precisely: the verifier consumes the authenticated imported BCH
spectrum but does not locally derive that spectrum from the concrete
generator matrix.

The current candidate subsequently removed all 176 zero rows. It has exact
power-of-two dimensions \((k,N)=(2^{20},2^{21})\). The old 26.19-bit receipt
remains a historical padded certificate and does not transfer automatically.
At 11%, the pointwise spectrum loss exceeds the approximately 188.38-bit
random-code budget. The exact 11% target therefore remains open for this
unpadded instance.

Relaxing the target to 10.9% yields a complete outward certificate with
RandomStepConv-M20. For

\[
 D=\lceil0.109N\rceil=228590,
\]

the verifier proves

\[
 \Pr[d_{\min}<D]<2^{-7.58436084118}<2^{-5}.
\]

Occupation one uses the exact BCH spectrum. The parity-pivot transfer covers
occupations 2 through 16319. An exact-complement recurrence covers the final
65 occupations and removes the old dense-bound artifact. See
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_CERTIFICATE.md` and
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_STATUS.md`. The manifest
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_MANIFEST.json` binds the
theorem, spectrum, witnesses, verifier, and outward receipt.
`BITTRANSPOSE_RANDOMSTEP_CONV_PROOF_TEMPLATE.md` freezes the reusable
construction interface, occupation decomposition, certificate condition, and
substitution obligations.

`SINGLE_RANDOM_CONSTITUENT_SIZE_PROXY.md` replaces the Random SPIN paper's
independently resampled block injections by one random rate-half constituent
reused in every outer row. Its exact occupation-one diagnostic at 10.9% and
RandomStepConv-M30 gives 18.7859 bits for B=128, 55.3387 bits for B=256, and
127.9284 bits for B=512. B=256 is the first plausible size for the proposed
two-BCH-block permutation--accumulation experiment. At B=256, memory 21 is
the first tested value above 40 occupation-one bits. Memory 22 gives 46.5928
bits and is the safer first full-proof target.

`SINGLE_RANDOM_CONSTITUENT_Q2.md` closes the first shared-code correlation
case in nearest-binary64 arithmetic. At B=256, M22, and 10.9% distance, the
occupation-two margin is 64.1559 bits. Equal local messages give the
rank-one bottleneck. Distinct local messages give 99.7046 bits. The combined
occupation-one and occupation-two margin is 46.5928 bits. Occupations three
and above, and outward rounding, remain open.

`SINGLE_RANDOM_CONSTITUENT_ALLQ_RAMP.md` records the first attempt to cover
all occupations for the single random [256,128] constituent. A uniform
factor-8 spectrum event has at least 76.55% one-sample probability, but its
distance wrapper fails at Q=L by 18,102 bits because the factor is raised to
the 8,192nd power. The ideal factor-1 endpoint passes by 6,474 bits. A new
three-band event has at least 94.07% one-sample probability, and its
factor-1.5 central class passes every occupation. A categorical cover of
mixed fringe/central compositions remains open. The corresponding receipts
are `single_random_constituent_B256_allq_factor8_s22_d109.json` and
`single_random_constituent_B256_three_band_event.json`.

At M19, the complete nearest-binary64 analysis has 0.759 bits of margin at
10.9%. This is not an outward certificate. A lower distance target could
increase its sparse margin, but M20 is the current certified memory.

## Shortened BCH near length 256

`FINITE_BCH256_SHORTENING_OPTIONS.md` selects the genuine shortened BCH code
([250,125,\ge38]) as the clean non-power-of-two candidate. Exact GF(2) rank
audits prove the parameters. The code is even, omits the all-one word, and
has nonzero weights only in ([38,218]). At the nearest admissible full-rate
point,

\[
 (B,L,k,N,d)=(250,8448,1{,}056{,}000,2{,}112{,}000,232{,}320).
\]

Padding 7,424 zero input bits gives an exact (k=2^{20}) wrapper. The
three-coordinate shortening ([253,128,\ge38]) is rejected: its rate exceeds
one half, and even an ideal random first-moment benchmark misses 11% by about
12,113 bits.

The exact constant-weight packing envelope closes the nearest-binary64
occupation-one calculation with 13.777 bits and occupations 2 through 100
with 19.044 aggregate bits. These calculations do not assume an exact BCH
spectrum. They are not outward certificates.

The dense range remains open. A pointwise half-Bernoulli shell comparison
fails by 446,724 bits at (Q=L). Exact total row mass would recover about
414,224 bits, but the current transfer would still miss by about 32,500 bits.
The next proof target is therefore a band-coupled coefficient transfer using
the exact total mass, shell packing bounds, and weight range. Do not describe
the shortened candidate as certified at 11%.

## One-sampled BCH250 parity-fanout ramp

`ONE_SAMPLED_BCH250_PARITYFANOUT_RAMP.md` studies one fixed composition of
ParityFanout maps applied to the shortened ([250,125,\ge38]) code. The
composition is sampled once and reused in every row.

The exact-mass shell optimizer shows that one fanout layer is insufficient to
erase the unknown source spectrum. With 64 independent ParityFanout-31x33
layers, the expected central envelope is nearly random-like, but the proposed
three-band transfer fails. Broad tails fail at their vertices; narrower tails
fail on thin mixed faces. Two low-tail rows and 8,446 central rows have a
diagnostic margin of -228.861 bits. Splitting the tail into singleton weights
does not repair the categorical-conditioning loss.

The revised diagnostic uses 128 layers. The expected number of transformed
codewords with weight at most 16 or at least 234 is at most approximately
(2^{-41.40363381}). The pure complementary-bulk all-active calculation
retains 177.390 bits. These figures do not form a certificate: the remaining
near-tail rows need a defect-sensitive transfer, and the expected bulk
spectrum needs a simultaneous concentration theorem for one sampled and
reused code.

The exact random-model lesson is quantitative. A uniform pointwise spectrum
slack of (\varepsilon) bits costs (8448\varepsilon) bits in the dense
calculation. Retaining 40 bits requires (\varepsilon<0.01626) before other
finite costs. The next proof target is therefore an endpoint setup-failure
bound, exact treatment of defect weights 17 through 36 and 214 through 233,
and at most 0.01 bits per row of simultaneous bulk slack on weights 37
through 213. A sampled BA outer may use the same fixed-code interface if its
spectrum event can be proved directly.

## BCH250-124 with independent row-local fanout

`FINITE_K20_BCH250_124_ROWLOCAL_FANOUT_STATUS.md` records a new finite 11%
route. It repeats one fixed $[250,124,\ge 38]$ subcode of the audited BCH250
constituent. Every row receives an independent 56-layer
ParityFanout-31x33 inner wrapper.

The exact parent dimensions are

\[
 (K_0,N,D)=(1{,}063{,}424,2{,}144{,}000,235{,}840).
\]

Zero-shortening 14,848 parent input coordinates gives exactly $2^{20}$
message bits. A shell-sensitive proof separates weights 1 through 12 and
238 through 250. The probability that any row wrapper contains such a word
is below $2^{-45}$. On weights 13 through 237, the expected pointwise density
excess is below $389/1250$ bits.

The outward occupation receipts give strict bounds of $2^{-53}$ for $Q=1$,
$2^{-52}$ for $2\le Q\le31$, and $2^{-41}$ for $32\le Q\le8576$ on the
central-shell event. Together with the tail event, their sum is below
$2^{-40}$. The displayed endpoints combine to 41.507281 bits.

The theorem and hash-bound manifest are in
`FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_CERTIFICATE.md` and
`FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_MANIFEST.json`. The earlier
256-layer certificate remains valid but is superseded for implementation.

The packed 56-layer action is now implemented. A pinned i7-13700H benchmark
measures 2.177353 ms for the wrapper and 0.006059 ms for the post-BCH
zero-layer copy. The isolated increment is 2.171294 ms for all 480,256
row-layer applications. The explicit schedule occupies 29.3125 MiB.
`FINITE_K20_BCH250_FANOUT56_IMPLEMENTATION_BENCHMARK.md` states the kernel,
checks, probability-space boundary, and receipt.

The replacement performance proxy preserves the frozen 256-by-8192 packed,
tiled, fused transpose encoder. A pinned Ryzen 9 7950X benchmark measures
10.650813 ms with zero layers and 20.454284 ms with 56 layers. Fanout-56 adds
9.803471 ms by the difference of medians and multiplies total time by
1.920443. Both modes pass staged correctness checks.

The replacement is a BCH256 performance proxy, not the exact BCH250 encoder.
The proxy uses 8,192 rows instead of 8,576. Scaling only the measured fanout
increment by the row-count ratio gives 10.263009 ms, but that value is an
engineering estimate. See `FINITE_K20_BCH256_FANOUT56_PROXY_BENCHMARK.md` and
`FINITE_K20_BCH256_FANOUT56_PROXY_MANIFEST.json`.

A proof-compatible setup sampler, ordinary forward encoder, exact BCH250
implementation, and complete implementation-equivalence proof remain open.
The unfused proof-model count is 3,528 scalar XORs per row, or 30,256,128
over all exact rows. A 52-layer candidate fails the current tail-plus-central
pointwise argument; a multi-band proof may do better.

## Repeated BA with random Toeplitz convolution

`TOEPLITZ_PREFIX_ARGUMENT.md` records the reusable theorem, probability
space, triangular suffix proof, all-message prefix identity, and its exact
scope.

FINITE_K20_RANDOM_CONV_PREFIX_ANALYSIS.md gives the complete all-message
analysis for random lower-triangular Toeplitz convolution. The convolution
has unit diagonal. For each fixed nonzero input difference, its activation
output is one and every later output is independent and fair.

For a fixed routed outer code, the analysis records the dimension of the
subcode whose first \(t\) coordinates are zero. A summation-by-parts identity
then evaluates the expected number of bad final words over all nonzero parent
messages. This prefix-rank calculation uses Gaussian elimination and does
not enumerate messages.

At \(B=720\), independent uniform random outer rows give 189.230985
binary64 bits. Independent BA rows give 20.244610 bits. Their occupations
\(Q\ge2\) contribute only \(2^{-41.491}\), so rare one-row BA words explain
essentially the entire loss.

Four fixed setups using one repeated BA code give 190.521692 bits each. Four
fixed setups at the current \(B=240\) design point give the same margin.
Thus increasing \(B\) suppresses rare bad BA samples but does not improve the
tested typical setups.

The random-outer and BA spectra differ by hundreds of bits at weights 2
through 5, despite their similarity near weight \(B/2\). At the critical
prefix, 78% of the routed input is zero and 22% remains. The random suffix
then has expected weight 11% of \(N\). Weights 2 through 5 supply 98.806% of
the BA survivor mass at that prefix.

The repeated-code ensemble requires high moments of local prefix counts.
Powers of the expected BA spectrum are invalid. Conditioning on one fixed
setup and checking its actual prefix ranks avoids those high moments.

The proof task is now closed for one explicit \(B=240\) setup.
`FINITE_K20_BA240_REPEATED_TOEPLITZ_CERTIFICATE.md` freezes the BA code and
all routes with a stable SplitMix64/Fisher--Yates specification. The
outward verifier recomputes every prefix rank and obtains

\[
 \log_2 \mathbb E[Z_{\mathrm{bad}}]
 <-190.5216915553.
\]

The certificate claims the conservative integer bound
\(\Pr[d_{\min}\le233164]<2^{-180}\), hence minimum relative distance above
11% for the \(k=2^{20}\) shortened code. This is a finite proof-model result,
not a linear-time claim. The next implementation task is to benchmark forward
and transposed binary Toeplitz multiplication at \(N=2{,}119{,}680\), or to
identify a bounded-state inner that preserves enough of the suffix law.
`FINITE_K20_BA240_REPEATED_TOEPLITZ_MANIFEST.json` hash-binds the proof,
verifier, and outward receipt.

## Fixed outer: bridge back to RM2Sub-S19

`FINITE_K20_RM2SUB_FIXED_OUTER_BRIDGE.md` records the first attempt to port
the prefix proof. The port is not direct. RM2Sub observes only 19 syndrome
bits per 128-bit epoch. Its complete syndrome map has rank at most 314,640.
The \(2^{20}\)-dimensional shortened outer therefore contains a silent
subcode of dimension at least 733,936. RM2Sub acts as the identity on this
subcode for every multiplier schedule.

This is not evidence of a low-weight word. It identifies the exact missing
property. The finite proof must certify the weight enumerator of the silent
subcode and weight-sensitive profiles of every suffix-silent subcode. Raw
coordinate-prefix ranks record neither quantity.

A random code at the forced silent dimension has a 326,094-bit binary64
benchmark at the 11% cutoff. Charging a \(B^2\) selected-BA factor in all
8,832 rows costs 139,667 bits and leaves a provisional 186,427 bits. This is
a budget diagnostic, not a comparison theorem for the structured syndrome
map.

The live state remains strong. A Walsh-positivity argument bounds every
coset moment by the exact selected-\(A\) weight enumerator. An optimistic
uninterrupted-live diagnostic gives 228,109 bits of single-word Chernoff
margin over 16,559 epochs and crosses 40 bits at 3,684 epochs. This diagnostic
does not include termination or deterministic silent segments.

The next proof task is a weighted syndrome-prefix evaluator for the fixed
outer. It should attach a fugacity to silent input weight and then join the
result to the existing zero/deterministic/uniform-live/punctured-live RM2Sub
transfer. `audit_rm2sub_fixed_outer_bridge.py` and
`rm2sub_fixed_outer_bridge_audit.json` authenticate the structural audit.

## Random bounded-memory convolution baseline

`FINITE_K20_RANDOMSTEP_CONV_BASELINE.md` defines the clean finite-state
baseline. For memory \(M\), every bit position receives an independent
uniform linear map

\[
\mathbb F_2^{M+1}\longrightarrow\mathbb F_2^{M+1}.
\]

The input is one bit, the output is one bit, and the other \(M\) coordinates
are the current and next state. The matrices are sampled once and shared by
all messages.

For a full random linear outer, an exact two-state first-moment calculation
covers all nonzero messages. At \(k=2^{20}\), \(N=2{,}119{,}680\), and the
strict 11% cutoff, nearest-binary64 margins are 1,422.7138 bits at \(M=10\)
and 11,421.6932 bits at \(M=19\). Memory 9 fails. The unshortened
\(N/2\)-dimensional parent retains 157.6932 bits at \(M=19\).

For the desired one-repeated Golay--BA-3 outer, the exact occupation-one
calculation closes with 21.5360 bits at \(M=19\) and 41.6332 bits at \(M=23\),
conditional on the sampled constituent having no nonzero word outside
weights 23 through 217. That good event has probability at least
0.9936608994938328.

These are not yet outward certificates. The repeated-BA result still needs
occupations \(Q\geq2\). The result nevertheless isolates the obstruction:
bounded memory is adequate, while RM2Sub's fixed low-rank syndrome interface
creates the silent-subcode proof problem. The unrestricted random-step model
is a comparator rather than a performance proposal; explicit step maps cost
about 101.1 MiB at \(M=19\).

The follow-up repeated-outer calculation is now a complete outward finite
certificate. It samples one random \([512,256]\) constituent, retains it when

\[
 A_w\le(13/2)\mathbb E[A_w]
 \qquad(1\le w\le512),
\]

and repeats that same constituent in every outer row. There are 4,140 parent
rows, of which 44 complete information rows are zero-shortened. Thus

\[
 (k,N,D,M)=(2^{20},2{,}119{,}680,233{,}165,30).
\]

The exact Markov--Cantelli--rank calculation proves that one sampled
constituent is admissible with probability greater than
0.6365388449160215. Conditional on any admissible constituent, the positive
outward transfer calculation covers all 4,096 nonzero-row occupations and
gives

\[
 \Pr[d_{\min}<233{,}165]<2^{-108}.
\]

Trying at most 28 independent outer candidates and retaining the first
admissible one gives a combined setup-abort-or-distance-failure bound below
\(2^{-40}\). This is a literal 11% result because
\(233{,}165=\lceil0.11N\rceil\). The exact spectrum test may enumerate
\(2^{256}\) messages, and the explicit \(31\times31\) step matrices occupy
about 242.83 MiB. The theorem therefore closes the mathematical random-inner
comparator but makes no competitive implementation claim.

`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md` states the
construction, probability space, theorem, and proof. The associated search,
outward verifier, and receipts are
`evaluate_repeated_random512_randomstepconv_g1.py`,
`repeated_random512_randomstepconv_g1_s30_allq_d11.json`,
`certify_repeated_random512_randomstepconv_g1_outward.py`, and
`repeated_random512_randomstepconv_g1_s30_allq_outward_d11.json`. The file
`FINITE_K20_REPEATED_RANDOM512_RANDOMSTEP_CONV_MANIFEST.json` hash-binds this
set.

The later one-shot variant removes the acceptance test from setup and is now
an outward certificate. It samples one uniform ([512,256]) generator and
repeats it in all 4,096 rows. The explicit integer spectrum event fails with
probability at most (2^{-42.2614729204}). Conditional on that event, the
all-occupation distance failure is at most (2^{-51.6439589890}).

The combined theorem proves

\[
 \Pr[d_{\min}<228{,}590]<2^{-42.2593129905}
\]

at ((k,N,M)=(2^{20},2^{21},22)). Thus the minimum relative distance is at
least 10.9% except with more than 40 bits of margin. The sampler draws
131,072 unbiased bits once; it does not enumerate codewords or reject a
sample.

The dense matrix can now be replaced by a compact pairwise-systematic
constituent. Over \(\mathbb F_{2^{256}}\), sample \(a,b\) and encode
\(u\mapsto(u,au+bu^2)\). For distinct nonzero messages \(u,v\), the
determinant \(uv(u+v)\) proves that their parity halves are independent and
uniform. Thus the shell counts retain the exact variance bound used by the
uniform-matrix certificate. The systematic half guarantees injectivity.

The original cap vector remains valid, so the sparse and dense transfers are
unchanged. The outward combined margin is 42.2593129905 bits. The outer
sampler uses 512 unbiased bits and the outer encoder uses two field
multiplications per row. The result is in
`FINITE_K20_PAIRWISE_SYSTEMATIC_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md`.
It strengthens the random comparator but does not prove the BA concentration
lemma.

The closing lemma observes that the RandomStepConv tilted moment decreases
with the Bernoulli input density. It therefore merges the symmetric spectrum
tails into one defect category. An exact shared-category recurrence covers
(Q=3,ldots,159), while a convex simplex cover certifies every dense integer
composition. The integration owner should use
`FINITE_K20_ONE_SHOT_RANDOM512_RANDOMSTEP_CONV_CERTIFICATE.md` as the current
random-inner finite baseline. `SINGLE_RANDOM_CONSTITUENT_ONE_SHOT_STATUS.md`
retains the failed categorical relaxation as proof-search history.

This certificate makes no RM2Sub or implementation-performance claim. The
next smallest useful task is to replace RandomStepConv by a structured inner
whose transfer is dominated by the certified two-state law, without changing
the one-shot outer interface.

A second algebraic bridge preserves the fixed EBCH direct sum. Identify its
512-bit output with \(x\in\mathbb F_{2^{512}}\), sample \(a,b\), and apply
\(x\mapsto ax+bx^2\). The same determinant argument gives pairwise-uniform
images. The restriction loses rank with probability below \(2^{-256}\).
The existing cap-conditional transfer therefore proves 42.2593129905 bits
overall with zero accumulator stages. This is a certified fallback with two
\(\mathbb F_{2^{512}}\) multiplications per row; it is not a BA-only result.
See `FINITE_K20_EBCH128X4_LINEARIZED_MIXER_RANDOMSTEP_CONV_CERTIFICATE.md`.

## Structured repeated-constituent continuation

The current structured continuation is documented in
`EBCH128X4_BA_REPEATED_CONCENTRATION_STATUS.md`. It uses the direct sum of
four fixed extended BCH \([128,64,22]\) codes followed by three independently
interleaved accumulator stages. One sampled \([512,256]\) constituent is
repeated in all 4,096 outer rows.

This route is not yet certified. Its exact missing concentration target is

\[
 \operatorname{Var}(A_w)\le512\,\mathbb E[A_w]
 \quad\text{for every nonzero shell }w.
\]

Conditional on this target, the spectrum-event failure is at most
\(2^{-42.255553}\). The RandomStepConv-M22 binary64 transfer covers all
occupations below 160; its worst aggregate slice is \(Q=3,\ldots,16\) with
46.114992 bits. Sampled dense compositions have at least 559.784886 bits.
The resulting projected combined margin is about 42.1594 bits if the complete
dense cover retains at least 46.11499 bits.

The variance target cannot be derived from the ordinary BCH spectrum alone.
It requires a joint pair-type enumerator or an association-scheme operator
bound. A complete dense convex cover and outward rounding are also still
required. Integration must not present this route as a theorem until those
obligations close.

The one-stage pair-type transition is no longer open. For accumulated pair
states \(s,s'\in\mathbb F_2^2\), the transfer monomial is
\(x_{s+s'}y_{s'}\). Coefficient extraction from the four-by-four matrix power
gives the exact common-permutation accumulator kernel. The independent
dynamic-program and exhaustive-permutation implementations agree in
`accumulator_pair_type_kernel_small_exact.json`. The open step is a tractable
three-stage contraction bound strong enough to imply the stated shell
variance inequality.

The exact MacWilliams transform gives base dual distance 22. Hence two
independent uniform base-code words have 21-wise independent pair symbols,
so all pair-type polynomial modes through degree 21 match the multinomial
stationary law. The rank-two conditioning correction has mass only
\(3\cdot2^{-256}-2\cdot2^{-512}\). The intended contraction proof should act
only on degree-22-and-higher modes. This constituent-specific cancellation is
the main reason the variance target may be substantially easier than a
worst-start mixing theorem.

A one-word singular-value proxy limits that optimism. After removing
polynomial degrees zero through 21, its three-stage generic L2 norm contracts
by only 15.2928 bits. Therefore dual distance alone is insufficient. The open
pair proof must also use the target shell or the specific EBCH-direct-sum pair
distribution. Tested tail splits and shoulder categories did not improve the
all-regime transfer; the original two-band \(F=512\) interface remains the
current target.

Small exact pair chains give a concrete representation clue. At lengths 8
and 16, the leading nontrivial pair-chain singular values match the one-word
values, each with multiplicity three. The length-16 direct-sum model has
exact worst shell variance-to-mean ratios 1.4622, 1.1628, 1.0415, and 1.0034
after one through four stages. This is diagnostic evidence, not a
length-512 theorem.

The most promising remaining proof interface is a triangle-convolution
factorization of the exact pair kernel. For independent messages \(U,V\), the
three messages \(U,V,U+V\) are pairwise independent. A pointwise majorant
\(\Phi(h_1,h_2,h_3)\le f(h_1)f(h_2)g(h_3)\) reduces the second moment to
the exact diagonal mean plus
\(M^2\mathbb E[f(H)^2]\mathbb E[g(H)]\), which uses only the ordinary base
spectrum. The length-8 exact model gives worst bounds 3.0755, 1.0330, and
0.9416 after one, two, and three stages. A length-16 direct-sum model gives
19.2869, 2.7092, and 1.2642, so two and three stages beat its
\(F=B=16\) analogue target. The open task is a non-materialized, outward
length-512 factorization for every output shell.

The exact-kernel probe now reaches length 24. For direct sums of one, two,
and three RM(1,3) constituents, the three-stage triangle-convolution ratios
are 0.9416, 1.2642, and 1.5681. At length 24, the two-stage bound is 7.7268 and
the exact ratio is 2.1517. An inverse finite sweep shows that BA-2 needs a
uniform variance factor below 2.8285 merely to pass the outer-event-plus-
\(Q=1\) gate, with almost no remaining budget. BA-2 therefore appears too
fragile for the present proof interface; BA-3 remains the primary target.

A reduction to the three binary character projections of the pair kernel was
tested and rejected. Exact rational checks at lengths 8 and 16 verify its
one-stage inequalities, including equality cases. Iterating the resulting
factored envelope gives a length-8 three-stage variance bound above 8,232;
jointly optimizing all intermediate factors only lowers it to 8,086. A
successful length-512 proof must retain the genuine four-symbol pair kernel or
an equivalent representation-theoretic compression.

There is now an exact candidate compression boundary. The span of functions
of the three character weights \(\operatorname{wt}(x)\),
\(\operatorname{wt}(y)\), and \(\operatorname{wt}(x+y)\) is reducing for the
pair operator; its restriction is determined by the one-word accumulator
chain. The complementary interaction block has binary64 norm 0.279680 at
length 8 and 0.185857 at length 16. This identifies precisely what the
length-512 verifier must retain. The norm alone is not a certificate: a
generic L2 bound still pays for the initial density of the fixed rate-half
code, so the shell-specific pointwise factorization remains necessary.

Genus-two MacWilliams positivity supplies another exact interface. The
unknown pair enumerator and its dual are nonnegative tables related by the
four-symbol MacWilliams transform. Their three character marginals are the
known ordinary primal and dual spectra. The BA shell second moment is linear
in the primal table, so any outward dual LP solution proves the required
variance bound without enumerating message pairs.

At length 16, marginal constraints alone allow variance inflation 22.9738.
Adding only nonnegativity of the transformed table reduces all tested one-,
two-, and three-stage optima to the exact values within binary64; the worst is
1.4622. The result identifies a promising certificate cone, but not yet a
length-512 certificate. The existing character-projection envelope is too
loose even on the actual pair table, reaching 115779.95 after three stages.
The next implementation must therefore expose the exact second-order
input-output distribution, or a substantially tighter shell-specific bound,
in a compressed form.

An optimal Christoffel-moment calculation also closes the primal/dual-distance
question negatively. At twelve accumulator stages, a 43.7852-bit event gives
primal support (44,ldots,468) and dual distance 38, but the resulting M64
(Q=1) bound is (2^{65.34}). At thirteen stages, dual distance 44 only
improves that bound to (2^{53.25}). Thus no argument based solely on primal
and dual distance can provide the required spectrum envelope, even with many
more stages and a substantially larger inner state.

An optimized three-marginal intersection inequality was also tested. It uses
the fact that two weight-(w) outputs have difference weight at most
(2\min(w,512-w)), but it retains no pair-overlap type. The resulting exact
Triangle--Holder reduction gives a worst variance factor (2^{256}) through
sixteen stages. This route is therefore no substitute for the four-symbol
pair kernel.

A two-stage transfer optimization does not close the 40-bit target under the
current \(F=512\) cap interface. Forcing all shells through weight 26 to zero
costs 40.5002 bits in the unbiased outer proof event, while the resulting
RandomStepConv-M40 \(Q=1\) term has only 38.1036 bits; their union has at most
37.8528 bits. Lower cutoffs are worse. This rejects the present two-stage
proof interface, not the construction.

### Conditional universal BA-128 closure

There is now a complete all-occupation fallback conditional on one finite
accumulator lemma. Let (C_0) be the direct sum of four fixed extended BCH
([128,64,22]) codes. Apply 128 independently permuted prefix accumulators
to (C_0), sampling the permutations once and repeating the resulting
single constituent in all 4096 rows.

Assume that one stage has second stationary-(L_2) singular value at most
(63/1000) on both nonzero words and ordered pairs of distinct nonzero
words. Stationary chi-square contraction then bounds every shell-count mean
and second factorial moment, uniformly over the fixed starting code. An
outward Markov--Cantelli calculation proves that BA-128 obeys the exact
random-([512,256]) cap vector except with probability

\[
 2^{-40.3474252482}.
\]

BA-127 has only 36.7702 bits at this gate, so 128 is the first passing count
under this universal bound. The existing outward cap-conditional sparse and
dense receipts already cover every (1\le Q\le4096) and contribute at most
(2^{-51.6439589890}). Consequently, under the singular-value lemma,

\[
 \Pr[d_{\min}<228{,}590]<2^{-40.3468518020}.
\]

This is a conditional outward finite certificate for relative distance at
least 10.9% with RandomStepConv-M22. It is not yet an unconditional theorem,
and 128 accumulator stages are not an implementation proposal.

The remaining lemma has a concrete decomposition target. The three scalar
character sectors are reducing. The length-512 one-word second singular
value is 0.0626241431 in binary64. On the interaction complement,
exact-kernel/binary64 probes at even lengths 4 through 20 lie below (3/B).
Proving the one-word bound, the interaction bound, and their complete
scalar-sector combination would remove the final condition. No further
finite-(Q) transfer proof is required for this fallback.

An outward 256-bit Arb calculation now proves the weaker one-word inequality
\(\sigma_2<0.11\) by bounding the squared centered Hilbert--Schmidt norm by
\(0.011892115495235<(0.11)^2\). If the complete ordered-pair operator obeyed
the same \(0.11\) bound, the cap event would first pass at BA-161 with 41.5690
bits. This does not close BA-128: its \(0.063\) one-word target and the pair
interaction bound remain open.

The supporting scripts and receipts are
`evaluate_ba_pair_singular_mixing_sweep.py`,
`ba_pair_singular_mixing_sweep_B512.json`,
`evaluate_ba_pair_singular_random_cap_event.py`,
`ba_pair_singular_random_cap_event_B512.json`,
`certify_ba128_pair_singular_random_cap_event_outward.py`,
`ba128_pair_singular_random_cap_event_outward_B512.json`,
`certify_ebch128x4_ba128_randomstepconv_combined.py`, and
`ebch128x4_ba128_randomstepconv_combined_outward_s22.json`. The weaker
one-word fallback is certified by
`certify_accumulator_oneword_singular_outward.py` and
`accumulator_oneword_singular_B512_outward.json`; the conditional BA-161
diagnostic is `ba_pair_singular_random_cap_event_B512_lambda011.json`.

### Pure sparse-mixer--accumulator concentration route

A separate approved route replaces the BCH--BA constituent by sparse linear
maps and accumulators. Each sparse-map output independently samples a fixed-
degree input neighborhood. One sampled constituent is reused in every SPIN
outer row.

The one-word and ordered-pair laws are exact Krawtchouk formulas. Four small
instances match exhaustive sparse-map enumeration in rational arithmetic.
At \((k,n)=(256,512)\), the degree-33 one-stage candidate passes the
kernel and extreme-shell first-moment gate with a 44.6379-bit diagnostic
margin. The lower-XOR candidates \((17,3)\) and \((15,2,3)\) also pass that
gate diagnostically. None of these first moments proves a realized-spectrum
event.

The one-stage variance now has an exact dual formulation. Its shell count is
a weighted sum of indicators that specified sparse-map rows XOR to zero.
The covariance kernel depends only on two row-subset weights and their
intersection, so it has an explicit symmetric-group block decomposition.
Two small instances reproduce exhaustive rational shell variance exactly,
and their block eigenvalues reproduce the full covariance spectra.

The three-index character moments also factor through a radial
\((k+1)\)-state walk:

\[
 \mathbb E[\beta(X)^a\beta(Y)^b\beta(X+Y)^d]
 =\sum_{h=0}^{k}\binom kh v_h(a)v_h(b)v_h(d).
\]

At \(k=256\), this reduces the message-pair calculation to a
\(257\times513\) table. An exact small verifier checks 165 instances of the
identity.

Two simpler relaxations are rejected. A global covariance norm loses to
two-row collisions, whose probability is
\(\binom{256}{33}^{-1}\approx2^{-138.18}\). A weight-level block norm gives
4611.74 where the exact moderate-size variance factor is 2.7906, and
15854.46 where the exact factor is 4.6473. Both relaxations discard essential
signed or sector-consistent cancellation.

The remaining target is to evaluate the irreducible sector quadratic forms,
or their sector-projection norm bound, with directed rounding at length 512.
Only then can the result be compared with the sufficient
\(\operatorname{Var}(N_w)\leq512\mathbb E[N_w]\) cap interface. No
high-probability realized-spectrum certificate for the sparse-mixer
constituent is claimed yet.

The current route and its open obligations are recorded in
pure_expander_accumulate/ONE_STAGE_VARIANCE_ROUTE.md. Supporting scripts and
receipts in that directory include the pair-kernel checks, finite first-
moment sweeps, moderate-size variance probes, dual-walk energy calculation,
block verifier, radial-factorization verifier, and the rejected weight-level
block probes.

The irreducible-sector projector energies also have an exact finite formula
through the Johnson association scheme. Small rational checks reconstruct
each row-subset level energy and the total Parseval energy. The remaining
accumulator-side input is the table of signed coefficient correlations for
two row subsets of fixed weight and intersection. This table is a memory-one
two-sequence transfer object, but it has not yet been evaluated outward at
length 512.

A Walsh conjugation now provides a more direct low-shell route. Writing
\(f_w(z)=\mathbf 1\{\operatorname{wt}(Az)=w\}\), the dual vector satisfies
\(c_w=Wf_w\). The normalized conjugate of the covariance operator has the
same symmetric-group sectors, while \(f_w\) has an exact run-count level
profile supported only through level \(2w\) for \(w\le n/2\). A valid
level-norm bound in this basis gives 1.07079 at the low-shell
\((32,64,7,5)\) analogue, against the exact ratio 0.999994, and 10.4285 at
the matched \((64,128,9,10)\) analogue. These values use nondirected
binary64 arithmetic and are not a certificate.

This route is relevant because the two-band finite transfer can bound the
central aggregate by the exact constituent size. Its first shellwise target
is weights 42 through 79 and their complements; weight 42 requires conjugated
primal levels only through 84. The next obligation is a compressed,
outward-rounded length-512 evaluation of those block entries. The exact
Johnson correlation table remains the fallback if the conjugated level bound
does not close.

The conjugated block has an exact direct form. For a fixed ordered message
pair with one-row output matrix \(P\), sector \(j\) of the product kernel is
\(\det(P)^j\operatorname{Sym}^{n-2j}(P)\). A small direct aggregation over
all message pairs agrees with the Walsh-conjugated dual blocks to relative
binary64 residual below \(1.4\times10^{-15}\). The open scaling step is now
precise: aggregate these symmetric-power entries over the target message-pair
types while retaining complement cancellation and outward rounding.

That aggregation is now implemented. An exact insertion inequality for
diagonal symmetric-power entries bounds every sector at fixed primal level by
the maximum of sectors zero, one, and two. Positive semidefiniteness converts
those diagonals into a bound on all cross-level entries. The target program
streams all 2,862,209 message-pair types and uses one-sided local sums, so the
mathematical upper bound needs no global cancellation.

At \((k,n,r,w)=(256,512,33,42)\), the resulting nondirected diagnostic
variance factor is 31.0225, below the required factor 512 by a factor of
about 16.5. This is the first target-size concentration result for the pure
sparse-mixer constituent. It is not yet a certificate because the one-sided
sums and expected spectrum are not outward rounded. Extending the diagonal
table through level 159 and certifying shells 42 through 79 and 433 through
470 are the next obligations.

The numerical wrapper is now partly outward. Directed binary64 intervals
certify the difficult sector-zero diagonals at levels 80 and 84 by 1.277282466
and 206.661452. The corresponding mixed-arithmetic shell-42 factor is
33.8177. This is not a shell certificate because lower levels still use
nondirected inputs.

The same pairwise-positive relaxation is rejected for later levels. Its
outward bounds grow to \(2.72\cdot10^6\) at level 90 and
\(6.97\cdot10^{12}\) at level 100, while stable signed values remain near
one. Four-message complement grouping reduces the level-100 positive mass to
1.0000032, but small models disprove universal orbit positivity.

Sector zero now has an exact positive replacement. Expanding around the
independent-marginal pair kernel writes its variance as a sum over orders
\(j\). Each order is a quadratic form in the Schur power of the one-row
determinant covariance matrix, hence is nonnegative. A radial Arb evaluation
of those forms now encloses the complete level-100 sector-zero diagonal as
1.000003554708789325353157978375139901... with radius below
\(2.14\cdot10^{-111}\). The merger checks one finite enclosure for every
order from 1 through 512. Repeating or uniformly bounding the calculation
through level 159 is the remaining sector-zero task. Sectors one and two
retain their existing positive outward formulas. An exact rational verifier
checks 2,304 pair-shell identities and 72 fixed-order aggregates in the
\((4,8,3)\) model.

The radial transform also has the exact polynomial form

\[
 F_{j,\ell,p}(t)=2^{k-(n-j)}\sum_{s=0}^{n-j}
 \binom{n-j}{s}K_p^{(n)}(j+s)v_t(s+j-\ell).
\]

The small verifier now checks 1,980 such transforms exactly. The current
polynomial Arb evaluator is correct but slower than the scalar evaluator for
one target level (81.4 versus 44.2 seconds through order 128), so scalar is
again the default. A directed complement-orbit prototype gives a tight valid
upper bound in the small model without assuming orbit positivity. Its Python
target path now quotients complement-and-swap type orbits and constructs the
determinant from its exact integer numerator. It proves
\(D^{(0)}_{100,100}\le2.45522037042628\), covering all 2,862,209 types with
361,985 representatives in 479.6 seconds. This independently confirms that
the sector-zero gate has ample factor-512 slack. The next implementation
target is a compiled or multi-level recurrence across the full band. A
sector-zero trace shortcut was
tested and rejected because low-dual-level collision energy is enormous and
the direct long-double trace calculation is numerically unstable.

### Final one-stage sparse-EA closure

The one-stage route described above is no longer open.  A cancellation-free
dominant-character likelihood argument replaces the failed diagonal shell
relaxation.  It directly certifies the defect shells 42 through 79 and 433
through 470, while an outward central calculation covers every shell 80
through 432.  The high shells are evaluated directly; the proof does not
assume a symmetry of the realized spectrum.

The complete local cap event has positive support exactly on weights 42
through 470.  Each integer cap is the minimum of an independently budgeted
Cantelli cap and the frozen transfer-band envelope.  Conditioning on a
successful rank test and allowing at most 16 independent setup attempts gives

\[
  \Pr[\text{setup failure}] < 2^{-41.3632966090}.
\]

For every accepted constituent satisfying those caps, the finite
RandomStepConv-M22 transfer gives

\[
  \Pr[d_{\min}<228{,}590\mid\text{cap event}]
  <2^{-51.6422972013}.
\]

Adding setup and transfer failure proves the final bound

\[
  \Pr[d_{\min}<228{,}590] < 2^{-41.3621359295}.
\]

The theorem and its exact probability space are stated in
`ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md`.  The hash-bound artifact
inventory is `ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json`, and
`audit_one_stage_sparse_ea_certificate.py` checks range coverage, independent
sector-two agreement, exact small-model tests, cap domination, receipt hashes,
and the final 40-bit gate.

The central-shell certifier now supports process-level parallelism through
`--workers`.  On Peach, eight worker processes with one FLINT thread per
process evaluated the 353 central shells in 84.7095 seconds, versus
771.2017 seconds for the serial local run, a 9.1041-fold speedup.  The two
runs agree exactly on all shell rows and the aggregate claim.  The validation
receipt is
`pure_expander_accumulate/dominant_character_deviation_peach_workers8_validation.json`.

The original Block Expand--5 factor-two route remains a useful historical
alternative, but its length-512 variance lemma is still conditional.  The
one-stage sparse-EA result is the closed finite theorem to carry forward.

An exact-circuit experiment materially lowers the implementation estimate.
On one sampled degree-33 map, verified common-subexpression synthesis reduces
the expander from 16,384 to 8,966 forward XORs. Including its accumulator
gives 9,477 XORs, versus the raw estimate of 16,895. The transposed total is
9,733. This is an exact circuit identity and a heuristic cost result, not a
wall-clock benchmark or a realized-spectrum certificate for the fixed
sample.

The matched Peach benchmark is now complete. At \(k=2^{20}\), both paths map
\(2^{21}\) 128-bit input blocks to \(2^{20}\) output blocks. Two independent
31-sample runs give mean-of-median times of 8.703 ms for sparse-EA and
4.896 ms for two optimized BCH constituents. The sparse path is 1.777 times
slower and adds 3.807 ms to the isolated outer transform. Both generated
circuits and their packed C++ implementations pass exact reference checks.
The benchmark does not include routing or the inner code.

A depth-two \((17,7)\) sample costs 9,300 optimized forward XORs including
both accumulators. It therefore saves only 177 XORs, or 1.87%, relative to
the optimized map from the already-certified one-stage ensemble. Its
realized-spectrum proof remains open. A degree-7 square-mixer diagnostic
improves the degree-5 covariance relaxation by 85.25-fold but still gives a
variance/mean upper bound of \(1.31\cdot10^{90}\) at shell 80. Conditioning
the square mixer on invertibility exactly preserves the dense uniform
reference, but a new conditional sparse-defect lemma is required.

The benchmark rejects the current degree-33 implementation as a direct
performance replacement for BCH. The smallest useful next task is now:

1. search the one-stage ensemble jointly for XOR count and transposed live
   pressure, then benchmark only the best materially different circuit;
2. resume depth two with an invertible-conditioned square mixer and prove
   contraction of the sparse-reference defect;
3. replace RandomStepConv-M22 by RM2Sub-S19 while preserving the proven outer
   cap interface; or
4. replace the uniform routing average by the frozen structured routing.

None of those extensions is part of the theorem proved here.

### Sparse-EA lane freeze

The completed sparse-EA work is preserved by `SPARSE_EA_FREEZE.md` and the
hash-bound `SPARSE_EA_FREEZE_MANIFEST.json`. The freeze contains the closed
10.9000205994%, 41.362-bit RandomStepConv-M22 theorem; the open depth-two
concentration obligation; the optimized XOR circuits; and the matched Peach
benchmark showing 8.703 ms for sparse-EA versus 4.896 ms for two optimized
BCH constituents. `audit_sparse_ea_freeze.py` verifies all 23 frozen files
and reruns the semantic certificate audit. No files were relocated.

### Active RM(4,9) route

The new deterministic outer is the direct sum of 4096 copies of one fixed
RM(4,9) \([512,256,32]\) constituent. Its exact authenticated spectrum has
SHA-256
`995aab561da18f22074b5c6f5413882f492084510aa1cd19b848355f1fcd4ed7`.
The probability space now contains only the uniform routing and
RandomStepConv maps; there is no outer sampling event.

The first outward certificate covers occupation one at the 10.9%, M22
target. It proves

\[
  Q_1 < 2^{-26.5608944032},
\]

dominated by the 52,955,952 weight-32 local words. The result is a proved
upper bound for Q1, not a full distance certificate. It also does not show
that the true Q1 probability is this large. Memory-only diagnostics rise to
only 34.5109 bits at M40, so the current transfer is unlikely to reach 40
bits by enlarging state alone.

The active record is `rm49_outer/README.md`; the outward Q1 receipt is
`rm49_outer/rm49_q1_randomstepconv_M22_d109_outward.json`. The next proof
decision is distance/margin versus a stronger weight-32 argument. Q2 and the
high-occupation transfer remain open and should not be run until that gate is
resolved.

### RM(5,11) follow-up

The next deterministic candidate repeats one fixed
RM(5,11) \([2048,1024,64]\) constituent 1024 times.  Its outer has no setup
failure event.  The transfer retains uniform row-coordinate and
transposed-region routing and RandomStepConv-M22.

The exact five nonzero shells below weight 128 are now checked by integer
formulas.  The same implementation reproduces the authenticated RM(4,9)
low coefficients.  Outward arithmetic proves that these RM(5,11) shells
contribute at most

\[
  2^{-83.3743314000}
\]

to Q1 at the 10.9% target.  Weight 64 dominates.  This is a partial Q1
certificate and must not be described as a distance certificate.

The remaining Q1 obligation is an RM-specific envelope for weights 128
through 380.  All constituent mass can be charged at weight 384 thereafter.
Generic Type-II/Gleason nonnegativity and distance-only Johnson packing are
both too weak.  Available asymptotic RM weight bounds have not yet yielded
explicit usable constants at \(m=11\), and sampled spectra are not proof
inputs.

The active record is `rm511_outer/RM511_Q1_ENVELOPE_STATUS.md`; the outward
receipt is
`rm511_outer/rm511_q1_low_randomstepconv_M22_d109_outward.json`.  The next
task is a finite RM-specific cumulative envelope.  Q2 and implementation
benchmarking should remain paused until complete Q1 closes.

### Direct combined-spectrum EA follow-up

`ea_combined_spectrum/` replaces simultaneous realized-shell caps by the
positive transfer-weighted enumerator actually consumed by the SPIN proof.
The probability space still samples one rank-tested sparse-EA constituent and
repeats that same map in all outer rows.

For the existing ([512,256]) block, direct occupation-one weighting reduces
the apparent threshold from degree 33 to degree 23.  The full diagnostic gives
44.826 bits at degree 23.  The equal-message Q2 sector then exposes the reuse
cost: degree 23 gives only 19.174 bits, and degree 33 gives 42.707 bits.

Larger blocks improve the rank-one Q2 sum while preserving leading raw work
(Nr).  Positive low-weight partial sums nominate degree 31 at block 1024,
degree 29 at block 2048, degree 27 at block 4096, and degree 25 at block 8192.
These are candidates, not certificates; high outer weights, distinct-message
rank two, and (Q\ge3) remain open.

The exact combined-functional first- and second-moment identities have been
checked against exhaustive rational enumeration in a small model.  The
recommended next target is ([2048,1024]), degree 29: sharpen and certify the
rank-one witness, then derive a scalable positive majorant for the rank-two
pair law.  Do not benchmark the larger-block variants until those gates pass.

This follow-up is paused by user decision.  Its artifacts are retained as a
negative-result record and a reusable combined-spectrum template.

### Small-message exact-spectrum replay

The exact extended-BCH and RM(4,9) spectra were replayed at smaller message
lengths with one fixed constituent repeated across rows, uniform routing, and
RandomStepConv. Extended BCH \([32,16,8]\) and \([128,64,22]\) do not clear
the 40-bit occupation-one gate in the tested range.

RM(4,9) at \(k=2^{13}\), \(N=2^{14}\), 32 rows, and memory 22 is the retained
candidate. Binary64 diagnostics give 40.71735 bits for occupation one and
84.27710 bits for occupation two. A random rate-half reference has 59.9674
bits at this length and less than 40 bits below it, making \(2^{13}\) the
first plausible tested length for the intended random-like proof route.

No all-occupation certificate has yet been obtained. Three positive high-occupation
majorants were too lossy. In particular, the banded Bernoulli method needs an
exact split of the all-zero candidate-input atom before applying the Chernoff
distance factor. The corrected split was run and remains vacuous because live
lowest-band compositions dominate.

The random model is now the requested one-map-reused constituent, sampled
uniformly subject to full rank and compared only at matched block sizes. At
\(B=32\) and 128, \(k=2^{13}\), and \(M=12\), its Q1 margins are -10.529 and
6.180 bits; neither closes, and the matching BCH constituents are slightly
better. At \(B=512\), Q1 has 91.019 bits and the exact Q2 split has 180.187
bits. That larger result is only a control for same-size RM(4,9), not evidence
for the smaller BCH regime. Higher random-512 occupations require
rank/relation-type accounting.
The active route is restricted to BCH, RM, and random outer models;
`small_k_replay/SMALL_K_EXACT_SPECTRUM_STATUS.md` records the probability
space, receipts, negative results, and remaining obligations.

A controlled 10% Q1 comparison now uses the common schedule
\(M(k)=\lceil\log_2 k\rceil+2\) from \(k=2^8\) through \(2^{20}\). The maximum
authenticated margins are 4.107 bits for BCH \([32,16,8]\), 30.014 bits for
BCH \([128,64,22]\), and 38.925 bits for RM(4,9). No structured curve reaches
40 bits. Same-size random controls confirm that the two BCH spectra outperform
their random ensembles, whereas RM(4,9) is far below a random 512-bit
constituent. These remain occupation-one diagnostics, not distance
certificates. `small_k_replay/MATCHED_CONSTITUENT_CURVES.md` defines the
comparison and identifies the receipts.

The completed rate-half extension adds BCH-derived block lengths 8, 32, 64,
and 128; RM block lengths 8, 32, 128, and 512; and random full-rank block
lengths 8 through 1024 by powers of two. The length-64 shortened-XBCH spectrum
was exhaustively reconstructed and matches its frozen SHA-256 manifest. The
published RM(3,7) spectrum was added with its source and integer consistency
checks.

No fixed structured constituent reaches the 40-bit Q1 screen under the common
10% and \(M(k)=\lceil\log_2 k\rceil+2\) schedule. RM(4,9) is best at 38.925
bits. The random ladder first clears at block length 256, reaching 58.745
bits. Its length-512 and length-1024 maxima are 129.461 and 274.879 bits.
These are shell-separated binary64 Q1 diagnostics, not outward-rounded or
all-occupation certificates.

The exact coverage and stopping rules are recorded in
`small_k_replay/SPECTRUM_SOURCES.md`. The BCH-derived ladder stops at length
128 for lack of a complete rate-half spectrum. The RM ladder has no
rate-half length-1024 member. The reproducible results are in
`small_k_replay/RATE_HALF_FAMILY_CURVES.md` and the
`rate_half_family_k_margin_d100` JSON/CSV receipts.

### RM2Sub calibration replay

The rate-half BCH, RM, and random-outer scan has now been repeated with an
explicit RM2Sub inner. An equal-persistence comparison selects 64-bit epochs
over the tested 128- and 256-bit alternatives. The neutral state schedule is
\(s(e)=\max(7,e-4)\) for \(k=2^e\); it matches the earlier
RandomStepConv-\(M(e)=e+2\) reset exponent for \(e\ge11\).

At \(k=2^{13}\), exact-spectrum \(Q=1\) and \(Q=2\) comparisons on BCH
\([32,16,8]\) and \([128,64,22]\) place RM2Sub within a few bits of the
persistence-matched random-step abstraction. Over longer regions, the full
RM2Sub transfer can be stronger because the audited live-output spectrum is
reused at every epoch; reset-rate matching alone does not equate the models.

The persistence-matched Q1 replay gives 33.966 bits for extended BCH
\([128,64,22]\) at its maximum and 41.462 bits for RM(4,9) at \(k=2^{16}\).
RM(4,9) stays above 40 diagnostic bits for \(2^{14}\le k\le2^{18}\), but falls
to 38.941 bits at \(k=2^{20}\). These are nearest-binary64 occupation-one
diagnostics, not distance certificates. Outward rounding, occupations
\(Q\ge3\), a realized random-outer event, an arbitrary-length wrapper, and an
XOR-cost measurement remain open.

The controlled record is
`small_k_replay/RM2SUB_CALIBRATION_AND_FAMILY_REPLAY.md`. The audit
`small_k_replay/audit_rm2sub_calibration_replay.py` rechecks \(BA=0\), all
locally generated spectra, exact MacWilliams duality, schedule admissibility,
receipt consistency, and numerical anchors; its receipt reports `PASS`.

### Ten-percent RM2Sub parameter study

The active optimization target is now 10% relative distance with 40 bits of
total failure margin. A 10.9% or 11% result is a bonus, not a selection
constraint. `small_k_replay/TEN_PERCENT_PARAMETER_STUDY.md` defines two
outputs: the fastest fully certified construction and the fastest conditional
construction under one explicit realized outer-spectrum cap.

The first tranche uses one fixed outer constituent repeated across rows. It
tests RM2Sub epoch lengths 128 and 256, persistence exponents 20 through 26,
and native message lengths `2^14` and `2^16`. Exact RM(4,9) at `k=2^16` is the
only exact-spectrum 40-bit Q1 survivor. Its complete 121-witness Q1 margins
range from 41.406 to 42.991 bits over the tested inner configurations.

The corrected Q2 screen gives 85.947 bits for `t128_s13`. An earlier split-grid
calculation reported -11.587 bits, but it combined completed aggregate margins
instead of taking the pointwise minimum over the union of tilt witnesses. That
elimination was invalid. The tested `t256` states from `s=12` through `s=18`
have 83.721 through 88.737 Q2 bits. Occupations at least three, outward
rounding, the final occupation sum, implementation equivalence, and runtime
measurements remain open.

At fixed `k=2^16` and persistence exponent 20, the corrected Q1/Q2 pairs are
42.583/87.424 for `t64_s14`, 42.147/85.947 for `t128_s13`, and
41.406/83.721 for `t256_s12`. Increasing the epoch length decreases both
screened margins. A nested `t=128` family also gives smooth adjacent-state
behavior. The current frontier is therefore `t64_s14`, `t128_s13`, and
`t256_s12`, with larger states retained as safety alternatives.
The consolidated machine receipt is
`small_k_replay/rm2sub_primary_frontier_summary_d100.json`. The exact map and
receipt-structure audit passes in
`small_k_replay/rm2sub_primary_tranche_audit.json`.

`small_k_replay/RM2SUB_PERSISTENCE_UNCONFOUNDED.md` gives the correction,
formalizes the pointwise witness order, and records the nested-family and
matched-persistence experiments.

### Fixed-RM RM2Sub occupation ladder

The matched-persistence replay now covers occupations 3 through 8. It uses
one fixed RM$(4,9)$ constituent in every row and sums exactly over the number
of rows containing the unique all-one codeword. A positive two-colour
hypergeometric recurrence replaces the numerically unstable repeated finite
difference. Its independent direct-enumeration audit passes with maximum
absolute error $1.6654\times10^{-15}$.

At 10% distance and $k=2^{16}$, `t64_s14` is stronger than `t128_s13` and
`t256_s12` at every tested occupation from 3 through 8. Its margins are
82.414, 110.168, 135.423, 166.847, 174.884, and 217.458 bits. The extended
sparse ladder remains positive through $Q=29$, with 30.683 bits there, and is
vacuous at $Q=30$. A 0.01-spaced local Chernoff grid confirms
that the boundary is not a coarse-grid artifact. The no-all-one face is
dominant.

The exact Holder reduction and recurrence are proved, but the numerical
margins are binary64 diagnostics. The construction still lacks a bridge for
occupations 30 through 256 and an outward-rounded final sum. The next target
is a finite three-/four-state dense transfer overlapping the sparse bound at
$Q=29$.

The first dense probes narrow that target. A one-shot full-spectrum Holder
bound is vacuous. With jointly optimized spectrum bands, every pure face of
weight at least 64 is positive at $Q=29$. The minimum-weight dense face has
-29.391 bits at $Q=52$, 4.780 bits at $Q=53$, and 39.057 bits at $Q=54$.
The exact sparse minimum-weight
face has 1224.573 bits at $Q=29$ and stays above 2,100 bits through $Q=54$.
The first mixed-face claim is now available diagnostically.
`small_k_replay/RM2SUB_TWO_BAND_BRIDGE.md` gives an exact positive reduction
for compositions supported on the minimum-weight band and one partner band.
Separate Bernoulli references at 0.25, 0.5, and 0.75 give 1172.882 bits for
the union of all displayed two-band families over Q=30 through 52. The direct
collapse audit passes with error below $5\times10^{-15}$.

This is not an all-occupation certificate. Compositions using three or more
bands remain open in Q=30 through 52, as does the dense union for Q=53 through
256. The calculation is nearest binary64 and is not outward rounded.

The later refined-band calculation supersedes that scope statement. It uses
weights 1--95, 96--416, and 417--512 with Bernoulli references 0.25, 0.5,
and 0.8677722630069483. The positive three-binomial reduction checks every
composition for Q=30 through 256. The consolidation audit passes over
2,857,249 compositions and 227 occupations. Their union has 1228.638
diagnostic bits, with Q=30 weakest.

This closes the missing multiband and dense region diagnostically. It does
not yet close an end-to-end certificate because the arithmetic is not
outward rounded. The controlled account is
`small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md`.

The diagnostic combination has since been completed. The full Q=1--256
union has 42.577982 bits and clears the comparison line by 2.577982 bits.
Q=1 is the bottleneck. The combined audit checks exact occupation coverage.
Only directed outward rounding now separates this finite instance from a
formal 10-percent, 40-bit distance certificate.

Supporting material:

- `small_k_replay/RM2SUB_Q_LADDER.md`;
- `small_k_replay/evaluate_rm2sub_q_ladder.py`;
- `small_k_replay/audit_rm2sub_two_colour_recurrence.py`;
- `small_k_replay/rm2sub_two_colour_recurrence_audit.json`;
- `small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md`;
- `small_k_replay/rm2sub_refined_band_bridge_q30_q256_diagnostic.json`;
- `small_k_replay/rm2sub_refined_band_bridge_audit.json`;
- `small_k_replay/rm2sub_full_occupation_q1_q256_diagnostic.json`;
- `small_k_replay/rm2sub_full_occupation_q1_q256_audit.json`;
- `small_k_replay/audit_rm2sub_q_ladder.py` and
  `small_k_replay/rm2sub_q_ladder_audit.json`; and
- `small_k_replay/evaluate_rm2sub_fixed_rm_dense.py`,
  `small_k_replay/probe_rm2sub_fixed_rm_dense_bands.py`, and their diagnostic
  receipts;
- the `small_k_replay/rm2sub_q_ladder_*_d100.json` receipts.
# Current addition: activation-aware parameter study

The parameter lane now has 715 preferred Q1 observations for exact-spectrum
BCH/RM constituents and random-outer references, across five message lengths
and thirteen nested-map RM2Sub configurations. See
`landscape_db/ACTIVATION_PARAMETER_STUDY.md` for the grid, findings, commands,
and next experiments. All changes are confined to finite_asymptotic_theory.

New artifacts include the activation-aware batched Q1 evaluator, sequential
pilot and grid-refinement drivers, exact-arithmetic tests, hashed map/CSV
receipts, comparison reports, preliminary operation counts, and plots.
The database schema, importer, query output, README, and historical fit status
were updated. Existing proof receipts and frozen implementation sources were
not changed. Database source IDs and original result IDs remain stable.

The old RM(4,9) certificate entries remain indexed but are excluded from the
current-certificate view pending activation review. BCH-256 is not a primary
fit input. No new full certificate or benchmark is claimed. Fourteen tests
pass, including exact support enumeration, a GF(4) first-moment comparison,
database coverage checks, and rejection of a changed screen CSV.

Next: refine smaller s within these chains, then compare activation-aware Q2
and higher occupations before fitting the scaling frontier. The dependency
on the separate BCH task is read-only review/import of its completed compact
certificate; no files or running sessions in that task were changed.
