# RM2Sub integration handoff — 2026-09-05

## Current direction: larger state, positive proof first

The user explicitly authorized larger state and requested a proof before
optimization. Read `LARGER_STATE_PROOF_FIRST.md` first. Independently audited
snapshots t128_s19 and t64_s20 are in generated/larger_state_inputs_v1.
The initial t64_s20 Q1 screen has margin50.439 bits. New numerical values
are discovery only until their outward producer and replay pass.
No old source or certificate was edited, and no other worktree was written.

**FULL t64_s20 DISTANCE TARGET CLOSED.** Read `T64_S20_FULL_CLOSURE.md`.
Every occupancy1..8192 is certified and replayed. The exact setup-failure
upper bound is below2^-50 (diagnostic margin50.43906854330953bits), stronger
than the requested2^-40. The final ledger is generated/larger_coverage_full.json;
the final audit is generated/larger_t64_s20_full_closure_audit.json.
The audit reaggregates every occupancy, rebuilds the1163-row outer joint
LP, replays all46 exact shell caps and OA29 caps, and passes nine tests
in six modules. Q1 has an independent coefficient implementation; range
replays use512-bit Arb plus the checked directed recurrence.

New map tests also check every generator's algebraic normal form has
degree at most2. Seven unit tests across the new endpoint tests and
existing polynomial/scaled/shared recurrence suites passed.
The exact-all-one-band screen (v2) sets p=1,cost=1. An earlier tilt8
adaptive bound failed at Q8192, but tilt6 closes8026..8192. The coarse
gap driver closed4097..8025 before failing at8026; all failed bounds are
retained and excluded. `finalize_gap_shards.py` replays the partial run
and extracts its complete passing interval. Constant-box fallback was
prepared and tested but not needed. No second moment enters this proof.

Next: independent proof review and compact proof packaging, then benchmark
the proved t64_s20 baseline before optimizing parameters. The theorem uses
fresh independent nonzero multipliers; PRG replacement and full SPIN
protocol security remain separate obligations. No performance is claimed.

All earlier sessions and this turn's34066,89943,71368,48657,97830,49563,55943
finished. No known job remains live. Generated evidence is kept compact
with dyadic rounding; nothing was committed or pushed.

The old fixed t128_s15 goal is superseded as the immediate direction.
The app goal API returned null during the latest status check; historical
statements below that it is active are not current app-state claims.
All old positive coverage remains specific to t128_s15.

The previously listed session42044 finished successfully: all six low and
three mid atlas tiles replayed at512bits. The final receipt is
generated/overlap_envelope_201_replay.json. Do not restart that session.
Retiled weighted screens remove one old numerical obstruction but leave
intermediate-total bounds inadequate. No full second moment was proved.

New sources: larger_state_maps.py, screen_larger_state.py,
certify_larger_state_q1.py. Once a source is hashed in a receipt, preserve
it and use a new version for changes. Keep computations sequential.

## Latest: sharper row law, wider atlas, and a proved limit of the old comparison

Previous and current goal turns made verified progress. Goal remains active,
fixed construction unchanged, complete positive coverage still Q1..1655
plus the retained dense partial classes. No actual failure-probability
lower bound or full positive guarantee has been obtained.

Read `JOHNSON_OVERLAP_LAW.md` and `CONVEX_ROUTE_OBSTRUCTION.md` first.
The latter changes the next action: optimizing unweighted h(k) further
cannot close the current comparison, even with perfect numerics.

- `certify_overlap_extremal_law.py`: all existing63 atom caps and mean25
  admit a greatest convex law with boundary atoms5,35. Variance225.584.
  Exact replay and all81 stop-loss chord duals passed.
- `johnson_overlap.py`: exact normalized Hahn kernels. Primary mathematical
  source is Chailloux--Debris-Alazard arXiv2405.07666v2, Proposition8 and
  section4.2. This adds valid Johnson-scheme positivity for the actual
  constant-weight T80 pair distribution. Small exact projector tests passed.
- Numerical conditioning mattered: the unscaled80-degree LP reported
  variance47.75, but its exact dual gave5375.74. Entries<=1e-9 omitted by
  the numerical solve explain the discrepancy; zeroing those entries
  reproduces47.75. This failed gate is retained and replayed in
  `johnson_overlap_variance_attempt01.json`; do not cite47.75 as proved.
  Cap scaling gives a checked degree20 variance bound43.1966; higher-degree
  optimistic values again have worse exact dual bounds. Sources frozen.
- `certify_johnson_convex_law.py --verify`: exact replay of316 Hahn duals
  for79 stop-loss functions, followed by their lower convex hull. The hull
  still dominates every actual convex stop-loss curve. The resulting
  mean25 comparison law nu_J has variance42.39813614903227 and bounds
  EVERY convex expectation, not just variance. All masses, slopes, mean,
  and81 stop-loss identities checked. No optimizer used during replay.
  `test_stop_loss_hull.py` and prior convex-reduction toys passed.
- Six low-overlap producers passed the zero-intercept quadratic gate:
  all even x,y740..900 and k0..39. Combined with the old atlas this covers
  1056321 types for k0..160. Zero overlap is included. Domain ledger passed.
- Three mid-overlap producers cover k161..200. Their zero-intercept gates
  FAILED, but they certify valid positive overlap-dependent intercepts,
  max1.4861787371337414. `certify_overlap_envelope_201.py` combines them
  into a1318761-type bound: log R<=b_k+(3/25000)(k-xy/8192)^2, b_k=0
  through160, and stored per-k b_k through200. Do not use the stricter
  `verify_extended_overlap_atlas.py --stage all`; its gate intentionally
  rejects these mid tiles. Use the201-envelope ledger/replay instead.

**Live job: exec session42044**, command
`python -B workstreams/bch_rm2sub_bridge/certify_overlap_envelope_201.py --replay`.
It sequentially replays low tiles0..5 and mid tiles0..2 at512bits.
On completion it writes `generated/overlap_envelope_201_replay.json`.
Poll this SAME handle; do not restart after an observation timeout.
All earlier sessions27311,3400,7621,48369 are terminal. No other job runs.

The most consequential result is `certify_convex_route_obstruction.py`:
256-bit producer and512-bit replay prove the old unweighted comparison
is >2^23000 (log2 lower~23189.64854), for ANY covering log-convex h,
including separate h_A for each total overlap. It uses identical core
supports and artificial nu_J composition1592*80+1000*61+18*60=189440,
then the multinomial atom with all256 counts740. Every unweighted h has
h(740)>=R(740,740,740)=1/beta740. This is an obstruction to the relaxed
BOUND, not evidence of actual setup failure and not completion of the goal.

One identifiable loss: repeating marginal type(740,740,740) violates the
actual fixed marginal totals208800 per core message. The old unweighted
comparison forgets those totals. A valid next repair is, for each A,

    R(x,y,k+d) <= exp(b_A*(x+y-2*815.625))*h_A(k), with b_A<=0.

The exponential factors multiply to<=1 because correction rows only add
weight. This permits lower h_A(740) without paying a positive global cost.
Use actual-MGF bounds via nu_J, not a direct distribution replacement
for the potentially nonconvex function A->f_{h_A}(A). The notes prove
that per-A envelopes plus actual-MGF Chernoff bounds are valid.

Next: finish/review session42044, then screen jointly optimized b_A and
h_A using retained full type/tilt proposals and the fixed marginal totals.
All overlaps must remain covered. Stronger actual pair caps at difference
weights38,40 are another possible repair: current shell caps exceed
random-spectrum counts by roughly21 and18 bits respectively. Those
random-spectrum values are diagnostics, not substitutes for proved caps.

Additional retained screens: `screen_convex_moment_fit.py` improves the
old LP envelope but is insufficient; `screen_johnson_second_moment.py`
uses the better row law and exact hypergeometric support distribution;
`screen_adaptive_overlap_envelopes.py` uses the new low atlas and separate
h_A, but still fails at large A, as the obstruction now explains.
`screen_adaptive_overlap_201.py` is prepared but deliberately NOT run:
it remains unweighted and cannot solve the identified obstruction.

## Latest: corrected weighted witness and a convex overlap reduction

The goal remains active. Complete positive occupancy coverage is still
Q1..1655, with the previously recorded dense partial classes. No new full
occupancy is closed. The current route tests whether the fixed target can
hold; a successful second-moment bound here would refute it, not prove it.

Read `INVERSE_WEIGHT_AND_PARITY_CORRECTORS.md`, then
`CONVEX_OVERLAP_REDUCTION.md`. The new ingredients are:

- Inverse-probability weighting cancels the two single-message kernel
  products exactly in the conditional second moment. The resulting
  normalized second moment contains only products of actual pair ratios.
- Three existing rows can cancel all256 column parities, except on a
  setup event of probability <2^-100. The exact rank bound uses degree18
  Christoffel caps for Cdual and a triple-intersection union bound.
  A fixed-order right inverse selects actual messages after row setup,
  without changing the encoder or using region randomness.
- A2610-row weight80 core plus those three rows has input weight<=209568,
  occupancy2610..2613, and corrected even region counts740..900 on the
  core window740..897. Normalized weighted mean>=the stored rational
  lower bound approximately0.7703919913290518, hence>3/4. Exact rank and
  512-bit mean replay passed. Independent basis/subcode/right-inverse
  tests passed, including the deficient identical-permutations case.
- Nine outward tiles cover793881 types: even x,y740..900, k40..160.
  All producer intercepts are negative, giving
  log R <=(3/25000)(k-xy/8192)^2. Exact coverage ledger passed.
  Full sequential512-bit replay passed all9 tiles. Session80165 is terminal;
  `generated/overlap_atlas_replay.json` records completion and all hashes.
- An exact convex-order upper law nu for a shared T80 row's overlap has
  mean25, variance~350.137, atoms0,39,40..61,80, and tail mass<.002564.
  `certify_overlap_convex_law.py --verify` passed. This is convex order,
  not ordinary stochastic domination.
- A proved coupling replaces uniform subset contributions by multinomial
  counts for a positive log-convex product envelope h. Its expectation
  f(A) is convex in the total A. This permits nu replacement and gives a
  finite one-dimensional reduction, including all supports and diagonals.
  Optional hypergeometric-to-binomial replacement gives the coefficients
  of (1-p+p sum_a nu(a)z^a)^2610, p=2610/8189. Positive coefficient bounds
  can evaluate both remaining layers with one-variable tilts.
  Four exhaustive toy tests passed, including U-shaped h.
- The earlier uncorrected central witness also improved: degree6 log-beta
  envelope, exact additive moments, and512-bit replay give diagnostic
  log2 mean lower19633.21535078605, about8.53 bits above the quadratic
  bound. This remains a fallback; the weighted witness avoids that cost.

All jobs are terminal. The following discovery screens also completed:

- `screen_global_overlap_envelope.py`: raw positive tilts cover all k0..900
  but are loose. Retained, not outward evidence.
- `pair_positive_tilt.py`: analytic-gradient optimizer, finite-difference
  test passed. Frozen by the following discovery receipts.
- `screen_optimized_overlap_tails.py` and
  `screen_global_overlap_optimized.py`:15 marginal blocks with optimized
  positive tilts cover all k0..900. All tilt proposals retained for a later
  outward pass. Worst screened log ratios at k0,39,161,400,600,740 are
  approximately8.03,10.61,11.99,39.35,193.11,658.76. The central atlas still
  supplies much sharper bounds for k40..160. Session61156 is terminal.
- `screen_convex_second_moment.py`: mean-log LP objective produces a bad h.
  At A20800 its combined coefficient-bound summand has about7027 bits.
  This is slack in a proposed upper bound, NOT evidence of actual failure.
  LP h(0) has log~90.5, and its final slope3 gives excessive tail costs.
- `screen_two_line_convex.py`: improved but still unsuccessful global
  convex proposal. A20800 bound is about127 bits, already above the
  roughly40-bit total target; high-A bounds are much worse. One affine
  line for the entire upper tail is too restrictive. Retain the failure.

Next: construct a better correction-robust log-convex h. A natural cheap
test is the LP envelope with endpoint slopes restricted to about-.7 and
1.8, or multiple upper-tail affine pieces, rather than one line. Optimize
an exponential moment, not expected log h. The present mean-preserving
nu law may also be too loose for large A: consider extending the exact
two-block Christoffel pair caps to difference weights82..160 and using
more row-overlap information. Do not certify a numerically failed global
proposal. After screens succeed, outward-certify both the global h and
the entire finite sum in `CONVEX_OVERLAP_REDUCTION.md`.

The threshold needed to refute the fixed failure target is normalized
second moment <(9/16)*2^40. Neither a large upper-bound screen nor the old
large actual first moment proves that failure probability exceeds2^-40.
Unhandled overlap tails must never be discarded.

Frozen certificate sources and previous receipts were preserved. New
receipts contain only compact tables; Fourier arrays were temporary.

## Latest: 73600 certified pair types and a central-count witness

Previous and current goal turns made verified progress. Goal remains active
and unchanged; no new complete positive occupancy coverage. Read
`PAIR_ALIAS_TABLE_AND_CENTRAL_WITNESS.md` for the new proofs and assumptions.

- `certify_pair_alias_table.py --verify`: 73600 actual shared-region types
  certified, x,y=780..858 even and overlap60..105, all ratio beta2/(beta_x
  beta_y)<41/40. Maximum stored upper1.0244557494297624. One256x256x128
  Fourier alias table, explicit outward radix2 butterflies and input errors.
  256-bit producer and512-bit root/final replay passed. Arrays temporary;
  only46 per-overlap maxima retained. `test_outward_radix2_fft.py` passed
  direct512-bit DFT checks, including zero input errors.
- `zero_state_central_witness.py --verify`: define Z_G to count only actual
  T80 bad messages whose every input region has even weight740..900.
  Construction/setup distribution unchanged. Its expected count exceeds
  2^19624 (diagnostic19624.688243732788), losing10.405423268675321 bits
  against unrestricted Jensen. Uses exact energy mean/variance, degree4
  parity mixing, conditional Cauchy--Schwarz, and integer-checked quadratic
  log(beta) envelope. 512-bit replay passed; no optimizer used in replay.
  This removes outside-window weights from THIS witness's second moment.
- `certify_pair_type_conditioned_parity.py`: conditional on ANY compatible
  single-region pair type with even marginal counts, global pair parity is
  2^-508 times(1 +/- epsilon), epsilon<2^-3200. Uniform over fixed occupied
  supports/bit patterns. Uses conditional character squares, close-pair
  Bayes bound, and Cauchy--Schwarz; actual shared conditional gap<.373748.
  Exact replay and nonuniform toy character-square/Cauchy tests passed.
  Arbitrary nonnegative SINGLE-type weights may be averaged against this
  estimate; products across different regions are still not covered.

Retained screens: small128x128x64 alias grid was too loose (max~15.94);
larger grid succeeds. Wider three-tilt atlas covers x,y740..900, k40..126
only diagnostically. Edge aliases and cancellation matter; imaginary
residual up to~.006 relative. Targeted probes: type6838:614:614:126 ratio
~1.302 remains genuinely larger; type6514:860:778:40 improves from~1.317
to~1.0115 after retilting. Do not promote wide-screen maxima to certificates.
All jobs are terminal. Certificate sources are now frozen by receipts.

Next: extend outward types across Z_G's full740..900 window, with an
overlap-dependent envelope and explicit tail control. Then integrate the
joint kernel factors across regions. Neither multiplying41/40 across256
regions nor dropping unhandled overlap tails is justified. Full setup
failure upper/lower probability remains unproved.

## Latest: weighted parity and a sharper actual denominator

The previous and current goal turns made verified progress. Goal remains
active, fixed parameters unchanged, no new complete occupancy closed.
Read `WEIGHTED_PARITY_AND_JENSEN.md` for complete derivations.

- `zero_state_jensen_lower.py --verify`: 256-bit producer and512-bit replay
  passed. Exact-weight80 actual zero-state probability lower log2 is
  -244964.3269344031; mean lower log2 is19635.093667001463. This gains
  41.80401552073644 bits over the earlier T80 mean. Use this denominator
  for the next second-moment attempt. It uses conditional one-column parity
  (uniform79/80 slices on255 coordinates), then Jensen on the actual joint
  region-count law. No region independence is assumed.
- `certify_weighted_parity.py`: exact polynomial weighted pair-parity lemma.
  Variables are averages of [0,1]-valued row-local functions; degree d marks
  at most d rows. For d<=512 and coefficient-L1/mean<=2^2048, the relative
  all-even pair-parity error is<2^-500. Receipt replay and exact toys passed.
- `certify_aligned_setup_cluster.py`: admissible identity-permutation setup
  contains >2^64299 T80-family bad messages atQ2620. The B kernel quartet
  {0,4,65,69} has32 directly verified disjoint translated quartets. Repeat
  one arbitrary T80 word within each of655 selected quartets. This explains
  why a universal max-count cap plus the retained mean lower cannot give
  useful failure probability; even the sharper mean/cluster is<2^-44000.
  NOT a probability refutation. Frozen receipt and replay retained.
- `test_conditional_column_parity.py`: exact exhaustive conditional Fourier
  identities. `test_weighted_parity.py`: shared-permutation marked-row checks.

No process is running. New certificate sources are frozen by their receipts.
Next: bound the normalized joint kernel factor with polynomial enclosures
or an exponential-moment argument with explicit remainder. Weighted parity
can now be retained for a defined approximation class. Its theorem does not
automatically cover the true kernel factor. Atypical types remain essential.

## Latest: actual pair coefficients, parity mixing, and fixed-weight moments

Read `PAIR_FOURIER_AND_PARITY.md` first for the new derivations and scope.
The goal is still active and unchanged. No new complete occupancy is closed.
The first-moment obstruction is now stronger, not removed.

Completed, replayed results:

- `certify_pair_type_fourier_tracked.py --tag central_tracked --verify`:
  beta_2(6638,736,736,82)/beta_818^2 <=1.0032168417224794.
  Discrete Fourier alias/L1 bound, with point-dependent outward binary64
  error and 512-bit root/final replay. No iid substitution.
- `extend_pair_type_box.py`: 125 exact integer types, weights814..822 even,
  overlaps80..84, all ratios<1.101. Maximum1.1004688048850855. The attempted
  1.1 threshold FAILED; preserve the false receipt flag.
- `certify_tail_parity_mixing.py`: joint all-even column parity probability
  is 2^-510 times (1 +/- epsilon), epsilon<2^-3100, uniformly over occupied
  supports of size2620. This is UNWEIGHTED and excludes other B constraints.
- `zero_state_affine_lower.py --verify`: original T38..80 actual mean
  exceeds2^20423 using an integer-checked affine envelope of log(beta_j).
- `certify_weight80_family.py`: actual T80 family, exact inputweight209600,
  shell lower log2=98.16756486701678; mean lower log2=19593.289651480736;
  close-pair upper diagnostic bits8.607435602369804. Repeated replay passed.
- `certify_weight80_fourth.py`: exact cap/dual majorant bounds the mixed
  fourth tensor of shared-permutation T80 rows. tau<=.003273915558 and
  normalized shared-support fourth correction<1/36000. Replay passed.
  This controls moments through degree4, NOT the full second moment.
- `verify_weight80_dependence.py`: independent exact-Fraction reconstruction
  of all63 difference caps and the dual, without optimizer or floating point.
  It also proves T80's own pair-parity relative error<2^-3300, uniformly over
  occupied supports, and checks the125-type ledger's <1.101 threshold.

New exact tests: Fourier alias identity and roundoff budgets, fixed-weight
fourth tensor. Complex512 spot checks supplement the Fourier error proof.
All passed. No process is running; preserve all frozen source bindings.
Current sources and generated receipts are isolated in this workstream.

Next: use the exact-weight80 family (p=5/16) to bound the weighted joint
kernel event across256 regions. Two message count vectors have zero cross
covariance and mixed third moments; the fourth correction is explicit.
Higher-order terms and atypical types remain uncontrolled. Do not multiply
the central type ratio256 times and claim a full variance bound. Do not
restart attempts to upper-bound the unconditional first moment below1.

## Latest: exact second-moment ingredients, not a probability conclusion

The previous goal turn and this turn both made verified progress. The full
fixed-s15 goal remains active; its unconditional first-moment route is
rigorously ruled out by `ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md`.
Read the new `ZERO_STATE_SECOND_MOMENT.md` for the actual next obligation.

New completed/replayed results:

- `inner_pair_spectrum.py`: exact fixed-A pair spectrum, 155 nonzero triples,
  all2^30 ordered pairs, from integer Walsh transforms. Receipt
  `generated/t128_s15_pair_spectrum.json`.
- `bch_translation_symmetry.py`: checks1048 basis images for8 additive
  translations. P and its p37 fibers are preserved. The fixed C and its
  weight-defined tail T are coordinate-transitive. Receipt
  `generated/bch_translation_symmetry.json`.
- `certify_bch_tail_overlap.py`: bivariate OA29 Christoffel caps prove
  Pr[wt(U+V)<=80]<1/128 for independent uniform U,V in T (weights38..80).
  256-bit producer and512-bit replay passed; margin7.42979623196436.
  Receipt `generated/bch_tail_overlap_outward.json`.

New exact toy tests all pass: `test_inner_pair_spectrum.py`,
`test_pair_kernel_iid.py`, `test_bivariate_christoffel.py`, and
`test_single_region_invariance.py`.

Important consequence of transitivity: for the FULL message family with
uniformly chosen Q occupied row positions, the fraction of messages passing
any ONE region's kernel checks is identical for every setup. Its one-region
pair probability is the square of that fraction. This does NOT factor over
all256 regions; exact toy counts demonstrate that distinction.

`pair_kernel_iid.py` evaluates the pair MacWilliams polynomial. It is only
an IID evaluator, not the actual shared fixed-type region probability.
The latter is the four-type coefficient in equation(5) of the new note.
No full second moment or setup-failure lower bound has been obtained yet.
No process is running. Keep all frozen source hashes and generated receipts.

Next: a useful bound for the shared-permutation pair coefficient, integrated
over the outer rows without discarding cross-region dependence. Do not
resume attempts to make the unconditional bad-message first moment small.

## READ FIRST: actual first-moment obstruction at Q2620

Read `ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md` before doing more upper-bound
sweeps. The actual fixed t128_s15 construction has expected bad-message
count >2^19600 at Q2620. This is a verified LOWER bound, not an upper-bound
failure or random-code approximation. It rules out the unconditional
first-moment method for the full range. It does NOT alone refute setup
failure <2^-40. The full goal remains active, unchanged.

Exact BCH tail lower count for weights38..80 is
443405397513809550622719915749 (log2 98.4845410838). With Q2620 even,
every such message has input weight<=209600. Fourier inversion lower-bounds
all-even region counts by2^-255. After explicit binomial tails, each region
has even weight64..1536 with sufficient probability. The exact kernel
polynomial gives per-region zero-syndrome probability>2^-959. Together,
the averaged zero-state probability is >2^-245760. On this branch Y=X.
Multiply by binom(8192,2620)*tail_count^2620 for the actual first moment.

The exact tail solve78132 and replay60703 are TERMINAL, successful.
`exact_tail_lower_bound.py --verify` passed rational primal/dual replay.
`replay_zero_state_lower.py` passed512-bit replay and wrote the combined
obstruction receipt. `verify_zero_state_obstruction.py` gives the simplified
exact ledger. Both exact parity tests pass. All sources/receipts are frozen.

Next recommendation: investigate a second-moment/failure-probability bound
for this zero-state family. Do not spend more turns merely trying to drive
an unconditional first-moment upper below1: it cannot happen. Do not silently
change s=15 or mark the goal complete; probability concentration is still open.

Other calculations are TERMINAL: serial middlefill60339 completed every
even shell82..126 that was missing. All even weights38..128 now have exact
caps, with complement caps too. Shared-mode35329 and finer-mode17806 ended.
The finer fixed-weight screen passes Q1792 (~2408bits) but fails Q2048 and
above; it was not certified. Full exact-cap screen50529 likewise completed
and remains vacuous at Q2048/2620. No computation is left running here.

## Newest: full occupancy coverage through Q1655

This continuation made verified progress, not a no-progress turn. The goal
is still the full Q1--Q8192 range, unchanged. `verify_constant_split_coverage.py`
now verifies Q1--Q1655 plus the previous dense classes, combined <2^-49.

New complete range Q1252--Q1655 contributes <2^-79, diagnostic79.9460151154.
`full_exactcaps_q1252_q1792_outward.json` stores all attempted rows. Keep
Q1656--Q1792, but exclude them from the covered sum. The producer (12765)
and 512-bit replay (31104) are both TERMINAL, successful. Replay command:

```powershell
python -B workstreams/bch_rm2sub_bridge/certify_full_exactcaps_range.py --screen full_exactcaps_middle01_screen.json --lower 1252 --upper 1792 --tag q1252_q1792 --verify
```

The full thirteen-band screen identifies pure-band loss around Q2048.
`single_shells_q2048_first_screen.json` further shows vacuous individual
shells even with no assignment cost (e.g. w96, tilt -8: about -15352 bits).
Thus shrinking bands alone cannot fix those witnesses. Screen71557 is
TERMINAL. Further improvements require stronger caps or a tighter moment
comparison, not a claim that the code itself fails.

Two calculations were launched; poll their exact handles for liveness:

- exec35329: `screen_shared_mode_slopes.py --tilt -10 --grid 32 --tag m10_g32`.
  Common negative rates on early segments avoid the previous activation-entry
  fit's growing modes. It is a SCREEN ONLY. Last observed: polynomial powering.
- exec60339: serial exact shell sweep, tag middlefill, weights
  82,84,86,88,92,94,96,98,104,106,108,112,114,116,118,122,124,126,
  300 seconds per solve, warm `joint_shell_chain01_w80`.

Next: inspect the mode-fit result, finish the exact middle caps, and use
either improvement for new shared range witnesses. Frozen producers and
receipts must remain untouched; all generated data stays ignored by Git.

## Newest continuation: constant rows separated; partial dense coverage

The goal remains active. Read `CONSTANT_ROW_SPLIT.md` first. Write only in
this worktree; preserve frozen source hashes and generated receipts.

New outward bounds, all replayed at 512-bit precision:

- Every Q1025--Q1251, including ALL row mixtures: combined <2^-64.
- d<=1200 ordinary rows, h>=513 all-one rows: combined <2^-5390.
- (d,h)=(8192,0), all rows ordinary: combined <2^-9382.

The older d<=512,Q>1024 bound (<2^-67182) also replayed, but is subsumed.
Q1--Q1251 plus these new classes still has combined contribution <2^-49.
The full goal is NOT closed: intermediate compositions remain open.
`verify_constant_split_coverage.py` checks the exact partial ledger.

All shell sweeps mentioned below are now TERMINAL. Successful exact caps
include every even weight 38 through 80, plus 90,100,102,110,120,128 and
their complements. Sessions 5875 and 20812 completed. The old weight-48
attempt in session 53646 failed; its fresh lowretry attempt succeeded.
Preserve both. Replay sessions 64775 and 54989 completed successfully.

The first middle-composition refinement (session 3756) completed without
positive margins. The weight-conservation single-exponential screen also
completed with vacuous bounds. The exponential-mode route remains a proved
conditional formula without a useful numerical envelope. Do not treat
floating screens or successful pure-band probes as certificates.

Screens 55369 and 42470 are TERMINAL. The lower tilts improve the ordinary
probes: d1025,h0 passes at -16 and d1536,h0 passes at -12. These are screens.
The full thirteen-band calculation with the old Q1024 tilt -18 witness
and improved caps then certified all Q1025--Q1251, subtotal margin64.8865.
`exactcap_shared_q1025_q1280_outward.json` retains failed Q1252--Q1280.
Its producer's embedded verifier fails a dictionary-key comparison because
weight zero is omitted from the snapshot. Do not edit the frozen producer.
`replay_exactcap_shared_range.py` fixes the comparison in a separate reader,
recomputes all rows at 512 bits, and has passed.
Next recommendation: seek shared witnesses covering ranges of compositions,
and diagnose the pure-band versus adaptive-mixture gap in the middle.

## Latest continuation: exact shell improvements and a fixed-weight route

Coverage remains Q1--Q1024, partial sum <2^-49; the full goal is active.
`DENSE_RANGE_ALTERNATIVES.md` records the new results and conditional bound.

New rationally certified shell caps: weights 50,52,54,56,58,60,62,128.
Examples: A50 <= about 2^60.775 (previous 2^69.129), A54 <= about 2^62.370
(previous 2^77.122), A128 <= about 2^124.673 (previous 2^126.335).
`exact_joint_shell_caps.py` reconstructs the retained model and checks every
primal/dual relation exactly. Weights 50,52,54,128 were also replayed.
The floating solver's infeasibility reports were numerical, not proofs.

The first sweep, exec session 66824, is TERMINAL: it timed out on weight
64 and has no solution there. Preserve its attempted folder. A retry sweep
was launched as exec session 5875 using:

```powershell
python -B workstreams/bch_rm2sub_bridge/continue_exact_shell_sweep.py --weights 64 66 68 70 80 90 100 102 110 120 --tag chain01 --seconds 300
```

Poll that exact handle to determine current liveness; do not restart from
this note alone. It uses a nearby successful basis and creates fresh folders.
Failure of one objective is retained and does not stop the other objectives.

New algebraic route: if R_j <= sum_i P_i (a+i*Delta)^j entrywise for all j,
then matrix-polynomial mode coefficients and fixed-weight symmetric means
bound every occupancy via the outer enumerator. See equation (2) in the
new note. This keeps row weights fixed and has no Bernoulli density factor.
`test_exponential_mode_bound.py` passes exhaustive exact toy checks,
including noncommuting matrix products. A useful numerical mode envelope
has NOT yet been fitted or certified. This is the recommended next focus.

The convex-sequence alternative was also proved/tested but is numerically
too loose at Q8192 (over 880,000 bits lost in one probe), despite little
loss at some Q2048 tilts. `screen_dense_scaled.py` on Q8192 at tilt 0 remains
vacuous even without convexification. Do not treat that as code failure.
Region logs at Q2048 agreed at 192,512,1024-bit Arb precision.

## Latest: every Q1--Q1024 verified; dense range still open

The active goal is unchanged: fixed BCH [256,128] with RM2Sub t128_s15,
K=2^20, N=2^21, H=209716, fresh independent nonzero multipliers, every
occupancy through 8192, total failure below 2^-40. Do not mark it complete.

`SHARED_RANGE_1024.md` records the newest proof and replay commands.
The added Q129--Q1024 sum is below 2^-77 (77.0174168827 diagnostic bits),
and the total covered Q1--Q1024 remains below 2^-49. All rows were replayed.
`verify_coverage_1024.py` checks the exact overlapping-certificate ledger.
The shared four-anchor batch alone has a bad Q339 bound; the separate
`oa29_gap339_outward.json` closes that one gap. Keep both receipts.

The successful new ingredients are exact OA29 Christoffel shell caps,
thirteen groups with independently selected probabilities, Arb polynomial
matrix powering, and an outward binary64 adaptive recurrence shared across
all depths. All thirteen groups are retained in certification.

Next diagnostics:

- `scaled_adaptive.py` handles dense exponents per degree; exact tests span
  thousands of binary exponent bits. The unscaled recurrence cannot be
  relied on for useful dense bounds because upward underflow becomes loose.
- `screen_dense_scaled.py` has Q2048 probes at tilts -12,-10,-8,-6,-4.
  Their best reported margin is negative. They are screens, not certificates.
- `diagnose_dense.py` finds vacuous pure groups as well, so restricting the
  adaptive maximum alone is insufficient for those witnesses.
- `syndrome_activation.py` uses a valid Fourier bound on each syndrome's
  probability to dominate activation by a multiple of the uniform law.
  It passes exact syndrome and all weighted-prefix tests, but the tested
  Q2048 pure-group bounds still do not all pass. No occupancy receipt uses it.

Recommended continuation: distinguish outer-shell-cap loss from the loss
of hard input-weight information in the Bernoulli comparison, then choose
a tighter dense-range bound. Do not change the fixed construction or claim
that a vacuous envelope is a counterexample to the code. No process was
left running at this checkpoint. Keep all frozen source hashes intact.

## Latest: Q1--Q128 verified; stronger OA29 shell caps available

The active goal is the complete fixed-construction range Q1--Q8192, with
total failure below 2^-40. It is not complete. The exact covered partial
sum Q1--Q128 is below 2^-49.

`tightened_occupancy.py` intersects the old termination-mass cap with the
full emission-moment cap for every input weight, including kernel weights.
The Q65--Q128 sum is below 2^-148 (diagnostic margin 148.8632719886).
All 64 occupancies were replayed at 512-bit precision. Run:

```powershell
python -B workstreams/bch_rm2sub_bridge/extend_adaptive_certificate.py --lower 65 --upper 128 --tightened --verify
python -B workstreams/bch_rm2sub_bridge/verify_coverage_128.py
python -B workstreams/bch_rm2sub_bridge/test_tightened_occupancy.py
```

`christoffel_caps.py` obtains exact stronger caps from the retained strength-29
proof: A_w <= 2^128 / sum_{j=0}^{14} K_j(w)^2 / binom(256,j), with the sum
interpreted as the entire denominator. Weight 70 drops from 128 to about
89.55 bits. The rank certificate is replayed, not assumed from a table.
The derivation is in `OA29_SHELL_CAPS_AND_TIGHTENING.md`.

The new caps are not yet used in an outward occupancy certificate. The
floating `screen_oa29_range.py` keeps a common log-odds shift; its 1024-row
bound remains vacuous. `screen_independent_probabilities.py` instead chooses
each band's probability using a pure-group optimization, then applies the
full adaptive mixture bound. These remain witness searches only.

Next: finish the independent-probability screen, then certify complete
intervals with shared witnesses. Preserve all frozen sources. Dense floating
computations need care about underflow; only outward replays count as proof.

## Latest: composition-free adaptive bound covers Q17–Q64

`ADAPTIVE_OCCUPANCY_RANGE.md` proves the new bound. Distribute each group's
counting cost as Gamma_g^(1/256), then maximize weighted Bernoulli updates
entrywise. The resulting adaptive process dominates every fixed ordered
group sequence. Multiplying by 5^Q covers all assignments without enumerating
compositions or assuming a vertex maximum.

Every Q from 17 through 64 is outward-certified and replayed at 512 bits.
Their sum is below 2^-142; the full covered range Q1–Q64 remains below 2^-49.
The receipt is `generated/adaptive_q17_q64_outward.json`, about 68 KB.
`certify_adaptive_range.py --verify` and `test_adaptive_range.py` are read-only.
The coarse `adaptive_q1_q64_screen.json` has failed diagnostic points; the
successful discovery file is `adaptive_q17_q64_refined_screen.json`.

Next: extend above 64, sharing witnesses across intervals. The composition
enumeration bottleneck is removed, but dense-range performance is not yet
established: 5^Q and adaptive freedom may become too loose. Earlier source
files and certificates were not changed; the other worktree is untouched.

## Latest: general kernel transfer and complete Q4–Q16 batch

`GENERAL_OCCUPANCIES.md` supplies the all-weight epoch envelope, including
zero-syndrome inputs, and a five-group maximum-density counting argument for
any Q in 1..8192. Numerical coverage is now every Q in 1..16, not the full
range. The added Q4–Q16 sum is below 2^-162, and the complete partial sum
Q1–Q16 is below 2^-49.

All 20,293 compositions in the new batch were evaluated with 256-bit Arb
and replayed at 512 bits. Exact toy tests cover every epoch weight, including
kernel weights, with 8,736 prefix comparisons. `general_batch_certificate.py
--verify` is read-only. The roughly 75 KB receipt is
`generated/general_q4_q16_outward.json`; its witness screens are hashed.
Prior certificates and source files are unchanged.

Next: certified bounds over ranges of compositions and occupancies. The
five-group formula has about 1.88e14 compositions at Q=8192, so extending
enumeration alone is not a practical closure strategy. Do not assume a
vertex maximum: some observed worst compositions are mixed. Coarse bands
may also need refinement in the dense range. The all-occupancy formula is
established, but the Q17–Q8192 numerical bound is still open.

## Latest: occupation three completed

`OCCUPATION_THREE.md` closes Q3 for t128_s15 with U_3 < 2^-126 and
diagnostic margin 126.5064403982 bits. U_1 + U_2 + U_3 < 2^-49 is checked
exactly. The method conditions auxiliary Bernoulli rows back to the exact
fixed-weight support law and pays all binomial costs. It does not use the
old M22 weight-deletion monotonicity or a cubic coefficient table.

All 455 weight-group cases passed independent 512-bit replay, including
their rational conditioning costs and final aggregation. The successful
receipt is `generated/t128_s15_q3_compact_outward.json`. The earlier
unoptimized `q3_conditioning_outward` receipt is a retained vacuous bound,
not the current result. Use `certify_q3_compact.py` for fresh generation and
`verify_q3_compact.py` for read-only replay. No Q1/Q2 source or certificate
was changed. Next is Q4, where the kernel's weight-four words must enter
the epoch transfer.

## Follow-up: occupation two completed

`OCCUPATION_TWO.md` records the new t128_s15 result: U_2 < 2^-91, with
diagnostic margin 91.5312913644 bits, and U_1 + U_2 < 2^-49. The region
calculation includes both same-epoch and different-epoch placements. Distinct
B columns exclude kernel events at input weights one and two. The directed
pair recurrence, exact toy tests, and independent dominant-pair Arb check
passed. New files remain confined to this directory. The next step is Q3;
the earlier Q1-only handoff below is retained as historical context.

## Ownership

Our worktree is `ba80/permute_conv`, on branch `codex/bch-m22-proof-notes`.
All integration changes are in `workstreams/bch_rm2sub_bridge/`.
The other agent owns `3061/permute_conv`; we only read files there.
There are no modifications to shared scripts, production code, or the old
BCH bundle. No commit or push was made during this integration.

## Reusable result

The README defines the model, proves the activation-aware Q1 transfer, and
explains how to reuse the existing BCH weighted LP bound by domination.
For K = 2^20 and cutoff 209716, the Q1 upper-bound margins are 49.8589 bits
at t64_s16, 49.4813 bits at t128_s15, and 48.8336 bits at t256_s14.
No RM2Sub claim for Q >= 2 is made.

The source selection and spectrum inputs were snapshotted from
`workstreams/finite_asymptotic_theory/small_k_replay/rm2sub_calibration_constituents/`
in worktree 3061. `inputs/manifest.json` records all source and snapshot hashes.
The bridge independently enumerates each selected A spectrum and verifies its
full B-kernel spectrum by exact MacWilliams arithmetic.

## Transfer issue for cross-review

The snapshot's reference checker is
`small_k_replay/certify_rm2sub_rm49_q1_outward.py`. Its two-state epoch transfer
does not distinguish the state immediately after zero-state activation.
For the documented output-before-update recurrence, q' = B(e_J) when q = 0.
That law is supported on t states, each of mass 1/t. For our configurations,
1/t exceeds 1/(2^s - 2), the asserted near-uniform live-class density bound.

This disproves that distribution invariant at activation. It does **not** by
itself disprove the final numerical inequality: a separate moment-domination
argument could make the same two-state numbers valid. Simple two-epoch checks
for the three selected maps did not find a violated moment bound. Do not
interpret the issue as a counterexample to the RM2Sub construction.

The other worktree already uses a deterministic class in
`analyze_rm2sub_dense_occupation.py`. Our three-state Q1 transfer takes the
same precaution, using the minimum A distance for its arbitrary-live class.
This costs only about 0.009–0.013 bits relative to our historical two-state
screens. We did not alter or recertify the other worktree's RM(4,9) theorem.

## Checks completed

- Exact A/B algebra, complete A spectra, and complete kernel spectra for all
  three selected map pairs.
- 2,976 exact rational prefix checks on a GF(16) toy instance, covering zero,
  every nonzero point mass, uniform-live, and every punctured-live start.
- Explicit enumeration of marked epoch positions and region supports on a
  small instance, compared with the coefficient recurrences.
- Log-domain versus positive transfer checks for the actual maps.
- 256-bit Arb evaluation of all 92 Q1 coefficients for each configuration.
- Independent 512-bit Arb replay of all 276 coefficients, using binary region
  powering and unnormalized support counts; exact rational aggregation agrees.
- The old `verify_bch_m22_closure.py` replay passed, including reconstruction
  of the 1,163-row outer LP and its algebraic constraints.

The proof concerns the precisely stated fresh-multiplier model. It does not
establish seeded-randomness equivalence, measured encoder performance, or the
full SPIN guarantee. The notes distinguish these interfaces explicitly.

## Next integration step

Start with occupation two for the same fixed BCH code and t128_s15. Derive
weight-two epoch bounds that preserve the activation class and account for
B(X) = 0. Use the actual support law induced by two shuffled outer rows;
do not substitute independent Bernoulli bits without a justified comparison.
Retain t64_s16 as a reference, since its Q1 bound is strongest here.

Only after this step should the full higher-occupation machinery be ported.
The old M22 higher-occupation certificate is not an RM2Sub tail certificate.
Keep generated data ignored, preserve write-once evidence, and coordinate any
later shared-interface change before editing the other agent's files.
