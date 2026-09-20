# Sol-to-Astra handoff: finite and asymptotic SPIN distance work

## Astra completion addendum

The principal open task in this handoff is now complete. The fixed-RM(4,9),
RM2Sub (t=64,s=14), (k=2^{16}) instance has a full outward Q=1--256
certificate. It proves (d_{\min}\ge13{,}108) at (N=2^{17}) except with
probability below (2^{-42.5779817562755}). See
`small_k_replay/RM2SUB_RM49_FINITE_CERTIFICATE.md` and
`small_k_replay/rm2sub_rm49_t64_s14_full_distance_outward.json`.

The old open questions and bootstrap steps below are retained as historical
state. They are superseded where they say that outward compilation is still
open. The next active task is the 10% certified/extrapolated parameter
frontier across message lengths and exact-spectrum outer choices.

Snapshot date: 2026-09-04. Repository root:
`C:\Users\peter\.codex\worktrees\3061\permute_conv`.

This document transfers the working state of the
`workstreams/finite_asymptotic_theory/` workstream. It is not a claim that
every file in the worktree is committed, nor that every numerical receipt is
a theorem. The labels **proved**, **audited diagnostic**, **experiment**, and
**open** below are part of the claim.

## 1. Objective

### Ultimate goal

Prove a concrete minimum-distance theorem for **Structured SPIN**, or for a
close variant that preserves the design's useful structure:

- rate near (1/2);
- one constituent code sampled or fixed once and reused at every outer-block
  position;
- the factored route consisting of row-coordinate permutations, bit
  transpose, and region permutations;
- a bounded-memory, linear-time inner such as RM2Sub;
- ordinary and transposed encoding in linear time; and
- competitive concrete performance, with about 11 ms on Peach as the current
  upper edge of an attractive (k=2^{20}) implementation.

The practical distance target is now 10%. A 10.9% or 11% theorem is a bonus,
not a condition for exploration. Forty bits of setup-failure margin is the
standard comparison line, but signed results below 40 bits must still be
retained during parameter exploration.

### Current subproblem

Compile the completed binary64 proof template for the fixed-RM(4,9), RM2Sub
small-(k) construction into an outward-rounded finite certificate. The
current instance is

\[
  k=2^{16},\quad N=2^{17},\quad
  d_{\min}\ge 13108,
\]

where the bad event is output weight at most
(13107=\lfloor0.10N\rfloor). It uses one fixed exact-spectrum RM(4,9)
([512,256,32]) constituent in all 256 rows and the fixed audited RM2Sub
((t,s)=(64,14)) inner.

### What counts as success

Immediate success is an independent outward checker that authenticates every
input, witness, and occupation and proves

\[
  \Pr[\exists m\ne0:\operatorname{wt}(\operatorname{Enc}(m))\le13107]
  \le 2^{-40}.
\]

The current diagnostic margin is 42.577981756604 bits, so the available
rounding budget is about 2.578 bits. Longer-term success is to use the same
proof interface to produce, for each requested (k), (a) the fastest proven
10%, 40-bit design and (b) the fastest extrapolated design under one narrow,
explicit, justified assumption.

## 2. Current state

### Repository state

- Current `HEAD`: `dbe8b59 Develop block expand constituent proof route`.
- Recent relevant commits include `7672e0d` (finite expander--convolute
  alternative), `c2075cf` (genus-two BA concentration route), and `9409414`
  (fixed EBCH outer with algebraic mixer).
- The worktree is very dirty. In this workstream,
  `FINITE_K20_ATTEMPT_LOG.md` and `MERGE_SUMMARY.md` are modified, while the
  entire `small_k_replay/` directory and the one-stage sparse-EA certificate
  package are presently untracked. Many other repository directories are
  also untracked. Do not reset, clean, or infer reproducibility from `HEAD`
  alone.
- At this snapshot, PID 87196 is running
  `expander_codes\scripts\prime_field_biregular_ec_singleton_certificate.py
  verify tmp\gf136-d24-m5-local-z.json`. This process is outside this
  workstream. Do not interrupt it and do not start a competing benchmark.

### Current small-(k) result: complete diagnostic, not a formal certificate

The authoritative note is
`small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md`. The consolidated audit says:

- status: `PASS` for receipt consistency and coverage;
- occupations covered: exactly (Q=1,\ldots,256);
- (Q=5,\ldots,29): 4,925 compositions checked;
- (Q=30,\ldots,256): 2,857,249 compositions checked;
- combined binary64 margin: 42.57798175660381 bits;
- weakest occupation: (Q=1).

The selected proof by occupation is:

| Occupation | Selected reduction |
|:--|:--|
| (Q=1) | exact-spectrum one-active evaluator |
| (Q=2) | exact fixed-spectrum two-active evaluator |
| (Q=3,4) | positive two-colour RM ladder |
| (Q=5,\ldots,29) | three-group bridge |
| (Q=30,\ldots,256) | refined three-group bridge |

No occupation is missing. The only remaining theorem work is directed
rounding, witness authentication, and theorem packaging.

### Three construction tracks that must remain distinct

1. **Frozen Structured SPIN.** The frozen candidate is
   `Structured SPIN (B=256,t=128,s=19)`. Its end-to-end implementation passes
   correctness, has checksum `0x95c9d722a9539fef`, and has a 21-trial Peach
   median of 10.823872 ms. Its displayed 11% margin of 55.864642 bits is only
   a diagnostic: it uses a modeled outer spectrum, nearest-binary64
   arithmetic, and an outer sampling model that does not match the desired
   one-reused-constituent interface.

2. **One-stage sparse-EA plus RandomStepConv-M22.** This is a complete outward
   finite certificate for a close variant at (k=2^{20}), (N=2^{21}),
   (d_{\min}\ge228590), with total failure margin
   41.3621359294 bits. It is not the frozen construction: its outer is a
   random degree-33 sparse map followed by an accumulator, and its inner is
   RandomStepConv-M22. Its optimized outer is 1.777 times slower than the
   optimized BCH comparator in the isolated Peach experiment.

3. **Fixed RM(4,9) plus RM2Sub-(t64,s14), small (k).** This is the current
   track. It has a complete, audited binary64 occupation proof template at
   10% and 42.578 bits, but no directed-rounding certificate yet. It is not a
   (k=2^{20}) result, not an 11% result, and not the frozen Structured SPIN
   instance.

### Scalable theorems already in the workstream

Two asymptotic reference theorems are stated as proved and have their own
certificate chains:

- `RANDOM_OUTER_RM2SUB_CERTIFICATE.md`: independent random rate-half local
  outers, RM2Sub-((128,19)), relative distance 0.11,
  (B=9\log_2N+O(1)), and (O(N\log N)) outer work.
- `SINGLE_SAMPLED_BA_RM2SUB_D11.md`: one sampled Golay--BA-3 constituent
  reused everywhere, RM2Sub-((128,19)), relative distance 0.11,
  (B=(39/4)\log_2N+O(1)), and (O(N)) ordinary and transposed work.

Neither theorem certifies the frozen ParityFanout outer or supplies the
desired concrete (k=2^{20}) implementation certificate.

## 3. Key conclusions established so far

### Proved or independently audited

- For any setup distribution, if
  \[
    Z_D=|\{x\ne0:\operatorname{wt}(E(x))<D\}|,
  \]
  then
  \(Pr[E\text{ noninjective or }d_{\min}(E)<D]\le\mathbb E Z_D\).
  This is the common finite first-moment interface.
- A factored Structured SPIN interleaver generally preserves more than total
  Hamming weight. The exact proof object is
  \[
    \mu_D=\mathbb E_O\sum_t A_t(O)p_{O,t}(D),
  \]
  not automatically \(\sum_h A_h p_h\). The weight-only formula is legal
  only after proving the required orbit transitivity or a pointwise transfer
  bound.
- The small-(k) three-group identity and its positive
  change-of-measure reduction are exact. Direct enumeration checks the
  three-binomial collapse to (3.7748\times10^{-15}) absolute error, and the
  binary matrix-power audit differs by at most
  (1.0658\times10^{-14}) in log units. These figures audit the binary64
  implementation; they are not outward error bounds.
- The RM(4,9) spectrum is complete and authenticated: it has mass (2^{256}),
  minimum weight 32, and SHA-256
  `995aab561da18f22074b5c6f5413882f492084510aa1cd19b848355f1fcd4ed7`.
- The refined spectrum partition closes every middle and dense occupation.
  The previous dense failure was proof slack, not a demonstrated low-distance
  codeword.
- The one-stage sparse-EA certificate is a genuine complete finite outward
  certificate. Its manifest status is `COMPLETE_CERTIFICATE_CHAIN_AUDIT`.
- The arbitrary-length wrapper is theorem-safe when it starts from a certified
  admissible instance: input restriction preserves distance, and puncturing
  (q) fixed coordinates lowers distance by at most (q).

### Experimentally supported observations

- In the small-(k) RM2Sub parameter scan, (t=64,s=14) dominates the tested
  persistence-matched (t=128,s=13) and (t=256,s=12) choices for
  occupations (Q=3,\ldots,8). Thus (s+\log_2t) is a useful tuning
  coordinate, not a theorem invariant.
- Exact-spectrum RM(4,9) is substantially better than the available smaller
  BCH constituents for the one-active 10% screen. This does not imply that
  RM(4,9) is close to a random length-512 code; its matched random control is
  much stronger.
- The degree-33 one-stage sparse-EA outer costs 9,733 scalar transposed XORs
  per (512\to256) capacity, compared with 5,868 for two optimized BCH
  constituents. Its measured isolated slowdown is consistent with both XOR
  count and greater instruction/register pressure.

### Plausible hypotheses, not facts

- The current small-(k) 2.578-bit slack may be enough for a carefully
  designed outward certificate. It is not large enough to tolerate crude
  rounding or a new global relaxation.
- The small-(k) proof architecture may extrapolate to a family of fastest
  proven designs. That claim still requires separate outward certificates at
  each (k), or a uniform analytic theorem.
- A pressure-aware sparse circuit or a depth-two sparse design might narrow
  the sparse-EA performance gap, but no result presently shows that it meets
  the 11-ms target.

## 4. Current approach

### Probability space for the small-(k) construction

Fix one linear RM(4,9) map

\[
  C:\mathbb F_2^{256}\to\mathbb F_2^{512}
\]

and reuse it in all (L=256) outer rows. The exact spectrum is fixed; the
constituent is not resampled by row. Setup independently samples:

1. a uniform coordinate permutation within each outer row;
2. a uniform permutation in each of the 512 transposed regions; and
3. a nonzero field scalar for each RM2Sub epoch.

The RM2Sub A/B maps are fixed and audited. The epoch length is (t=64), and
the state dimension is (s=14). All failure probabilities refer only to the
route permutations and field scalars just listed.

### Refined change of measure

Partition nonzero local-codeword weights into

\[
  W_0=[1,95],\qquad W_1=[96,416],\qquad W_2=[417,512].
\]

Use Bernoulli references

\[
  (p_0,p_1,p_2)=(0.25,0.5,0.8677722630069483).
\]

For the exact constituent spectrum (A_w), define

\[
  \eta_j=
  \max_{w\in W_j:A_w>0}
  \frac{A_w}{\binom{512}{w}p_j^w(1-p_j)^{512-w}}.
\]

The current binary64 maxima are

| Group | \(\log_2\eta_j\) | Maximizing weight |
|:--|--:|--:|
| low | 119.966741 | 32 |
| central | 260.318572 | 416 |
| high | 104.761160 | 420 |

For occupation (Q), split the active rows into a composition
((q_0,q_1,q_2)), (q_0+q_1+q_2=Q). Their location factor is

\[
  \binom{L}{q_0,q_1,q_2,L-Q}.
\]

Let (F_u(z)) be the positive RM2Sub transfer for a uniform support of (u)
live positions in one region. The reference-region transfer is

\[
R_{q_0,q_1,q_2}(z)=
\sum_{u_0,u_1,u_2}
\left(\prod_{j=0}^2
\binom{q_j}{u_j}p_j^{u_j}(1-p_j)^{q_j-u_j}\right)
F_{u_0+u_1+u_2}(z).
\]

Conditioned on the total live bits, the region permutation makes their union
a uniform support, which is the reason this collapse is valid. The evaluator
raises the region transfer across all 512 regions with a positive binary
log-semiring matrix power, multiplies the pointwise density cost
(\eta_0^{q_0}\eta_1^{q_1}\eta_2^{q_2}), chooses one Chernoff witness per
composition, and sums only nonnegative terms.

### Certificate-compilation invariant

The outward checker must reproduce this exact positive computation. Discovery
may continue to use binary64 to select witnesses, but the checker must:

- parse the fixed decimal (p_2) as an exact rational, not recompute an
  optimizer;
- upper-bound every (eta_j) with directed rounding;
- authenticate the finite witness table and bind each witness to its
  composition;
- outward-evaluate matrix products, powers, multinomial factors, and final
  log-sum-exp;
- cover each (Q=1,\ldots,256) exactly once; and
- report the bad event as weight at most 13107, equivalently
  (d_{\min}\ge13108).

Do not redesign the bands while compiling the certificate unless a verified
outward bound actually consumes the 2.578-bit slack.

## 5. Important files and code locations

### Read first

- `small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md` — authoritative description
  of the current construction, reduction, result, and remaining obligations.
- `small_k_replay/rm2sub_full_occupation_q1_q256_diagnostic.json` — compact
  authoritative binary64 occupation ledger. SHA-256:
  `5f8ba4b40b4e0cba8efbd09d5772d37c0442488f86709566ea9533943931d86e`.
- `small_k_replay/rm2sub_full_occupation_q1_q256_audit.json` — compact `PASS`
  audit. SHA-256:
  `9d33d3695ea86dd68e8b17e07ed9b05329dcd518db391addcea497782d01a269`.
- `FINITE_K20_ATTEMPT_LOG.md` — long authoritative decision log for tried,
  failed, deferred, and closed finite routes. Read targeted sections rather
  than replaying it chronologically.
- `MERGE_SUMMARY.md` — workstream-wide inventory and proved/conditional/
  diagnostic claim ledger. It is modified and contains older sections; prefer
  the current small-(k) note for the active result.

### Current small-(k) implementation

- `small_k_replay/probe_rm2sub_three_group_bridge.py` — per-composition
  evaluator and witness search. It is exploratory binary64 code.
- `small_k_replay/audit_rm2sub_three_group_bridge.py` — independent checks of
  the three-binomial collapse and binary matrix power.
- `small_k_replay/consolidate_rm2sub_refined_band_bridge.py` — verifies the
  raw Q30--256 chunks and emits the compact dense receipt and audit.
- `small_k_replay/consolidate_rm2sub_full_occupation_diagnostic.py` — combines
  the strongest Q1, Q2, Q3--4, Q5--29, and Q30--256 receipts; this script
  itself writes `rm2sub_full_occupation_q1_q256_audit.json`.
- `small_k_replay/rm2sub_refined_band_bridge_q30_q256_diagnostic.json` and
  `small_k_replay/rm2sub_refined_band_bridge_audit.json` — compact dense
  result and coverage audit.
- `small_k_replay/rm2sub_three_group_bridge_q03_q52_probe_d100.json` — source
  for Q5--29.
- `small_k_replay/rm2sub_epoch_geometry_rm49_t64_s14_q1_d100.json` and
  `...q2_d100.json` — exact-spectrum Q1 and Q2 source ledgers.
- `small_k_replay/rm2sub_q_ladder_rm49_t64_s14_d100.json` — Q3--4 source.
- `small_k_replay/SPECTRUM_SOURCES.md` — authentication ledger for every BCH,
  RM, and random-control spectrum used in the comparisons.
- `scripts/rm512_256_spectrum.csv` — authenticated RM(4,9) spectrum source;
  read-only from this workstream.

### Formal framework and scalable results

- `FINITE_LENGTH_FRAMEWORK.md` — exact first-moment and type-orbit theorems.
- `ARBITRARY_LENGTH_WRAPPER.md` — admissibility, padding, shortening, and
  puncturing rules.
- `ASYMPTOTIC_SCALING.md` — parameter schedules and sparse/bulk proof scales.
- `STRUCTURED_SPIN_THEOREM_TARGET.md` — exact theorem interface and frozen
  obligations.
- `RANDOM_OUTER_RM2SUB_CERTIFICATE.md` — proved random-outer 11% reference.
- `SINGLE_SAMPLED_BA_RM2SUB_D11.md` — proved one-sampled Golay--BA-3 11%,
  linear-time asymptotic variant.

### Closed close-variant certificate and performance evidence

- `ONE_STAGE_SPARSE_EA_RANDOMSTEPCONV_CERTIFICATE.md` — theorem and exact
  probability space for the formal finite sparse-EA certificate.
- `ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json` — complete chain manifest;
  SHA-256
  `052db81f0fb2b3eab5987ecfbd9e9e0b911ba5d448d6c9b92653d54c6ef54197`.
- `audit_one_stage_sparse_ea_certificate.py` — independent semantic entry
  point for that chain.
- `ONE_STAGE_SPARSE_EA_VS_BCH_PERFORMANCE.md` — isolated optimized outer
  comparison; it does not benchmark a complete encoder.
- `SPARSE_EA_FREEZE.md` and `SPARSE_EA_FREEZE_MANIFEST.json` — frozen lane
  map and integrity boundary for this alternative.

### Frozen implementation

- `constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/MAIN_CODE_FREEZE.md`
  — authoritative frozen status; do not edit.
- The source manifest is actually nested at
  `.../frozen_source/SOURCE_MANIFEST.json`, not at the construction root.
- `.../PROOF_STATUS.md` — diagnostic proof ledger.
- `.../PERFORMANCE.md` — 10.823872-ms Peach result and correctness description.

## 6. Experiments and evidence

### Refined RM2Sub bands

Setup: fixed RM(4,9), (k=2^{16}), (N=2^{17}), (t=64), (s=14), 10%
bad-weight cutoff, three reference bands. The Q30--256 calculation checked
2,857,249 compositions. Selected margins include:

| (Q) | Binary64 margin (bits) |
|--:|--:|
| 30 | 1228.638 |
| 53 | 2134.716 |
| 128 | 4866.548 |
| 192 | 6297.419 |
| 232 | 7173.271 |
| 248 | 4983.658 |
| 256 | 1611.249 |

The Q30--256 union margin is 1228.638 bits. The full Q1--256 union is
42.577982 bits and is limited by Q1. Therefore more dense-band search has no
value for this instance unless outward rounding unexpectedly fails.

Useful local recheck commands, run one numerical job at a time, are:

```powershell
python workstreams/finite_asymptotic_theory/small_k_replay/audit_rm2sub_three_group_bridge.py
python workstreams/finite_asymptotic_theory/small_k_replay/consolidate_rm2sub_full_occupation_diagnostic.py
```

The second command rewrites the compact receipt and audit. Do not run it
while another certificate job is active. The raw Q30--256 per-composition
chunks were deleted after consolidation because they occupied about 800 MB;
their hashes remain in the compact receipt, and the probe can regenerate
them.

### Parameter calibration

The primary tranche compared exact BCH [128,64,22], exact RM(4,9), and
matched random controls under RM2Sub. At (k=2^{16}), RM(4,9) was the only
exact-spectrum candidate to clear the 40-bit Q1 screen. In the Q ladder,
(t64,s14) gave margins 82.414, 110.168, 135.423, and increasing values for
Q3--Q8; it exceeded the tested persistence-matched larger-epoch variants at
every such occupation. These are exploration results, not formal distance
certificates.

### One-stage sparse-EA certificate

Setup: (k=2^{20}), (N=2^{21}), local 256-to-512 degree-33 sparse map,
one accumulator, 16 bounded full-rank setup attempts, the same accepted local
map reused in 4,096 rows, structured route, and RandomStepConv memory 22.
The cap/setup failure margin is 41.3632966090 bits; the conditional transfer
margin is 51.6422972013 bits; their union is 41.3621359294 bits. The audit
entry point is:

```powershell
python workstreams/finite_asymptotic_theory/audit_one_stage_sparse_ea_certificate.py
```

This is the cleanest existing example of the proof-engine packaging standard.

### Performance evidence

On Peach (Ryzen 9 7950X, CPU 0, GCC 15.2, separate processes, 31 trials), the
one-stage sparse-EA outer averaged 8.703 ms across two medians, versus 4.896
ms for two optimized BCH constituents. The result is isolated outer work,
not end-to-end. The frozen Structured SPIN encoder's independent end-to-end
median is 10.823872 ms on CPU 15. Do not add the two numbers as if they came
from one identical harness; the inference is only that the current sparse-EA
replacement is unlikely to meet the 11-ms target.

## 7. Failed approaches and dead ends

- **Treating the frozen modeled BCH-like spectrum as the actual spectrum.**
  This produced attractive 9--11% margins but is not a proof. Revisit only
  after authenticating the actual constituent or proving a sufficient
  envelope.
- **Resampling a constituent for every row.** This simplifies first moments
  but violates the desired interface. A theorem for independent-row outers is
  a diagnostic, not the requested construction.
- **Pure bit transpose without the region permutations.** Exact obstruction
  calculations show damaging adjacency survives. The project term “bit
  transpose” means the complete factored route, not a literal transpose.
- **Unmodified repeated Golay--BA-3 at finite (k=2^{20}) with RM2Sub-S19.**
  Extensive attempts did not close a complete 11% or 10% finite certificate.
  This is not a distance counterexample; revisit only with a new transfer or
  construction idea. The separate asymptotic BA theorem is nevertheless
  proved.
- **BCH250-124 plus 56 independent row-local fanout layers.** This produced a
  complete 11%, greater-than-40-bit mathematical certificate, but it changed
  the desired one-reused-constituent interface and had unattractive overhead.
  The “DO NOT RUN” implementation was removed. Do not reconstruct it casually.
- **EBCH32--ParityFanout--BA.** It produced a complete 10% mathematical
  checkpoint, but its bounded setup requires an impractical exact test and
  used independent rows. It is not deployable as stated.
- **Repeated EBCH128 plus RM2Sub-S19.** Q1/Q2 diagnostics gave about 26.335
  combined bits at 11%, but dense all-one and mixed-tail classes remained
  open. The old “complete 26-bit certificate” wording was corrected.
- **Random Toeplitz inner.** Repeated Golay--BA plus Toeplitz has a complete
  high-margin finite proof model, but Toeplitz is not the user's random
  bounded-memory convolution and does not establish a linear-time inner.
- **Enumerating all messages of a 128- or 256-dimensional outer.** This is not
  a viable setup test. Use exact known spectra, analytic moments, or
  efficiently checkable spectrum caps; never propose enumerating (2^{128})
  or (2^{256}) words.
- **Old three-group dense bridge with bands 1--47, 48--416, 417--512 and
  references 1/4, 1/2, 3/4.** It failed badly, for example -3665.254 bits at
  Q224. This was repaired by moving weights 48--95 into the low group and
  optimizing the fixed high reference. Do not rerun the old partition as if
  it remained the active approach.
- **Optimizing only the low-band density envelope.** A sparser Bernoulli
  reference improves the pointwise outer envelope but weakens the inner
  low-output bound. The current (p_0=1/4) is deliberate.
- **Concluding that any failed relaxed bound is a low-distance codeword.** The
  finite log records several invalid or overly lossy Hölder, marginal, and
  categorical relaxations. They are failed proof techniques, not construction
  refutations.

## 8. Things that are easy to misunderstand

- **SPIN names are fixed.** Use Accumulator SPIN, Random SPIN, Structured
  SPIN, and Linear SPIN. SPIN means Single-Permutation INterleaved. Do not use
  “Fast SPIN” or lowercase the family name.
- **One reused constituent is a hard requirement for the desired outer.** It
  may be sampled once during setup, but it is not independently resampled at
  every block position. Row-coordinate and region permutations may remain
  random.
- **Expected spectrum is not automatically a high-probability realized
  spectrum.** If the same sampled constituent is reused, the proof needs a
  good-spectrum event, a pointwise envelope, or a joint product bound. Paying
  an average independently in every row is invalid.
- **The outer/transfer product is the object to bound.** Separate average
  bounds on (A_t(O)) and (p_{O,t}) do not bound
  (mathbb E[A_t(O)p_{O,t}]) when they are dependent.
- **RandomStepConv, Toeplitz, and RM2Sub are different inners.** A random
  bounded-memory step convolution uses a fresh random linear state/output map
  at each step over `(state,input)`. Toeplitz is not that distribution.
  RM2Sub uses fixed A/B maps plus sampled nonzero epoch multipliers.
- **“Bit transpose” is shorthand for the full factored route.** The region
  permutations are essential to the current support-uniformity argument.
- **Distance cutoffs have an off-by-one convention.** The framework defines
  (Z_D) using output weight (<D). For the small-(k) 10% target, bad
  weight is at most 13107 and the certified minimum distance would be 13108.
  For frozen (N=2^{21}), bad weight at most 230686 corresponds to
  (D=230687). `PROOF_STATUS.md` labels `D=floor(0.11N)=230686`, while
  `STRUCTURED_SPIN_THEOREM_TARGET.md` uses the theorem-safe 230687 convention;
  normalize this before any formal statement.
- **Frozen implementation units need care.** `MAIN_CODE_FREEZE.md` describes
  (2^{20}) message and (2^{21}) code lengths in 128-bit blocks, whereas
  many theorem documents use binary coordinates. Never transfer a distance
  or rate without stating the unit.
- **The frozen source manifest is nested.** The correct path is
  `.../frozen_source/SOURCE_MANIFEST.json`. A construction-root
  `SOURCE_MANIFEST.json` does not exist.
- **A large displayed margin may still be diagnostic.** “PASS” can mean
  receipt consistency, not outward theorem arithmetic. Check each JSON
  `status` and its prose scope.
- **Full local distance is not the only design metric.** Combined spectrum,
  sparse-class suppression, XOR count, register pressure, and setup
  checkability matter more than maximizing constituent minimum distance in
  isolation.
- **Do not silently change the construction to make a proof close.** The user
  explicitly requires approval before reinterpreting a proposal. State any
  new sampler, per-row randomness, padding, fanout, or inner map as a new
  variant.

## 9. Open questions

1. **Can the 42.578-bit small-(k) diagnostic be converted to a rigorous
   outward certificate without losing 2.578 bits?** The mathematical cover is
   complete. The unknown is numerical compilation and authentication.
2. **What is the cleanest outward representation?** Candidates are exact
   rationals for combinatorial factors plus Arb/MPFI intervals for logarithms
   and matrix powers, or a fully rational positive recurrence with certified
   exponential enclosures. The one-stage sparse-EA checker is the repository's
   best packaging model.
3. **After small-(k) closure, how do the fastest proven parameters scale
   with (k)?** The desired output is a (k\mapsto) best proven design curve
   at 10% and 40 bits, retaining sub-40 results rather than discarding them.
4. **Which narrow extrapolation is acceptable for the “fastest extrapolated”
   curve?** It must be stated explicitly, probably as a calibrated transfer or
   spectrum-envelope assumption, and must never be presented as proved.
5. **Can the frozen Structured SPIN outer receive an actual spectrum or
   transfer-weighted envelope under a one-reused-constituent probability
   space?** This remains the direct route to certifying the 10.823872-ms code.
6. **Can a proof-compatible alternative beat optimized BCH in XORs and live
   state?** Sparse-EA is certified but too slow in its current circuit; no
   pressure-aware or depth-two replacement is closed.
7. **How should the public arbitrary-length interface prioritize exact output
   length versus exact message dimension?** The wrapper proves the available
   tradeoff but does not choose the API.

## 10. Recommended next steps

1. **Freeze the current diagnostic inputs before editing arithmetic.** Create
   a workstream-local manifest containing the hashes of the RM spectrum, fixed
   RM2Sub A/B maps, Q1/Q2/Q3--4/Q5--29/Q30--256 receipts, and this handoff.
   This prevents certificate work from silently changing the candidate.
2. **Write a read-only outward checker for Q1 first.** Q1 is the 42.578-bit
   bottleneck and therefore determines whether the 40-bit theorem is viable.
   Parse the existing witness rather than re-optimize it. If Q1 retains more
   than 40.5 bits, continue; if it falls below 40, improve only Q1 before
   touching dense bands.
3. **Compile one representative refined-band composition outward.** Include
   exact multinomial factors, exact rational band probabilities, outward
   (eta_j), the positive RM2Sub recurrence, matrix powering, and the
   Chernoff coefficient extraction. Compare its interval against the stored
   binary64 value. This isolates arithmetic-design errors before processing
   2.86 million compositions.
4. **Run the complete outward Q30--256 checker serially or with deterministic
   disjoint chunks.** Never run it concurrently with another numerical
   certificate or benchmark. Bind every chunk by hash and make the final
   consolidator verify exact coverage before summing.
5. **Compile Q2--29 and the final union.** Emit one certificate manifest with
   the exact bad event, setup distribution, all hashes, outward margin, and an
   independently rerunnable audit command. Only then change the prose status
   from `BINARY64_DIAGNOSTIC` to a theorem.
6. **Resume the parameter study after closure.** Use the certified engine to
   scan chosen (k), exact-spectrum BCH/RM constituents, and RM2Sub (t,s)
   values at 10%. Produce both the fastest proven and explicitly extrapolated
   frontiers. Record XOR/performance metadata alongside margins.
7. **Return to frozen Structured SPIN only with an authenticated outer
   interface.** The next direct frozen task is not another float sum; it is a
   proof of the actual fixed constituent spectrum or a sufficient
   transfer-weighted envelope under the correct reuse model.

## 11. User preferences and project constraints

- Work only inside `workstreams/finite_asymptotic_theory/` for this workstream.
  Frozen construction sources are read-only evidence.
- Treat performance as a first-class design constraint. The user cares about
  total XORs and concrete Peach performance more than a cosmetically optimal
  distance percentage.
- About 11 ms end-to-end at (k=2^{20}) is the current attractive upper
  limit. Do not restart or replace an optimized implementation without a
  specific reason.
- Never run two benchmarks at the same time. For this project, also avoid
  concurrent long certificate jobs because their timings and outputs become
  hard to interpret.
- Use 10% as the primary distance during exploration. Keep 40 bits as the
  standard comparison line, but report signed margins below it.
- A 10.9% or 11% theorem is welcome when it comes naturally; do not distort a
  practical construction merely to make the number land at 11%.
- The desired outer uses one constituent repeated everywhere. Sampling that
  constituent once with high probability is acceptable. Sampling a different
  constituent per row is not the target.
- Do not use brute-force message enumeration at dimensions 128 or 256.
- Distinguish proof, outward certificate, binary64 diagnostic, benchmark, and
  hypothesis in every result.
- Explore simpler baselines to understand a mechanism, then translate back to
  RM2Sub and the structured route. Do not claim that a random baseline is the
  desired final code.
- Ask before materially reinterpreting a construction. Do not “YOLO” replace
  a permutation by an expander, change the sampler, or add per-row randomness.
- Preserve fixed-width, explicit, low-overhead implementation structure in hot
  paths. Do not trade predictable performance for generic abstractions.
- At the end of each substantial continuation, recommend the next concrete
  turn.

## 12. Astra Bootstrap

- **GOAL:** certify a fast rate-half Structured SPIN or close linear-time
  variant; current practical target is 10% distance and 40 setup-failure bits.
- **CURRENT BEST APPROACH:** compile the fixed-RM(4,9), repeated-constituent,
  structured-route, RM2Sub-(t64,s14), (k=2^{16}) binary64 proof template
  into outward arithmetic.
- **ESTABLISHED FACTS:** exact Q1--256 coverage; audit `PASS`; diagnostic union
  margin 42.577981756604 bits; Q1 is weakest; dense bands have more than 1,200
  bits and need no redesign.
- **CURRENT PARAMETERS:** (k=65536), (N=131072), bad weight (le13107),
  256 repeated outer rows, RM(4,9) [512,256,32], (t=64), (s=14), band
  references (1/4,1/2,0.8677722630069483).
- **IMPORTANT FILES:** `small_k_replay/RM2SUB_REFINED_BAND_BRIDGE.md`,
  `small_k_replay/rm2sub_full_occupation_q1_q256_diagnostic.json`,
  `small_k_replay/consolidate_rm2sub_full_occupation_diagnostic.py`,
  `FINITE_LENGTH_FRAMEWORK.md`, `FINITE_K20_ATTEMPT_LOG.md`,
  `ONE_STAGE_SPARSE_EA_CERTIFICATE_MANIFEST.json`.
- **DO NOT REPEAT:** old 1--47/48--416/417--512 dense bands; independent outer
  resampling; pure transpose; brute-force (2^{128})/(2^{256}) enumeration;
  unmodified finite Golay--BA dense proof search; treating float diagnostics as
  certificates.
- **BIGGEST UNCERTAINTIES:** whether outward rounding preserves the 2.578-bit
  slack; best certificate arithmetic; how the certified parameter frontier
  scales; actual frozen-outer spectrum under one-code reuse.
- **NEXT 3 ACTIONS:** (1) hash-freeze small-(k) inputs, (2) outward-certify Q1,
  (3) outward-certify one refined composition before the full Q30--256 run.
- **WORKTREE WARNING:** critical current artifacts are untracked; do not clean
  the worktree. An unrelated expander-code verifier was running at snapshot
  time.
