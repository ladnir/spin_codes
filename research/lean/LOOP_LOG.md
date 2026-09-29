# Loop log

Append one entry per iteration. Newest last.

## Iteration 0 — setup

- Built the milestone-1 skeleton (see `ROADMAP.md`): 1,013 lines, builds
  clean, no `sorry`, axiom closure `[propext, Classical.choice, Quot.sound]`.
- Added `SpinCodes/Pin.lean`, a kernel-checked statement pin. Negative-tested
  it: adding a field to `ScalableCertificate` breaks `mkPinned`; weakening the
  sparse constant breaks the build.
- Added `scripts/check.sh` (build + sorry/axiom scan + axiom closure). Passes.
- Queue in `TASKS.md`, T1 next.

## Iteration 1 — T1, concrete family

- Added `SpinCodes/Structured/Instantiation.lean`: `Family` (messages pinned to
  `𝔽₂^K`, codewords to `𝔽₂^N`, randomness spaces abstract with bundled
  `Fintype`), and `ScalableCertificate.ofFamily`.
- `probBad_le`, `probBad_nonneg`, `EZ_nonneg` are now **proved**, not assumed:
  `Family.probBad_le` is `Setup.prob_bad_le_cond_sum` applied to the family.
  `Family.distance_whp` concludes about the family's own `probBad`.
- **Design call, flagged.** T1 said "remove `probBad_le` from the certificate".
  Deleting a field breaks `mkPinned`, which is a pin change, which is a stop
  condition. Resolved by keeping the field and supplying it by proof through
  `ofFamily`. Same effect, pin untouched.
- Pin additions (monotone — never relaxes the guard): `probBad` is pinned to
  equal `Pr[∃ nonzero message encoding to weight ≤ d_m]` via `one_le_Z_iff`,
  so it cannot later be redefined into something trivially small; plus
  `probBad_le` and `exists_nonzero_light`.
- **Deferred:** the concrete randomness parametrization (two BA permutations,
  route shuffles, IMT transvections). `prob_bad_le_cond_sum` is uniform in the
  randomness spaces, so nothing is lost; writing down a parametrization before
  the outer and inner maps exist would be guessing. Queued as T1a.
- `scripts/check.sh`: PASS. Next: T2.

## Iteration 2 — T2, Jensen and the native schedule

- `Structured/Majorant.lean`. `â_BA` is modelled the way the certificate
  actually represents it — `ConcaveMajorant` = a `Finset (ℚ × ℚ)` of affine
  supports, `toFun` = their lower envelope. Concavity and
  `eq:structured-ba-jensen` follow from the *representation*, so they are
  unconditional: whatever the 39 rational pairs are, the Jensen step holds.
  What the certificate still owes is `a_BA ≤ â_BA` on `𝒲` (T7) — a different
  statement, not touched here.
- `Structured/Schedule.lean`. `Lsched`, `SchedOK`, `bsched`, `Nsched`.
  `exists_sched` proves a positive multiple of 24 satisfying
  `b ≥ (39/4) log₂(L_m b)` exists — worth having, since without it the
  schedule could be empty and the theorem vacuous. `bsched_tendsto` and
  `Nsched_tendsto` give divergence.
- Existence avoids any numeric bound on `log 2`. The linear comparison would
  need `(39/4)/log 2 < 24`; using `log x ≤ 2√x` instead makes the constraint
  `A√k ≤ 24k`, settled at `k ≥ (A/24)²` for any `A`.
- `ScalableCertificate.ofNativeFamily` discharges `N_eq`, `b_pos`,
  `b_tendsto`, `N_tendsto`, `L_le_N`. With T1, a native-schedule family now
  needs only 6 hypotheses instead of 24 — the selection event and the three
  regime envelopes.
- Index shift recorded: `L_m = 128(m+1)`, since the paper takes `m` positive.
- Pin additions (monotone): concavity, Jensen on `𝒲`, schedule existence,
  `bsched` spec/divergence.
- Split off T2a (zero-extension wrapper) — a corollary, not load-bearing.
- `scripts/check.sh`: PASS. Next: T2a, then T3 (sparse Collatz).

## Iteration 3 — T2a and T3

**T2a.** `Wrapper.lean`: zero extension preserves injectivity and every
codeword weight, hence minimum distance. The companion rate claim is about
*consecutive schedule lengths*, not the extension, so it was split off as T2b.

**T3 — the sparse Collatz certificate is now kernel-checked.**

- Corrected a wrong lead first: `sparse.py` is the *diagnostic* ("not a
  theorem", status `..._NOT_CERTIFICATE`). The real generator is
  `certify_imt_sparse.py`. TASKS.md had pointed at the wrong file.
- `SPARSE_EXACT.json` stores only a SHA-256 of each residual's coefficients,
  not the coefficients. Regenerated all seven from `certify_imt_sparse.py`
  and **all seven SHA-256 digests match the stored values** — so the Lean
  data is the certified polynomial, not a lookalike recomputation.
- All 7 rows use `positive_power_tail` with `zero_order: 1`, degree 257,
  exactly as the appendix describes. (The verifier also has a Bernstein
  fallback; it is unreachable for these rows. Whether the 693 dominance
  comparisons use it is unchecked — noted for T9.)
- Denominator LCM per row is ~1016 digits, effectively a single common
  denominator, so clearing it gives integer coefficients ~1090 digits. Sign
  is preserved (`rat_poly_neg_of_int_scaled`).
- Checked by the kernel via `decide`. Compiled evaluation was avoided on
  purpose — it would enlarge the axiom closure, and `scripts/check.sh` now
  has a check that would catch it. `sparseRows_cert` closure: `[propext]`.
- **Performance note worth keeping.** First attempt had each data module
  import Mathlib; seven parallel builds exhausted memory and the build
  failed. Splitting the definitions into a dependency-free
  `PolyCertDefs.lean` (plain `Int` recursion, no `List.map`) fixed it *and*
  cut the per-row kernel check from 75s to 6.5s.
- `scripts/check.sh` hardened: step 2 also rejects compiled-evaluation
  tactics, and a new step 4 asserts the certificate's axiom closure contains
  no compiled-evaluation axiom.
- **What T3 does not establish**, recorded as T3a: that these polynomials
  really majorize the residuals of `T_occ`. That needs `T_occ` itself (T9).
  The numbers are certified; the bridge is owed.
- `scripts/check.sh`: PASS. Next: T2b, then T4 (exact BA spectrum).

## Iteration 4 — T2b, schedule gap and wrapper rate

- The paper derives "consecutive direct lengths have relative gap o(1)" from
  `b_m = (39/4) log₂ N_m + O(1)`. A cheaper route works and is what I proved:
  `bsched` is **monotone** (a longer outer makes the constraint strictly
  harder, `schedOK_of_succ`) and each step raises it by **at most one multiple
  of 24** (`bsched_succ_le`). Since `b_m → ∞`, the ratio is squeezed between
  `1` and `1 + 24/b_m`. No asymptotic characterization of `b_m` is needed at
  all.
- The step bound is the only real content: `L` and `b` each at most double
  across a step, costing `log₂ 4 = 2`, and `(39/4)·2 = 19.5 ≤ 24`. That the
  schedule constant leaves exactly this much slack is what makes the argument
  go through.
- `mOf n := Nat.findGreatest (fun m => Nsched m ≤ n) n`, with
  `Nsched_mOf_le`, `lt_Nsched_mOf_succ`, `mOf_tendsto`, and
  `wrapper_rate : N_{m(n)}/n → 1`. Together with `zeroExtend_preserves` from
  T2a this is `eq:structured-requested-length-wrapper` in full.
- Protocol note: `lt_Nsched_mOf_succ` was initially stated with a hypothesis
  `Nsched 0 ≤ n` that the proof never used. Dropped it rather than leave an
  unnecessary hypothesis on a lemma.
- Pin additions (monotone): gap limit, `bsched_mono`, `bsched_succ_le`,
  `wrapper_rate`.
- `scripts/check.sh`: PASS. Next: T4 (exact BA spectrum).

## Iteration 5 — T4, exact accumulator transition

- Verified the formula numerically first (exhaustive, b ≤ 13) before
  formalizing. Worth it: it surfaced that the naive `ℕ` translation of
  `T_b(a,c)` is **wrong at c = 0** — truncated subtraction turns `C(-1, r-1)`
  into `C(0, r-1) = 1` when `r = 1`, against the paper's "binomial
  coefficients outside their natural range are zero". `accT` guards `c = 0`
  and `b < c` explicitly, and the docstring records why both are load-bearing.
- Avoided run-length combinatorics entirely. Splitting on the **last** letter
  gives `N_{b+1}(a, c + a mod 2) = N_b(a, c) + N_b(a-1, c)`, because appending
  a letter extends the accumulation by the parity of the whole word
  (`accWtL_append_mod`). The closed form satisfies the same recursion in both
  parities by **one application of Pascal's rule each**. That turned what
  looked like a hard counting argument into a routine induction.
- `accCount_eq_accT` is the result: `T_b(a,c)` counts exactly the words it is
  claimed to count, so `P_b(a,c) = T_b(a,c)/C(b,a)` is the exact
  weight-transition law of a uniformly interleaved accumulator.
- Added semantic regression checks (hand-computed values, `by decide`). These
  pin the *definitions*, which `accCount_eq_accT` alone does not — if
  `accWtL` were wrong, both sides would move together. **One of my
  hand-computed rows was wrong and the check caught it** (`b=4, a=2` is
  `[0,3,2,1,0]`, not `[0,1,2,2,1]`). Exactly the failure mode the checks exist
  for.
- Worked over `List Bool` with words built by appending on the right, so the
  recursion is structural. An earlier attempt over `Fin b → Bool` would have
  needed `Fin.snoc` index wrangling throughout.
- **T4 split**: the accumulator half is done; composing two transitions with
  the Golay enumerator (`eq:structured-exact-ba-spectrum` proper) is queued as
  T4b.
- `scripts/check.sh`: PASS. Next: T4b.

## Iteration 6 — T4b, the uniform interleaver

- `AccTuple.lean`: carried the accumulator count from `List Bool` to
  `Fin b → Bool`, where permutations of positions actually act. The whole
  bridge is `List.ofFn (Fin.snoc v x) = List.ofFn v ++ [x]` — the append
  lemmas transport verbatim and only the cardinality split needed new work.
  The list representation was the right choice for the recursion (structural
  append) and the tuple representation is the right one for permutations;
  neither was wasted.
- `Interleaver.lean` proves the step the paper asserts in a single sentence,
  "the second permutation makes the word uniform on its Hamming slice":
  - `exists_perm_comp` — the permutation action is transitive on each Hamming
    slice (built from `Equiv.sumCompl` on support/complement).
  - `card_fiber_eq` — all fibers of `τ ↦ u ∘ τ` over a slice have equal size,
    by right-translation. This is orbit-stabilizer in concrete form.
  - `sliceCount_eq_choose` — the slice has `C(b,a)` elements, by the same
    snoc split and Pascal.
  - `card_perm_accWt_eq` — **`P_b(a,c) = T_b(a,c)/C(b,a)` is exactly the
    transition law of a uniformly interleaved accumulator**, stated
    multiplicatively to stay in `ℕ`.
- **T4b split further.** What remains of the original T4b: composing two
  independent stages (T4c), the direct-sum enumerator convolution (T4d), and
  the Golay weight enumerator itself (T4e). Recorded T4e as *computational*:
  the extended Golay code has only 2¹² = 4096 codewords, so a concrete
  generator matrix plus a kernel computation should discharge it the same way
  T3 discharged the sparse certificate. That is worth knowing — it was the
  piece I expected to be a blocker.
- `scripts/check.sh`: PASS. Next: T4c.

## Iteration 7 — T4c, two-stage composition

- Added the accumulation as a *map* (`accL` on lists, `accF` on tuples), not
  just its weight: the second permutation acts on the first accumulator's
  output. `wtL (accL v) = accWtL v` falls out term-by-term, since the two
  definitions build the same sum.
- Probability layer: `FinPMF.uniform`, `uniform_prob`, `expect_comp`
  (expectation of a function of a finite-valued variable), and
  `prob_prod_eq_expect` (a joint probability is the expectation over the first
  factor of the conditional probability under the second).
- `prob_two_stage`: **two interleaved accumulator stages compose by the
  Chapman-Kolmogorov sum**, `Σ_c P_b(a,c) P_b(c,w)`.
- Worth recording precisely, because the paper's phrasing is terse and easy to
  misread: the first stage's output is **not** uniform on its Hamming slice
  (the accumulator is a bijection, so it is uniform on the *image* of the
  slice, which is a different set). What makes the composition valid is that
  `τ₂` is independent of `τ₁`, so `prob_accWt` — which holds for *any* fixed
  word of weight `c` — applies to whatever the first stage produced. The paper
  says exactly this ("conditioned on that weight"); the formalization just
  makes it unmissable.
- `sum_Pt`: the transition probabilities out of a fixed weight sum to 1. Not
  needed downstream yet, but it is a real check on the whole `accT` edifice —
  if the closed form or the guards were wrong, this would almost certainly
  fail.
- Split off T4f: summing the per-word law over input words to get `Ā_b(w)`
  itself, which needs T4d first.
- `scripts/check.sh`: PASS. Next: T4d (direct-sum enumerator).

## Iteration 8 — T4d, direct-sum enumerator

- `Enumerator.lean`: `card_tuples_weight` — the number of `k`-tuples of block
  messages with total weight `a` is `[u^a] W(u)^k`, where `W` is the per-block
  weight enumerator. This is the general fact behind the paper's
  `G_b(a) := [u^a] G(u)^{b/24}`, with no Golay and no coding theory in it.
- Stated over an arbitrary finite block type with an arbitrary weight
  function, matching how the paper actually sets it up ("an even base length
  `a` and an injective rate-one-half map `H_a` ... its weight enumerator
  `W_a(z)`"). The Golay instantiation is then a separate, clearly-scoped
  obligation (T4e) rather than being baked in.
- Uses `Polynomial ℕ` so the statement reads as the paper's coefficient
  extraction literally. The induction is `Fin.cons` splitting on the first
  block against `Polynomial.coeff_mul`; the two sides meet after regrouping
  the antidiagonal sum fiberwise by the first block's weight.
- Added a definitional sanity check: for a single bit the enumerator is
  `1 + u`, so the theorem specialises to the binomial count. As with T4, the
  general theorem alone would not catch a wrong `enumerator` definition.
- `scripts/check.sh`: PASS. Next: T4e (Golay weight enumerator, computational).

## Iteration 9 — T4e, the Golay weight enumerator

**The paper's `G(u) = 1 + 759u⁸ + 2576u¹² + 759u¹⁶ + u²⁴` is now computed, not
assumed.** All 4096 codewords are evaluated by the Lean kernel.

- Verified the generator matrix in Python first (systematic `[I₁₂ | B]`, the
  standard quadratic-residue `B`), confirming both the distribution and
  minimum distance 8, before writing any Lean.
- **Measured before committing to an approach**, which changed the design. A
  *trivial* predicate over `Fin 12 → Bool` already takes the kernel 62s,
  because `Fintype.piFinset` reduction is slow; the Golay predicate would have
  been hopeless there. Over the list enumeration `allWords 12` the same work
  takes ~2m20s including 4096 encodings. So the computation lives on the list
  side and `card_filter_ofFn_eq_countP` — a new general bridge — carries it to
  the tuple side where `enumerator` lives.
- Second design point from measurement: the histogram is accumulated in a
  **single pass** (`bumpAt` fold), so the 4096 encodings happen once rather
  than once per weight class. The 25-pass version hit the heartbeat limit.
- `Golay.lean` holds the expensive computation; `GolayEnum.lean` holds the
  assembly, so editing the latter does not re-run the former. Worth doing —
  the assembly took five iterations to get right.
- Bonus: `min_distance` — no nonzero codeword has weight below 8, so the code
  really is `[24,12,8]`. That follows from the computed distribution for free.
- Axiom closure of both `hist_eq` and `enumerator_eq` is clean (no
  compiled-evaluation axiom).
- `scripts/check.sh`: PASS. Next: T4f (sum the two-stage law over input words
  to get `Ā_b(w)`), which now has all its inputs.

## Iteration 10 — T4f, the exact BA spectrum (lemma complete)

**`lem:structured-exact-ba-spectrum` is now fully proved**, Golay-instantiated,
with no ingredient assumed:

    Ā_b(w) = Σ_{a≥1} G_b(a) Σ_c P_b(a,c) P_b(c,w),
    G_b(a) = [u^a] G(u)^{b/24},
    G(u)   = 1 + 759u⁸ + 2576u¹² + 759u¹⁶ + u²⁴.

- `BASpectrum.lean` does the general case: concatenate blocks
  (`concatBlocks`, `wtF_concatBlocks`), apply `prob_two_stage` per message,
  group the nonzero messages by outer weight, and identify the class sizes via
  `card_tuples_weight`.
- The one step with real content is why the outer sum starts at `a = 1`: it
  needs the base map to send only the zero message to weight zero. The paper
  states this as injectivity of `H_a`. Here it is a *hypothesis of the general
  theorem* and a *consequence* in the Golay instantiation — `W_eq_zero_iff`
  falls out of the computed distribution (`A 0 = 1` together with
  `W 0 = 0`), so nothing is assumed about Golay that was not computed.
- `BAGolay.lean` instantiates and rewrites the enumerator by `enumerator_eq`,
  producing the paper's displayed formula with the literal polynomial.
- Both new modules built first try, which is what I would expect now that the
  hard pieces (interleaver, composition, enumerator, Golay) are in place.
- `scripts/check.sh`: PASS.

**Queue state.** T4 and all its parts are done. The remaining queue is T5
(interval-arithmetic decision), T6/T7 (BA tails and the concave majorant),
T8 (route domination — the Poisson-binomial Chernoff bound, flagged from the
start as the first genuinely research-grade step), T9-T12. Next: T5.

## Iteration 11 — T5, the interval-arithmetic decision

A decision task, so the deliverable is `DECISIONS.md` D1 plus evidence.

**Decision: build a purpose-built rational `ln` enclosure; do not depend on
`girving/interval`.**

- Evaluated `girving/interval` properly rather than guessing: Apache-2.0,
  actively maintained (last push 2026-08-10), software floating point (Lean's
  `Float` is untrusted so they wrote their own), `Approx` soundness
  typeclasses, conservative `exp`/`log`/`pow`/`sqrt`/`sincos` and complex
  boxes. It is genuinely the right library for this class of problem.
- The blocker is version: it is pinned to `v4.27.0-rc1` and this project is on
  `v4.34.0`. Adopting means downgrading onto a release candidate seven
  versions back, or porting across seven versions of Mathlib churn.
- Against that, the scope needed here is narrow — `ln` on rationals in a
  bounded range, no `exp`/trig/complex — and **Mathlib already has the hard
  part**: `Real.abs_log_sub_add_sum_range_le` is a proved remainder bound for
  the log series. With argument reduction to `|x| ≤ 1/3` the series converges
  like `3^{-n}`, so the rest is arithmetic on `ℚ`.
- Recorded the condition that would reverse this: if T9 or T11 need `exp`,
  powers or trig enclosures at scale, porting `Interval` becomes the better
  trade. Flagged to revisit at T9.
- `Numeric/LogBounds.lean` is the evidence, not just an assertion:
  `log_mem_Icc`, `logRem_le`, and a worked bound pinning `log(9/10)` to within
  `1e-11` with 11 terms. The certificates need about ten digits (the
  majorant's largest residual is `-1.70e-9`), so the approach reaches the
  required precision with room to spare.
- Also refreshed `ROADMAP.md`, which still described milestone 1.
- `scripts/check.sh`: PASS. Next: T6 (BA tails, sparse half — pure
  inequalities, no boxes).

## Iteration 12 — T6, BA sparse tail (combinatorial core)

- `BATails.lean` proves the paper's chain
  `p_{2ℓ} ≤ C(b,ℓ)C(D,ℓ)/C(b,2ℓ) = C(2ℓ,ℓ)C(D,ℓ)/C(b-ℓ,ℓ)` up to the middle
  equality, entirely in `ℕ` and division-free, so none of it depends on the
  numeric value of `δ_o`. The two ingredients are a hockey-stick sum and
  `Nat.choose_mul`; both are in Mathlib, and `accCount_eq_accT` is what lets
  them apply to `accT`.
- **A lemma I wrote was actually false, and the build caught it.** My first
  hockey-stick statement summed `C(h-1, ℓ-1)` over `h ∈ range (D+1)`. At
  `h = 0` truncated subtraction gives `h - 1 = 0`, contributing a spurious
  `C(0,0) = 1`, so the claim fails at `D = 0, ℓ = 1`. `accT` itself is
  correctly zero there (its `c = 0` guard), so the fix was to split the
  `h = 0` term off rather than weaken anything. Recorded in the file header —
  this is the second time the `c = 0` truncation trap has bitten, and both
  times it was in a *bounding* term rather than in `accT`.
- Also proved `central_le` (`C(2ℓ,ℓ) ≤ 4^ℓ`, the source of the factor of four
  in `c_o = 4δ_o/(1-δ_o)`) and `sum_accT_even_eq_zero` (the paper's "if
  `ℓ > D` the probability is zero").
- Sanity checks against hand-computed values, as usual.
- **T6 split.** The combinatorial core is done; T6b is the numeric step to
  `c_o^ℓ` plus the odd case, T6c is the `ζ`-moment series and `S_b`.
- `scripts/check.sh`: PASS. Next: T6b.

## Iteration 13 — T6b, the numeric tail bound

- `sum_accT_even_le_pow`: **`eq:structured-ba-one-acc-tail` for even input
  weight**, `p_{2ℓ} ≤ (13/28)^ℓ`, stated as
  `Σ_{h≤D} T_b(2ℓ,h) · 28^ℓ ≤ 13^ℓ · C(b,2ℓ)` — cross-multiplied, so the whole
  proof stays in `ℕ` with no rational arithmetic.
- The substance is `choose_mul_pow_le`: `C(D,ℓ)/C(m,ℓ) ≤ (D/m)^ℓ` for `D ≤ m`,
  proved by induction on `ℓ` through `Nat.choose_succ_right_eq`, with the
  per-step inequality `(D-k)·m ≤ D·(m-k)`. That is where the factor-by-factor
  decrease of the binomial ratio lives.
- The hypotheses are exactly the paper's: `ℓ ≤ D` (otherwise the tail is empty,
  already covered by `sum_accT_even_eq_zero`) and `D ≤ δ_o b`. `tail_ratio`
  derives `112D ≤ 13(b-ℓ)` from them by `omega` — the constants `δ_o = 13/125`
  and `c_o = 13/28` enter only there.
- Checked numerically that the bound is not vacuous before trusting the
  algebra: at `b=125, ℓ=2` the certified bound is within a factor of about
  four of the true tail mass. Added the concrete instance as a Lean check.
- My first assembly had a redundant factor in the calc chain (multiplying
  through by `112^ℓ (b-ℓ)^ℓ` when cancelling `C(b-ℓ,ℓ)` alone sufficed); the
  type error caught it and the restructured proof is about half the length.
- Split off T6d for the odd-weight case, which has the same shape with
  `r_v = ℓ+1`.
- `scripts/check.sh`: PASS. Next: T6c (the `ζ`-moment series and `S_b`).

## Iteration 14 — T6c, the negative-binomial bound

- `NegBinom.lean`: `Σ_{h ≥ ℓ} C(h-1,ℓ-1) z^h ≤ (z/(1-z))^ℓ`, the identity the
  paper uses to evaluate the `ζ`-moment.
- **The inequality for partial sums is all that is needed, so no convergence
  argument is required.** The moment is a finite sum over `h ≤ b`; proving the
  identity as a convergent series would have meant `tsum` and summability
  side-conditions for no gain. Nested induction instead: outer on `ℓ` through
  Pascal, inner on the truncation point.
- The bound closes with **no slack**: `z(A_{k-1} + A_k) = A_k` holds exactly
  for `A_k = (z/(1-z))^{k+1}`, so the induction loses nothing.
- Reindexed the sum as `Σ_{j<H} C(j,k) z^{j+1}` (`h = j+1`, `k = ℓ-1`) to
  remove every truncated subtraction. The paper's `h = 0` term is precisely
  where `ℕ` truncation injects a spurious `C(0,0) = 1` — the same trap as in
  T4 and T6, now avoided by construction rather than patched afterwards.
- Sanity check: at `z = 1/2` the bound is exactly `1`, i.e. a negative
  binomial distribution sums to one — so the constant is tight, not an
  artefact of the estimate.
- Split off T6e for the remaining analytic estimates of the sparse half.
- `scripts/check.sh`: PASS. Next: T6d (odd-weight case) then T6e.

## Iteration 15 — T6d, the odd-weight tail case

- `sum_accT_odd_le_pow`: `p_{2ℓ+1} ≤ (13/28)^ℓ`. With T6b this makes
  `eq:structured-ba-one-acc-tail` hold for **every** input weight, which is
  what the downstream argument actually consumes (`p_v ≤ c_o^{⌊v/2⌋}`).
- Structure mirrors the even case with `r_v = ℓ+1`: the tail sum picks up
  `C(D,ℓ+1)` instead of `C(D,ℓ)`, the reindexing is
  `C(b,ℓ)C(b-ℓ,ℓ+1) = C(b,2ℓ+1)C(2ℓ+1,ℓ)`, and Mathlib's
  `Nat.choose_middle_le_pow` supplies `C(2ℓ+1,ℓ) ≤ 4^ℓ` exactly.
- **The bound I proved is stronger than the paper's.** The paper says the odd
  case carries "an additional factor at most `2δ_o/(1-δ_o)`" `= 13/56`; the
  chain here gives `13/112`, a factor of two smaller. Both are `< 1` so the
  same `c_o^{⌊v/2⌋}` conclusion follows either way — recorded because a
  discrepancy in the *other* direction would have been a finding.
- Built first try. Verified numerically before trusting it (ratios 0.06 down
  to 0.003 across `b ∈ {125,1000}`), and separately re-checked that `accT`
  rows sum to `C(b,a)` — a distribution check that would catch a
  mis-transcribed transition.
- `scripts/check.sh`: PASS. Next: T6e, the remaining analytic estimates of the
  sparse half.

## Iteration 16 — T6e, the sparse moment ratio

- `choose_ratio_le`: `C(b,ℓ)·(b-2ℓ+1)^ℓ ≤ 4^ℓ·ℓ^ℓ·C(b,2ℓ)`.
- **The `ζ` cancels.** The paper's step is
  `(C(b,ℓ)/C(b,2ℓ))(ζ/(1-ζ))^ℓ ≤ (κℓ/(b-2ℓ+1))^ℓ` with `κ = 4ζ/(1-ζ)`;
  writing `κ = 4·(ζ/(1-ζ))` makes `(ζ/(1-ζ))^ℓ` appear on both sides. So the
  whole content is the `ζ`-free `ℕ` inequality above, and no real arithmetic
  is needed at all. Noticing this before writing any Lean saved the iteration
  a real-valued detour.
- Supporting lemma `pow_le_pow_mul_choose`: `C(m,ℓ) ≥ (m-ℓ+1)^ℓ/ℓ^ℓ`, a
  three-step chain of existing Mathlib results
  (`pow_sub_le_descFactorial` → `descFactorial_eq_factorial_mul_choose` →
  `factorial_le_pow`). Worth the grep: the descending-factorial route is much
  shorter than bounding the binomial directly.
- Verified numerically across all `b ≤ 60`: no violations, and tight at
  `ℓ = 1` (ratio exactly `1/2`).
- Split the remainder: **T6f** (monotonicity check and `S_b = O(b^{-3/2})`),
  **T6g** (the sparse upper tail).
- `scripts/check.sh`: PASS. Next: T6f.

## Iteration 17 — T6f, the split-and-bound structure

- `sum_split_geom`: a sum dominated by one geometric series below a threshold
  and another above it is bounded by the two geometric tails. This is the
  abstract form of `eq:structured-ba-sparse-sum`'s split at `j = √b`, and it
  is also the shape of the regime sums from T1, so it is now stated once and
  reused.
- `pow_mul_rpow_sqrt_tendsto`: `x^k · q^{√x} → 0` for `0 < q < 1`. This is why
  the *high* half of the split is still `O(b^{-3/2})` even though its base
  (`C_*η³ < 0.657`) is a fixed constant rather than something shrinking with
  `b` — a `√` in the exponent still beats every polynomial. Reduced to
  `u^{2k} e^{-cu} → 0` by substituting `u = √x`.
- Mathlib has no `Tendsto Real.sqrt atTop atTop`, so that is proved inline
  (three lines via `Real.sqrt_sq`).
- **T6f split.** What is left of the `S_b` estimate is T6h: the monotonicity
  of `(κℓ/(b-2ℓ+1))^ℓ` in `ℓ` — the paper checks it by a log-derivative
  argument, which is the first place in this whole formalization that needs
  actual calculus rather than algebra or combinatorics — plus the numeric
  `C_*η³ < 0.657`.
- `scripts/check.sh`: PASS. Next: T6g (sparse upper tail) then T6h.

## Iteration 18 — T6g, the upper tail is a lower tail

- `accWtL_flipHead`: **`Acc(V + e₁) = ¬Acc(V)`**, so
  `accWt(flipHead v) = length v - accWt v`. Adding the first basis vector
  flips the entering parity at every position, hence complements the whole
  accumulation. `accWtL_ge_iff` states the consequence directly: the upper
  tail of `V` is the lower tail of `V + e₁`, which is how
  `lem:structured-ba-tails` reuses `eq:structured-ba-one-acc-tail` instead of
  proving a second bound from scratch.
- Proof is two short inductions: `accLAux_not` (complementing the entering
  parity complements the output) and `wtL_map_not`.
- Verified exhaustively for `b ≤ 12` — both the identity and the tail-swap
  `iff` — before pinning.
- **Build-cost note.** Touching `Accumulator.lean` invalidates
  `Golay.lean`, whose kernel computation takes ~2.5 minutes, so the full
  check ran past ten minutes. Foundational files are now expensive to edit;
  worth putting new low-level lemmas in a leaf module where possible.
- Split off T6i for the counting step (conditioning on the first bit, and the
  two binomial ratio coefficients).
- `scripts/check.sh`: PASS. Next: T6h.

## Iteration 19 — T6h, the BA sparse-tail constants

- **Checked both of the paper's numeric claims at 50 digits before
  formalizing anything.** Both hold, but tightly:
  - log-derivative at `x = 12η`: `-1.2275267…` against the claimed `< -1.22`,
    margin `0.0075`;
  - `C_*η³ = 0.6565638…` against the claimed `< 0.657`, margin `4.4·10⁻⁴`,
    a relative margin of `0.07%`.
  These are precisely the constants worth verifying rather than trusting, and
  they are correct.
- Formalized the second: `Numeric/BAConstants.lean` encloses
  `ζ = √(13/28)` to seven places, derives `8.5543 < κ < 8.5544`, and proves
  `Cstar_eta_cubed_lt`. Only rational arithmetic plus `Real.exp_one_lt_d9` —
  no logarithms, so `LogBounds` is not needed here.
- **My first rational chain was too tight and failed.** I used
  `4κ/(1-8η) < 34.4934` when the true value is `34.493548…`; `nlinarith`
  rejected it. Recomputed the whole chain with exact rationals before
  retrying, and recorded in the docstring that the chain proves `< 0.6565855`
  rather than the true `0.6565638` — the rounding loss is real and worth
  stating, since the total margin is only `4.1·10⁻⁴`.
- Put the module under `Numeric/` (a leaf) so it does not invalidate the
  expensive Golay computation — following the build-cost note from iteration 18.
- Split: **T6j** is the monotonicity claim (needs a log enclosure with
  argument reduction, since the argument `≈0.105` is far from 1);
  **T6k** assembles `S_b = O(b^{-3/2})`.
- `scripts/check.sh`: PASS. Next: T6i.

## Iteration 20 — T6i, the upper-tail counting step

- `card_upper_tail_le`: upper-tail words of weight `v` inject, via
  `V ↦ V + e₁`, into the union of the lower-tail words of weights `v-1` and
  `v+1`. Stated as a count — **no probability and no division**, so the
  paper's two binomial-ratio coefficients do not appear at all at this stage.
- Worked over `Fin (b+1) → Bool` rather than lists, because there flipping the
  first bit is `Function.update v 0 (!(v 0))`, visibly an involution and hence
  injective — which is the only property the counting argument needs. Over
  lists I would have had to show `flipHead` permutes `allWords b`, which needs
  the enumeration to be duplicate-free, a fact I never proved and do not need.
- The three transport lemmas (`ofFn_flipHeadF`, `accWtF_flipHeadF`,
  `wtF_flipHeadF`) carry the list-level identity from T6g across.
- Verified exhaustively: 570 cases over `b ≤ 11`, all `D` and all `v`, no
  violations.
- New leaf module, so the expensive Golay computation was not invalidated —
  the check ran in normal time again, confirming the iteration-18 note.
- Split off T6l for the probability form (dividing by `C(b,v)` and applying
  the `c_o^{⌊v/2⌋}` bounds).
- `scripts/check.sh`: PASS. Next: T6j.

## Iteration 21 — T6j, the monotonicity check

- `logDeriv_lt`: `ln(κx/(1-2x)) + 1 + 2x/(1-2x) < -1.22` at `x = 12η`. This is
  the paper's monotonicity criterion for `(κℓ/(b-2ℓ+1))^ℓ`, and the **first
  use of the `LogBounds` enclosure built in T5** — so the D1 decision (build
  rather than import) is now validated end to end on a real claim.
- The argument of the logarithm is `≈0.1052`, far from `1`, where the log
  series barely converges. **Reducing by a factor of 8** (`ln A = ln(8A) -
  3 ln 2`) moves it to `≈0.8424`, after which **four terms suffice**. Planned
  the numerics before writing Lean: with the true value `-1.2275267` against
  the claimed `-1.22`, the margin is `0.0075`, so only about three decimals of
  the logarithm are needed — the reduction turns an intractable series into a
  trivial one.
- Chain: `κ < 8.5544` (T6h) → `logArg ≤ 0.1053` → `ln(0.1053) ≤ -2.2507` via
  four series terms plus `Real.log_two_gt_d9` → total `< -1.2261`.
- Built first try.
- `scripts/check.sh`: PASS.

**T6 status.** Of the sparse half of `lem:structured-ba-tails`, what remains is
assembly only: T6k (`S_b = O(b^{-3/2})` from the pieces) and T6l (the upper
tail in probability form). Every analytic and numeric ingredient is proved.

## Iteration 22 — T6k, the `S_b` asymptotic

- `sparse_sum_mul_tendsto`: **`b · S_b → 0`**, which is what the downstream
  union actually needs (the paper states `S_b = O(b^{-3/2})`; the consumer
  only uses `b·S_b = o(1)`, so that is what is proved).
- Both halves handled separately: `mul_div_sqrt_tendsto` for the low half
  (`b·C/(b√b) = C/√b → 0`) and `mul_qpow_sqrt_tendsto` for the high half
  (`b·q^{√b} → 0`, from T6f's `pow_mul_rpow_sqrt_tendsto`).
- The hypothesis is stated in exactly the shape `sum_split_geom` produces,
  with `r/(1-r)` absorbed into `2r` — valid once `r ≤ 1/2`, which holds
  eventually. Keeping the interface aligned with the producing lemma avoids a
  glue step later.
- Avoided `rpow` on the low half by writing `b^{-3/2}` as `1/(b√b)`, which
  keeps everything in ordinary division and `Real.sqrt`.
- Extracted `tendsto_sqrt_atTop` from the inline proof in T6f, since it is now
  used twice. Mathlib still has no such lemma.
- `scripts/check.sh`: PASS. Next: T6l, the last piece of T6.

## Iteration 23 — T6l, upper tail against the neighbouring lower tails

- `upper_tail_le_sums`: the upper-tail count at weight `v` is bounded by the
  two `accT` lower-tail sums at `v-1` and `v+1` — the exact form
  `eq:structured-ba-one-acc-tail` bounds.
- **The paper's two coefficients are the same Mathlib lemma.**
  `C(b,v-1)/C(b,v) = v/(b-v+1)` and `C(b,v+1)/C(b,v) = (b-v)/(v+1)` are both
  `Nat.choose_succ_right_eq` read at `k = v-1` and `k = v`. Since each is at
  most `b`, the composite only needs `C(b,v±1) ≤ b·C(b,v)`, which is what
  `choose_pred_le_mul` and `choose_succ_le_mul` give.
- `sum_accT_eq_card` bridges the two representations the tail work has been
  using: `Σ_{h≤D} T_b(u,h)` as a sum, and the same quantity as a cardinality
  over `Fin b → Bool`. Both were needed — the tail bounds are proved about the
  sum, the injection argument about the cardinality.
- Split off T6m for the final numeric substitution. Worth noting for it: `v-1`
  and `v+1` always share a parity, so it is a single case split rather than
  four, but `v = 2` hits `ℓ = 0` where the tail lemmas require `ℓ ≥ 1`.
- `scripts/check.sh`: PASS.

## Iteration 24 — T6m, the upper tail at odd weight

- `upper_tail_odd_bound`: `U_{2m+1}·28^{m+1} ≤ 41·13^m·b·C(b,2m+1)`. This is
  the paper's `b(ζ⁻²+1)ζ^v` with `ζ² = 13/28` made explicit — the `41` is
  `28 + 13`, collecting the two neighbouring contributions.
- Chose the odd-`v` case deliberately: `v = 2m+1` has neighbours `2m` and
  `2m+2`, i.e. `ℓ = m` and `ℓ = m+1`, both `≥ 1`, so the tail lemmas apply
  with no edge case. Even `v` meets `ℓ = 0` at `v = 2`, queued as T6n where
  the bound is trivial anyway (`S_1 = D ≤ b`).
- Built first try.
- **Verification needed a second attempt to be meaningful.** My first
  numerical check brute-forced words for `B ≤ 16` and reported "0 violations"
  — but also 0 cases, because the hypotheses require `D ≥ 2` and hence
  `B ≥ 20`. A vacuous pass looks identical to a real one in the output, so I
  recomputed combinatorially via `accT` and got 259 genuine cases across
  `B ∈ [20, 400]`, all passing. Worth remembering: always print the case
  count, not just the violation count.
- `scripts/check.sh`: PASS.

## Iteration 25 — T6n, the even-weight upper tail (T6 complete)

- `upper_tail_even_bound`: `U_{2(k+1)}·28^{k+1} ≤ 41·13^k·b·C(b,2(k+1))`.
  With T6m this covers **every** input weight, so the sparse half of
  `lem:structured-ba-tails` is complete.
- **Removed the edge case rather than special-casing it.** `v = 2` was going
  to need separate treatment because its neighbour `v-1 = 1` has `ℓ = 0`,
  where the odd tail lemma requires `ℓ ≥ 1`. Instead I proved
  `sum_accT_odd_le_pow'`, a uniform version covering `ℓ = 0` — where the
  claim is trivial (`T_b(1,h) = 1` on `1 ≤ h ≤ b`, so the mass is `≤ D ≤ b =
  C(b,1)`). One extra lemma, and the assembly has no case split at all.
- Writing the weight as `2(k+1)` rather than `2m` also kept both neighbours
  free of truncated subtraction. Between that and the uniform odd lemma, the
  even case ended up no harder than the odd one — both built first try.
- Verified 314 cases across `B ∈ [20,400]`, **including 55 instances of the
  `k = 0` (`v = 2`) case** that motivated the uniform lemma. Reporting the
  case count explicitly, per the iteration-24 lesson.
- `scripts/check.sh`: PASS.

**T6 is complete.** `lem:structured-ba-tails`, sparse half: the one-accumulator
lower tail for all weights, the moment ratio, the negative-binomial bound, the
constants `κ` and `C_*η³`, the monotonicity check, `b·S_b → 0`, and the upper
tail for all weights.

## Iteration 26 — T7, and a finding that changes the plan

Started T7 (the box covers) by building the computable layer the covers need:
`Numeric/RatLog.lean` — `ℚ`-valued `logUpper`/`logLower` with power-of-two
argument reduction and soundness proved against `Real.log`. That part is done
and checked.

**But the architecture it was built for does not work.** The plan was: make
the enclosures computable over `ℚ`, then discharge each 6,747-box cover with a
single `decide`, as T3 did for the sparse certificate. Measured directly:

| expression | `decide` |
|---|---|
| `Nat.gcd 123456 789012 = 12` | reduces |
| `(123456789 : Int) * 987654321 - 1 < …` | reduces |
| `Rat.blt 1 2 = true` | reduces |
| `(1/3 + 1/6 : ℚ) = 1/2` | **stuck** |
| `(2 : ℚ) + 3 = 5` | **stuck** |

The kernel cannot reduce `ℚ` arithmetic. Notably it is *not* `Nat.gcd` that
blocks — that reduces fine — so the cause is `ℚ`'s arithmetic instances
themselves. Integer arithmetic reduces without trouble, which is exactly why
T3 worked: it was scaled to `ℤ` on purpose.

**Recorded as `DECISIONS.md` D2**, with the re-plan: the covers need integer
fixed-point arithmetic with directed rounding. T7 is re-scoped into T7a (that
layer), T7b (dense tail cover) and T7c (majorant cover).

**This partially re-opens D1.** D1 judged `girving/interval` not worth its
version cost because "Mathlib already has the hard part". That was true of the
*mathematical* layer — `LogBounds` and `RatLog` deliver it, and `norm_num`
handled the one-off constants of T6h and T6j comfortably. But the covers need a
*computable* layer, and building one means fixed-point arithmetic with directed
rounding, which is the bulk of what `girving/interval` is — and precisely why
it ships software floating point rather than using `ℚ`. D1's cost estimate for
the covers was wrong; the import-vs-build question should be revisited before
T9 with this on the build side of the ledger.

`RatLog` is kept: it is sound, it is the basis for the fixed-point layer, and
it is usable by `norm_num` for individual constants.

`scripts/check.sh`: PASS.

## Iteration 27 — T7a: the fixed-point layer

`Numeric/FixedDefs.lean` (99 lines, imports nothing beyond `Init`) and
`Numeric/Fixed.lean` (266 lines, soundness). Split for the reason
`Structured/PolyCertDefs.lean` was split: a data-side module that pulls in
Mathlib makes the certificate checks an order of magnitude slower, and at the
sizes here it exhausted memory outright.

`Fix.mk lo hi` stands for the reals `x` with `lo/scale ≤ x ≤ hi/scale`, with
`scale = 10^30` written as a literal so kernel reduction never runs `Nat.pow`.
Operations: `sc`, `ofInt`, `ofFrac`, `add`, `neg`, `sub`, `mul`, `divInt`,
`pow`, and the verdicts `isNeg` / `isNonneg`. Every one rounds outward.

Directed rounding comes free from `Int`'s division being Euclidean: for a
positive divisor it *is* floor division, so `fdiv n d = n / d` and
`cdiv n d = -((-n) / d)`. The one-line proof is `Int.mul_ediv_add_emod`
plus `Int.emod_nonneg`.

Interval multiplication is the hull of the four corner products, rounded back
down from scale `scale^2`. The corner lemma needs a two-level sign split
(`mul_le_of_corners`); the lower bound then comes from the upper one by
negating the second factor, which avoids writing the whole case analysis twice.

**`Fix.Mem` is pinned by `Iff.rfl`, not just by existence.** This is the one
pin in the file that really matters. Weaken `Mem` to `True` and every
soundness lemma still typechecks, every box cover built on them still
compiles, and every one of them is vacuous — with a green build and a clean
axiom closure the whole way. The `Iff.rfl` pin is what makes that a build
error.

### The measurement D2 asked for

`decide` reduces all of it. Confirmed directly: `ofFrac 1 3` evaluates to
`⟨333…333, 333…334⟩` (note the outward rounding), `pow (ofFrac 1 3) 20`
reduces, and a 2,000-box list check reduces. So the D2 re-plan works.

Throughput, on a list of 8,000 boxes each doing 5 multiplications and 3
subtractions: **2m49s**, about 21 ms per box.

That number reframes T7b and T7c. At 8 operations a box the covers would cost
a few minutes each, which is fine. But a log enclosure is more like 75
operations (≈25 series terms, each a `pow`, a `divInt` and an `add`), and at
that rate 6,747 boxes is over half an hour — per cover, on every build.

The way out is that the log *arguments* repeat heavily across boxes: they are
grid points, not 6,747 distinct reals. So T7b/T7c should tabulate the
distinct log enclosures once and make each box a ~10-operation lookup and
comparison — the same shape as the T3 sparse rows, which are data modules
checked once. Recorded as the note on T7a2.

`scripts/check.sh`: PASS. Axiom closure of `distance_whp` unchanged at
`[propext, Classical.choice, Quot.sound]`.

## Iteration 28 — T7a2: the logarithm in fixed point

`Numeric/FixedLog.lean`, plus the computable half appended to `FixedDefs`.
`flog p q k n` encloses `Real.log (p/q)`, and `flog_mem` proves it against
`LogBounds.log_mem_Icc`.

Two choices worth recording.

**Horner, not term by term.** The obvious definition mirrors
`logSeries n y = ∑ y^(i+1)/(i+1)` and so recomputes `y ^ i` at every term —
`O(n^2)` multiplications. At the forty-odd terms needed that is seconds per
argument. Horner is `O(n)`, and the real-side mirror `hornerR` plus
`hornerR_eq_sum` identifies it with `logSeries` in one induction, so the
soundness proof costs nothing for the speedup.

**The window is `[3/4, 3/2]`, so `|1 - z| ≤ 1/2`.** A narrower window would
converge faster, but a power-of-two shift can only guarantee a window of
multiplicative width `2` — `[3/4, 5/4]` has ratio `5/3 < 2` and so does not
cover every argument. With `|t| ≤ 1/2` the remainder is `≤ 2^{-n}`, which is
why `n` sits in the forties rather than the twenties. Using `1/(1-|t|) ≤ 2`
also keeps the remainder bound `2·rho^(n+1)`, computable by the same `Fix`
operations, with no interval reciprocal needed anywhere.

`log 2` comes from the same series at `y = 1/2` with 64 terms. It has to:
`Real.log_two_gt_d9` gives only nine digits, and the covers need better than
`1.7e-9`, so the Mathlib constant would have dominated the error at `k ≥ 1`.

### Checks against the truth

    log2          ∈ [0.693147180559945309362, 0.693147180559945309470]
    flog 9 10 0 48 ∈ [-0.105360515657829854, -0.105360515657822749]

against `0.69314718055994530942` and `-0.10536051565782630123`. Both contain
the true value, width `7.1e-15`, which is the predicted `2·2^-48`.

The out-of-window case is genuinely unsound, as it must be: `flog 1 10 4 48`
shifts `0.1` to `1.6`, outside `[3/4, 3/2]`, and returns
`[-2.30258509299422, -2.30258509299421]`, which does **not** contain
`log 0.1 = -2.302585092994046`. At `k = 3` (shift to `0.8`) the enclosure is
correct. So the window hypotheses in `flog_mem` are doing real work, not
bookkeeping — they are pinned along with the conclusion.

### Cost

**~280 ms per log** at `n = 48` (100 distinct arguments, 31.6 s). Slower than
the O(n) count suggested, but fine: `n` is the caller's parameter and `n = 36`
already gives `1.5e-11`, comfortably inside the `1e-9` the certificates need.

This confirms the tabulation plan from iteration 27 rather than displacing it.
Recomputing a log per box would cost 6,747 x 0.28 s = half an hour per cover;
tabulating the distinct grid arguments once brings it to under a minute. Noted
on T7b and T7c, which are now both unblocked.

`scripts/check.sh`: PASS.

## Iteration 29 — T7b, split; T7b1 done

Started T7b and found it is three tasks, not one.

### What the certificate actually is

`eq:structured-ba-dense-tail` is checked by
`certificates/single_sampled_ba_rm2sub/linear_time_audit/verify_outer_interval.py`
in `--mode low-tail`. It is a *live* branch-and-bound at 60-digit mpmath
interval precision: it reports 6,747 processed boxes and no unresolved ones,
but it writes no box list. So there is nothing for Lean to replay yet — the
script has to be extended to dump its leaves, plus the `r` witness it picks
for `g` on each. That is T7b3.

Two things about the script matter for the Lean design.

**The paper's `π` and the script's `p` are the same function**, which I
checked by expanding

    pi(a,c) = c h(a/2c) + (1-c) h(a/(2(1-c))) - h(a)

The `ln c` and `ln(1-c)` terms recombine, the `ln a` terms cancel exactly, and
what is left is

    a ln 2 + c ln c + (1-c) ln(1-c)
      - (c - a/2) ln(c - a/2) - (1 - c - a/2) ln(1 - c - a/2) + (1-a) ln(1-a)

which is `p_point` verbatim. Worth having checked: the script's form is the
one that can be evaluated without dividing by `c`, and the cancellation of
`ln a` is exactly what makes it stable near `a = 0`.

**The script does not rely on the naive interval extension.** It takes the
`min` of the range bound and a mean-value ("Taylor") bound that keeps the
derivative cancellations, and the same again at the level of the whole
objective. On boxes near the feasibility boundary the naive bound is very
likely too loose, so Lean will probably need the mean-value form too. That is
not a problem for the arithmetic: the derivatives are
`ln(2 sqrt(xy)/z) = ln 2 + (ln x + ln y)/2 - ln z` and
`ln(by/((1-b)x))`, so they need logs, products and quotients but no square
root.

### T7b1 — the missing primitives

`inv`, `div`, `flogI`, with soundness and pins.

`div a b = mul a (inv b)`. Writing division as multiplication by the
reciprocal means the four-corner sign analysis is not repeated — `mul_mem`
does it once. `inv` needs only `0 < lo`: then `1/y` runs over
`[scale/hi, scale/lo]`, whose numerators at scale `scale` are `scale^2 / hi`
and `scale^2 / lo`, rounded outward.

`flogI` is two `flog` calls, not a series over intervals: `Real.log` is
monotone, so the enclosure on `[lo, hi]` is pinned by the endpoints, and each
endpoint is the *rational* `lo / scale` that `flog` already takes. Two logs per
interval logarithm instead of one interval series.

Checked numerically: `div (1/3) (1/7)` gives `2.33333333333333333333333333331…`
to `…341` (width 2.4e-29), `div 1 (1/3)` gives `3` to within 1e-29, and
`flogI` on `[0.8, 0.9]` returns exactly
`[log 0.8, log 0.9] = [-0.223143551314213, -0.105360515657823]`.

`scripts/check.sh`: PASS.

## Iteration 30 — T7b2: the BA exponent on `ℝ`

`Numeric/BAExponent.lean`. The real-side definitions exactly as
`structured_appendix.tex` writes them, plus the two facts the interval
verifier relies on without stating.

### The `π` identity

The paper defines

    π(a,c) = c·h(a/(2c)) + (1-c)·h(a/(2(1-c))) - h(a)

and the verifier evaluates something else entirely: a division-free expression
with the `log a` terms already cancelled. `piBA_eq_piEval` proves they agree
for `0 ≤ a < 1` and `a/2 < c < 1 - a/2`.

The cancellation is the whole point of the verifier's form. Expanding
`c·h(a/(2c))` splits `log(a/(2c))` into `log a - log 2 - log c`, and the two
halves contribute `-(a/2)log a` each while `-h(a)` contributes `+a log a` — so
`log a` vanishes exactly. Without that, the expression has a `log a` blowing
up as `a → 0` cancelling against another one, which is precisely the corner of
the box cover where precision matters. The `a = 0` case needs its own branch
(both sides are zero, but the derivation divides by `a`).

Cross-checked outside Lean on 2,005 cases — the six corners of the certificate
region plus 2,000 random feasible points — with maximum discrepancy
`3.7e-16`, i.e. rounding.

### `g` is an infimum and needs a floor

`gBA a ≤ gObj u a` is the step that lets one witness `u` per box discharge
`g`, and `csInf_le` needs `BddBelow`. The family is **not** bounded below for
`a > 1`: `gObj u a ~ (1-a) log u → -∞`. So `gBA_le` carries `0 ≤ a ≤ 1`, which
is not padding — it is the exact condition.

On `[0,1]` the floor is clean and needs no analysis: `G(u) ≥ 1` handles
`u ≤ 1` (then `-a log u ≥ 0`) and `G(u) ≥ u^24` handles `u ≥ 1` (then
`gObj ≥ (1-a) log u ≥ 0`). So `gObj u a ≥ 0` throughout, and `gBA a ≥ 0` too.

### The Golay constants are not retyped

`GolayG_eq_sum` proves

    GolayG u = ∑_{m : Fin 12 → Bool} u ^ W m

by evaluating `Golay.enumerator_eq`, which was itself checked against the
actual generator matrix in `Structured/Golay.lean`. So the `759, 2576, 759`
the verifier hard-codes are now derived from the code rather than asserted
alongside it.

All five definitions are pinned by `rfl` against the paper's formulas. A
drifted definition here would have every box in the cover check the wrong
thing with a green build, so `rfl` pins are the right guard.

`scripts/check.sh`: PASS.

## Iteration 31 — T7b2b: the `Fix` evaluators, and a boundary finding

`Numeric/BAEvalDefs.lean` (evaluators) and `Numeric/BAEval.lean` (soundness):
`flogKJ`, `flogQ`, `flogIW`, `fhEnt`, `fpiEval`, `fGolayG`, `fgObj`.

### Shifts: searched, then checked

Threading a reduction shift through `h`, `π` and `g` as a parameter would mean
a dozen shift arguments per box. Instead `shiftPair` searches for one, and the
result is wrapped in `Option`: the window conditions are *checked by the
definition*, not assumed.

Nothing proves the search succeeds and nothing needs to. If it finds a good
pair the guard passes and the enclosure is sound; if it does not, the answer is
`none` and the box check fails visibly. So `flogQ_mem` and everything built on
it carry no hypothesis about the search at all — which is exactly why this
shape was worth the `Option` plumbing.

### Two things testing caught

**The shift search was one-sided.** `fgObj` returned `none` on its first test.
The reduction window is `[3/4, 3/2]` and my search only scaled *up*, but
`G(u) ≥ 1` always, so the Golay argument needs scaling *down*. Fixed with a
two-sided `shiftPair` and `log(p/q) = log((p·2^k)/(q·2^j)) - k·log2 + j·log2`.
After the fix, `fGolayG (ofInt 1)` is exactly `4096` and
`fgObj 1 0.5 = [0.34657359027948, 0.34657359028046]` around
`log 4096 / 24 = 0.34657359027997`.

**`fpiEval` returns `none` on the feasibility boundary.** At `c = a/2` the term
`(c - a/2)·log(c - a/2)` is `0·log 0`, which the paper and the Python verifier
both read as `0` by continuity — but `flogIW` correctly refuses, because
`log 0` is not what the enclosure machinery can produce. Tested at
`π(0.1, 0.05)`: `none`, where Python gives `-0.1292005`.

This is a real gap, not a formalization artifact: the cover certainly contains
boxes touching `c = a/2`, so the term has to be handled as `t·log t` rather
than as a product of a factor and a logarithm. The verifier does exactly that —
its `xlogx_range` special-cases the zero endpoint and uses the minimum at
`1/e`. Recorded as T7b2c, which now blocks T7b3.

Mathlib has both pieces already: `Real.convexOn_mul_log` on `Ici 0` gives the
upper bound through `ConvexOn.le_max_of_mem_Icc`, and the global minimum
`t log t ≥ -1/e` follows from `Real.log_le_sub_one_of_pos` applied at `1/(e t)`.

### Checked values

    fGolayG 1        = 4096 exactly
    fgObj 0.6 0.05   ∋ 0.14945743750550097
    fpiEval 0.5 0.4  ∋ -0.021005925701837
    fhEnt 0.25       ∋ 0.5623351446188083
    fpiEval 0.008 0.5 ∋ 0 (exact; the true value is 0, Python's 1.1e-16 is noise)

`scripts/check.sh`: PASS.

## Iteration 32 — T7b2c: `x log x` as a primitive

`xlogx` on the real side (`BAExponent`), `xlPoint` and `fxlogx` on the `Fix`
side, and `fhEnt` / `fpiEval` rebuilt on top. This clears the boundary failure
from iteration 31.

**Upper end** is the maximum principle: `t log t` is convex, so on a box its
largest value is at an endpoint. Mathlib supplies both halves —
`Real.convexOn_mul_log` on `Ici 0` and `ConvexOn.le_max_of_mem_Icc`. The
endpoint value at `0` is exactly `0`, which is the whole reason the term has to
be a primitive rather than a factor times a logarithm.

**Lower end** is the tangent at the box midpoint:

    x·(log m + 1) - m  ≤  x log x      (x ≥ 0, m > 0)

This is Gibbs' inequality `x log(x/m) ≥ x - m` rearranged, which in turn is
`log y ≥ 1 - 1/y` — so it needs only `Real.log_le_sub_one_of_pos`, no
differentiation and no constant for the minimum at `1/e`. It is exact at
`x = m`, so it sharpens quadratically as the branch-and-bound subdivides,
which is the right behaviour for a cover.

### Two bugs found by evaluating

**The degenerate box.** `fxlogx` on `{0}` took the midpoint `0` and then asked
for `log 0`. That box is not exotic — it is exactly the feasibility boundary
`c = a/2`. Now its own branch: `0 ≤ lo ≤ hi = 0` forces the value to `0`.

**Dependency.** The tangent was first written as `x·log m + x - m`, which
evaluates `x` twice as independent intervals and loses badly. Factoring to
`x·(log m + 1) - m` uses `x` once. Measured on three boxes, lower-bound error
against the true minimum:

    box          before     after
    [0.1, 0.2]   0.108      0.0076
    [0, 0.2]     0.239      0.039
    [0.3, 0.5]   0.190      0.0070

14x to 27x tighter, from writing the same expression differently. Worth
remembering for the mean-value bound in T7b3, which has the same hazard.

### End to end

At `α = 0.05, β = 0.10, x = 0.05`, witness `u = 0.6` — a point *on* the
boundary, since `x = β/2`:

    g + π(α,β) + π(β,x)  ∈  [-0.0077873258664, -0.0077873258558]

against the Python verifier's `-0.0077873258611`. Width `1.05e-11`, comfortably
inside the `7.68e-8` the dense tail needs.

`scripts/check.sh`: PASS.

## Iteration 33 — T7b3 scoped by measurement

No Lean changes this iteration. Four measurements, each of which would have
been an expensive thing to get wrong by guessing.

### 1. The verifier reproduces exactly

Driving `verify_outer_interval.py` in `low-tail` mode: **6,747 processed,
6,748 accepted leaves, `max_accepted_upper = -7.686951035770771e-08`, depth
26** — matching `artifact/PAPER_MAP.md` to the last digit. So the target is
well defined and the script is a usable oracle.

### 2. The naive range bound does not converge

I disabled both mean-value bounds (`p_taylor_upper` and
`objective_taylor_upper`) and re-ran to a wall-clock limit. It reached
**144,000 boxes processed with 58,920 still queued** and a best remaining
upper bound of `0.00117` — against `6,747` boxes and termination for the full
bound. Twenty-one times the work and still not close; the queue was still
near its peak.

So Lean cannot get away with the naive interval extension: the mean-value form
that keeps the derivative cancellations is load-bearing, not an optimization.
That is T7b3c. Its gradients are `ln 2 + (ln x + ln y)/2 - ln z` and
`ln b + ln y - ln(1-b) - ln x` — logs, products and quotients, no square root,
so the existing arithmetic covers it.

### 3. The cost model

Measured `flogQ` with the value actually forced (an earlier attempt used
`Option.isSome`, which the kernel satisfies without evaluating the payload —
a vacuous measurement, and the second time this project has produced one):

    n = 12   58 ms      n = 24  104 ms      n = 48  245 ms

so roughly `12 + 4.9·n` ms — the series dominates.

A box needs about 30 point logarithms (5 for `p` at the midpoint, 8 for the
gradients, per `π`; two `π`s; plus `G(u)` and `u`). At 245 ms that is 7 s a
box and **19 hours** for the cover. Not viable as it stands.

### 4. Where the factor comes from

**11,210 distinct log arguments across all 6,748 leaves** — about 18-fold
sharing against the ~200,000 requests. Their denominators are dyadic with at
most 19 bits (the bisection grid), except 176 arguments with 46–57 bit
denominators, which are the `choose_r` witnesses. At `scale = 10^30` every one
of those is *exactly* representable, since `2^19 | 10^19` and `1/125` and
`13/125` are exact decimals — so the boxes can be given as exact `Fix` points
with no rounding at the data boundary.

Two changes then bring the cover into range:

* **Odd series** (T7b3a). `log x = 2·atanh z` at `z = (P-Q)/(P+Q)`. The
  reduction window widens to `[2/3, 3/2]` — still reachable by a power-of-two
  shift, since its ratio `9/4` exceeds 2 — and there `|z| ≤ 1/5`, so each term
  buys a factor of 25 rather than 2. Eight multiplications replace forty-eight.
  No new analysis is needed: `2·∑_{k<m} z^(2k+1)/(2k+1)` is exactly
  `logSeries (2m) z - logSeries (2m) (-z)`, so the existing
  `log_one_sub_mem_Icc` applied at `z` and at `-z` gives the enclosure, with
  the two remainders adding. `log 2` itself comes from the same series at
  `z = 1/3`, which is not circular.
* **Tabulation** (T7b3b). 11,210 entries verified once, then looked up. The
  lookup must be `O(log n)` — a flat list at `O(n)` would cost more than
  recomputing — so a binary search tree keyed by the numerator.

Estimated after both: table `11,210 x ~55 ms ≈ 10 min`, boxes
`6,748 x ~55 ms ≈ 6 min`. Around a quarter of an hour of kernel time per
build. Heavy but tolerable — `Golay.lean` already costs 2.5 minutes.

### 5. The covering obligation

Worth recording before it becomes a surprise: a flat list of 6,748 boxes
would leave "these tile the initial box" to be proved, which is real
combinatorial work. Emitting the **split tree** instead makes it structural —
each node is a bisection, and `Icc lo hi = Icc lo mid ∪ Icc mid hi` at each
step. That is T7b3d, and it costs nothing extra on the generator side because
the tree is what the branch-and-bound already builds.

`scripts/check.sh`: PASS (tree unchanged).

## Iteration 34 — T7b3a: the odd series

`log x = 2·atanh z` at `z = (P-Q)/(P+Q)`, replacing the `log(1-x)` series
everywhere the box path uses a logarithm.

**No new analysis was needed.** `two_oddSeriesR` proves

    2 · ∑_{k<m} z^(2k+1)/(2k+1)  =  logSeries (2m) z - logSeries (2m) (-z)

by induction — the even terms cancel because `(-z)^(2k+2) = z^(2k+2)`. So the
enclosure is just `log_one_sub_mem_Icc` applied at `z` and at `-z`, with the
two remainders adding, and `logRem` is unchanged under `z ↦ -z` because it
depends only on `|z|`. Everything reduces to what `LogBounds` already had.

**The window widens to `[2/3, 3/2]`**, whose ratio `9/4` still exceeds 2, so a
power-of-two shift always reaches it. There `|z| ≤ 1/5` — the two window
inequalities `2Q ≤ 3P` and `2P ≤ 3Q` are *exactly* `5(Q-P) ≤ P+Q` and
`5(P-Q) ≤ P+Q`, so the bound is tight rather than incidental. Each term then
buys a factor of 25 instead of 2.

`log 2` cannot go through this path — `2/1` is outside the window and reducing
it would need `log 2`. It comes from the same series at `z = 1/3`, since
`(1+1/3)/(1-1/3) = 2` exactly, with 21 terms.

### Measured

    per logarithm    245 ms  (old, n = 48)  ->  42 ms  (new, m = 9)

5.8x, and the results are *tighter*, not a precision trade. The end-to-end
objective at the same boundary point as iteration 32:

    old (n = 48):  [-0.0077873258664, -0.0077873258558]   width 1.1e-11
    new (m = 10):  [-0.0077873258611431, -0.0077873258610702]  width 7.3e-14

against Python's `-0.0077873258611`. A hundred and fifty times narrower for a
sixth of the work.

That brings a box from ~7 s to ~1.3 s. With T7b3b's tabulation (11,210
distinct arguments rather than ~200,000 requests) the cover should land around
eight minutes of table plus six of boxes.

`scripts/check.sh`: PASS.

## Iteration 35 — T7b3b: the lookup table

### The evaluators now take their logarithms from an oracle

`LogFn := Int → Option Fix`, with `Oracle lg` saying an answer is positive at
the numerator and encloses `log (p/scale)`. Every evaluator is generalized
over it, and the old `n`-indexed names are kept as the instance at
`directLog n` — so the existing pins continue to name exactly what they named
before, and there is one soundness proof rather than two.

### Soundness of the table needs no search-tree invariant

This is the part worth recording. `LogTree.find` returns a value only at a
node whose key it has *compared equal*, so if `check` has verified every node
then whatever `find` returns is a verified entry for that exact key. A
mis-ordered tree can make a lookup fail; it cannot make one lie. The ordering
exists purely for speed and never appears in `find_check`.

Tested both ways: a table built from `flogQ` gives `check = true`, and
perturbing one entry to `⟨0,0⟩` gives `check = false`.

### Measured

    tree elaboration (8,191 literal nodes, 456 KB)   ~45 s, once
    marginal cost of one lookup                       2.3 ms
    cost of computing one logarithm                    42 ms

so a lookup is **18 times cheaper** than a computation — which happens to
match the 18-fold sharing found in iteration 33, and is the whole reason the
table is worth its weight.

Projected for the cover: verifying 11,210 entries at 42 ms is about 8 minutes,
and 202,000 lookups at 2.3 ms is about another 8. Against 142 minutes if every
request recomputed. Add the per-box interval arithmetic and the cover should
land near twenty minutes of kernel time — paid once, since a `.olean` caches
it until the data changes.

`Oracle` is pinned by `rfl`, for the reason `Fix.Mem` is: weakened to `True`
it would leave every evaluator lemma typechecking and saying nothing.

`scripts/check.sh`: PASS.

## Iteration 36 — T7b3c: which bound is actually needed, and the derivatives

### Two things ruled out, one ruled in

Before building the multivariate mean-value machinery I checked whether
something cheaper would do.

**A one-dimensional "difference form".** The two large terms of `p` have the
shape `f(t) - f(t - a/2)` with `f = t log t`, and a *1-D* mean value theorem
turns each into `(a/2)(ln ξ + 1)` for some `ξ ∈ [t - a/2, t]` — capturing the
same cancellation with a fraction of the machinery. It is sound (4,000 random
boxes, no violations) and it does help. It is not enough: at 12,000 boxes the
queue was still 11,062 and the best remaining bound 0.0194.

**Which of the verifier's two mean-value bounds carries the weight.** It runs
one at the level of `p` and another at the level of the whole objective. With
the `p`-level one alone:

    processed 13,512   accepted 13,513   remaining 0
    max_accepted_upper  -7.686951035770771e-08   depth 28

i.e. it terminates, with exactly the endpoint the paper reports, at twice the
box count. So the **three-variable centered form is not needed** — only the
two-variable one for `p`. That is a large simplification for a cost of 2x in
boxes, and the measurement took ten minutes rather than the days the
three-variable version would have cost to build and then discover unnecessary.

### The derivatives

    ∂π/∂a = ln 2 + (ln(c - a/2) + ln(1 - c - a/2))/2 - ln(1 - a)
    ∂π/∂c = ln c - ln(1-c) - ln(c - a/2) + ln(1 - c - a/2)

proved as `HasDerivAt` from `Real.hasDerivAt_mul_log` composed with the affine
inner maps. The `+1` terms from each `t log t` cancel exactly — `1/2 + 1/2 - 1`
in the first, `1 - 1 - 1 + 1` in the second — which is why both derivatives are
pure logarithms and why the verifier can write the first as `ln(2√(xy)/z)`.

Checked against the verifier's own `p_gradient` (agreement 1e-15) and against
central differences (1e-9, which is the difference step). They are pinned as
`HasDerivAt` statements rather than as definitions: a wrong formula here would
make the centered bound *unsound*, not merely loose.

Two Lean notes worth keeping. `HasDerivAt.comp` cannot infer the outer function
from a beta-reduced goal, so the composition has to be supplied explicitly and
unfolded with `Function.comp_def`; going through `simpa [Function.comp]`
instead produces a term at a different instance path and fails to typecheck.
And `HasDerivAt.add` builds `f + g` in the Pi sense, not a lambda, so the
assembly ends with `convert ... using 1` and `funext`, not `rw`.

`scripts/check.sh`: PASS.

## Iteration 37 — T7b3c2: the centered bound

`piEval_centered` on the real side, `fpiCenteredL` and `fpiCenteredL_mem` on
the `Fix` side.

### Two 1-D mean value theorems, not a multivariate one

    π(a,c) - π(am,cm) = [π(a,c) - π(am,c)] + [π(am,c) - π(am,cm)]
                      = ∂_a π(ξ,c)·(a-am) + ∂_c π(am,η)·(c-cm)

Hold `c` fixed and move `a`; then hold `am` fixed and move `c`. Both `ξ` and
`η` land inside the box, which is all the enclosure needs, and Mathlib's
`exists_hasDerivAt_eq_slope` does each step. The helper `mvt_between` absorbs
the three-way case split on which endpoint is larger, so the main proof never
sees it — the degenerate case `u = v` is where the derivative is irrelevant
and any point of the box will do.

The five hypotheses are exactly the verifier's feasibility guard, and on the
`Fix` side they are five `Int` comparisons a certificate discharges by
`decide` alongside the bound itself.

### It does what it was needed for

On `[0.3, 0.31] × [0.4, 0.41]`, where the true maximum of `π` is `-0.007116`:

    naive interval extension   upper = +0.005664   (positive: box unusable)
    centered bound             upper = -0.006576   (within 5.4e-4 of the truth)

That is the whole point. The naive bound cannot certify this box at any depth
of subdivision that matters, and the centered one certifies it immediately.
On the wider `[0.3, 0.4] × [0.4, 0.5]` the two are `0.1206` and `0.0557`, so
the gain grows as boxes shrink — which is the right direction for a
branch-and-bound.

### Notes

`Set.Icc_subset_Icc` plus `Set.Ioo_subset_Icc_self` keep the mean value point
inside the *enclosing* box rather than only inside the segment, which is what
makes the derivative enclosure apply.

The midpoint is `fdiv (lo + hi) 2`, floored, and `omega` proves it lies
between `lo` and `hi` directly — no rounding case analysis needed, because the
expansion point may be any point of the box, not specifically the centre.

`scripts/check.sh`: PASS.

## Iteration 38 — T7b3d: the split tree

### The covering argument, in nine lines

`checkTree_sound` walks the branch-and-bound's own tree. At a split the two
children agree with the parent except in one coordinate, where one takes
`[lo, m]` and the other `[m, hi]`; `mem_setHi_or_setLo` says any real in
`[lo, hi]` is in one of them, by `le_or_gt`. That is the entire covering
proof.

A flat list of 6,748 leaves would instead have left "these tile the initial
box" as an open obligation — genuine combinatorial work, and one that would
have had to be redone whenever the partition changed. The tree costs nothing
extra on the generator side, since it is what the search already builds.

Both the leaf test and the property proved are parameters, so nothing about
`π` or `g` appears in the lemma.

### Boxes on the feasibility boundary

The centered bound needs the box strictly inside the region where all five
logarithms have positive arguments — its derivatives do not exist otherwise.
Boxes touching `c = a/2` are therefore not covered by it, and those boxes
exist in the cover.

`fpiBestL` takes whichever bounds are available and intersects them when both
are:

* interior box `[0.3,0.31] x [0.4,0.41]` — plain gives `[-0.0221, +0.0057]`,
  centered `[-0.0097, -0.0066]`, and the meet is the centered one;
* boundary box `a = 0.4`, `c ∈ [0.2,0.21]`, where `c - a/2` starts at exactly
  zero — `feasOK` is `false`, the centered bound is absent, and the plain one
  still answers `[-0.2368, -0.1703]`, containing the true range
  `[-0.2231, -0.1858]`.

So the two bounds are complements rather than alternatives, and the leaf test
does not have to know in advance which kind of box it has.

`fpiCenteredG` wraps the centered bound with its guard *checked* rather than
assumed, so its soundness lemma carries no side conditions — the same shape as
`flogQ` and for the same reason. The unguarded `fpiCenteredL` is kept as it
was, so the pin written for it in iteration 37 still names exactly what it
named.

`scripts/check.sh`: PASS.

## Iteration 39 — T7b3e: the leaf test and the bridge

`DenseClaim thr a b w` is the statement the cover exists to prove: *if* the
point is strictly inside the accumulator feasibility region, *then*
`g(a) + π(a,b) + π(b,w) < thr/scale`. The paper's `-7.68·10⁻⁸` is
`thr = -76800000000000000000000` at `scale = 10^30`.

### A leaf passes for one of two reasons

The verifier accepts a box either because the bound is below the threshold or
because the box misses the feasibility region entirely — `π` is `-∞` off
`a/2 ≤ c ≤ 1-a/2`, so there the claim is vacuous. Both had to be expressible.

`infeasible` is the negation of the verifier's `intersects_feasible_region`,
four `Int` comparisons, and `infeasible_vacuous` turns it into a contradiction
with the antecedent of `DenseClaim`. It uses `A.lo` rather than `A.hi` in both
`a`-clauses, which is the direction that makes the implication hold for *every*
point of the box rather than some.

`leafOK` is then `infeasible ∨ (0 < u ∧ bound < thr)`, and `leafOK_sound`
chains: `gBA a ≤ gObj (u/scale) a` by the witness bound (which is where
`0 ≤ a ≤ 1` is needed, and why `gBA_le` carries it), `piBA = piEval` twice by
the T7b2 identity, and the enclosure's upper end below the threshold.

### Tested, both branches

On `a ∈ [0.032,0.035]`, `b ∈ [0.19,0.195]`, `w ∈ [0.099,0.102]` with witness
`u = 0.33`, the enclosure is `[-0.147657, -0.118343]` against a true range of
`[-0.144688, -0.123231]`, and `leafOK` returns `true`. Moving `w` to
`[0.002,0.003]` — below `b/2` — makes `infeasible` fire and the leaf passes
vacuously.

### The bridge

`denseTail_of_checkTree` composes `checkTree_sound` with `leafOK_sound`:
a tree that checks out proves `DenseClaim` at every point of the root box.
Everything specific to this certificate lives in `leafOK`; the covering
argument still knows nothing about `π` or `g`.

So the Lean side of T7b is complete. What remains is data: run the
branch-and-bound, emit the tree and the log table, and instantiate the bridge
at the root box `[1/125,1] × [0,1] × [1/500,13/125]`.

`scripts/check.sh`: PASS.

## Iteration 40 — T7b3f: the mirror, and a bound that was missing

### The search must be driven by the bound Lean checks

The original verifier bounds with 60-digit mpmath intervals; Lean bounds with
`Fix` at scale `10^30` and the odd series. Both are sound and neither
dominates, so a partition chosen by one is not guaranteed to satisfy the other
— and the mismatch would only surface at the end of a twenty-minute kernel
run, box by box.

So `scripts/fixmirror.py` is a line-by-line port of the `Fix` layer, and
`scripts/check_mirror.py` compares it to Lean numerator by numerator:
**22 cases, 0 mismatches**, covering `log2`, the odd series, `fxlogx` including
the `[0, x]` case, both `π` bounds, the `meet`, `fgObj`, and `leafOK` itself.

### What the first search run found

It stalled immediately: depth 60 at 4,000 nodes. The cause was not slow
convergence but an unavailable bound — `leafValue` returned `None`.

On a box straddling the feasibility boundary, `c - a/2` runs negative over
part of it, and `xlogx` of a negative enclosure does not exist. The root box
has `β ∈ [0,1]`, which straddles at almost every `α`, so nearly every early
box was unbounded and the search could only subdivide.

The original verifier clamps — `max(ZERO, b_lo - a_hi*HALF)` — and I had read
that line without registering that it was load-bearing rather than defensive.

`clampLo` raises the lower end to zero, and `fpiEvalClampL` uses it for the two
derived quantities. Soundness needs the feasibility hypothesis, which is
exactly what `DenseClaim` supplies: a point the claim says anything about has
`c - a/2 > 0`. If the box has no feasible point at all, the clamped enclosure
may be wrong and it does not matter — the claim is vacuous there.

`leafOK`'s *body* changed to use the clamped bounds; its statement and
`leafOK_sound`'s statement did not, so every pin still names what it named.

After the fix the search behaves: depth 29 rather than 60, and leaves closing
steadily.

### Where it stands

The run is in progress. At 34,000 nodes it had 16,993 leaves and 18,173
distinct logarithm arguments — already more than the mpmath partition's 13,513
leaves and 11,210 logs, which is expected since the `Fix` bound is looser in
places. The final counts decide the kernel cost, and are the first thing to
record next iteration.

`scripts/check.sh`: PASS.

### Iteration 40, continued — the tree closes, and two tuning findings

**The certificate exists under the Lean bound.** The first full run closed:

    status closed   nodes 100,525   leaves 50,263   depth 31
    distinct logs 41,024
    max_accepted_upper  -1.552e-07   against threshold  -7.68e-08

So `DenseClaim` is provable from a `Fix`-driven partition. But 50,263 leaves
against the mpmath partition's 13,513 is roughly two hours of kernel time,
which is not a per-build cost worth paying. Two things turned out to be wrong,
and measuring told me which.

**The tangent point.** Comparing my `π` bound to the verifier's on 400 random
boxes: median gap `1.4e-14`, but maximum `4.1e-2`, and every large gap was on
a *wide* box. That is the signature of the tangent lower bound for `x log x`
being taken at the midpoint, where the true minimum is at `1/e`.

The fix needed no new proof at all. `xlogx_tangent` holds at *every* positive
expansion point, so the point is a free choice, and clamping `1/e` into the
box makes the bound exact in all three cases: right of `1/e` the tangent at
`lo` bottoms out at `f lo`, left of it the tangent at `hi` bottoms out at
`f hi`, and if `1/e` is inside, the tangent there is horizontal at the true
minimum. So the exact range comes out without any monotonicity lemma, and
`einv` needs no proof because it is a heuristic, not a bound.

After the change the same 400 boxes give median gap `1.3e-14` and **maximum
`2.0e-14`** — the two bounds now agree to rounding everywhere.

**The split heuristic.** That change alone did not shrink the tree, which was
the clue that the `π` bound had not been the binding constraint. I had been
normalising each axis by its own initial width; the verifier normalises `α`
and `β` raw and `ω` by `125/99`. Since `ω` spans only `0.102`, my version
weighted its width about ten times too heavily and subdivided it that much
harder — multiplying the whole subtree.

Also worth recording because it misled me for a while: in a full binary tree
the leaf count is always `(nodes+1)/2`, so mid-run leaf counts carry no
information about whether one heuristic is beating another. Only the total at
termination does.

Re-running with the verifier's normalisation.

`scripts/check.sh`: PASS.

### Iteration 40, concluded — the partition matches the verifier's

With the exact tangent and the verifier's split normalisation:

    status closed   nodes 27,025   leaves 13,513   depth 29
    distinct logs 24,841
    max_accepted_upper  -7.686076567388105e-08   threshold  -7.68e-08

**13,513 leaves — the same number the mpmath run produced** when driven by the
`p`-level mean-value bound alone, and `max_accepted_upper` agrees with its
`-7.686951e-08` to four significant figures, the difference being the `Fix`
layer's coarser precision. Two independently written bounds, in different
arithmetic, partitioning the same region into the same number of pieces is
about as good a cross-check as this step admits.

Against the first run: leaves `50,263 → 13,513`, logs `41,024 → 24,841`.

Estimated kernel cost now: table `24,841 × 42 ms ≈ 17 min`, lookups
`13,513 × ~45 × 2.3 ms ≈ 23 min`, box arithmetic `≈ 16 min`, plus about five
minutes elaborating the two literals. Call it an hour — far better than two,
still a lot to pay on every build, so T7b3g has to decide where it lives.

One more measurement is running: `fpiBestClampL` currently computes *both*
bounds and intersects them, which roughly doubles the logarithm requests. If
taking the centered bound alone where it exists does not inflate the tree,
that halves the dominant cost.

### Iteration 40, tuning

**The `meet` is free.** `fpiBestClampL` computes both bounds and intersects
them, which looked like it should double the logarithm requests. It does not:
on a sample box, 13 requests with the meet against 12 without. The centered
bound's gradient enclosures already ask for logarithms at `C`, `1-C`,
`C - A/2`, `1 - C - A/2` and `1 - A`, which is exactly the set the plain
bound's `fxlogx` needs; only the tangent point is extra. Re-running the whole
search without the meet gave byte-identical totals. So the meet stays — it is
strictly tighter at no cost.

**`fGolayG` was costing sixty multiplications a leaf.** `pow u 8`, `pow u 12`,
`pow u 16`, `pow u 24` are 60 muls between them. In `v = u^4` the enumerator
is `1 + 759v² + 2576v³ + 759v⁴ + v⁶`, and Horner in `v` needs six. `u` is
always a point at this call site, so repeating `v` costs no width.

The partition is unchanged by it — 27,025 nodes, 13,513 leaves, 24,841 logs,
the same `-7.686076567388105e-08` — which is the right outcome: the rounding
difference between sixty roundings and six is far below anything the cover
depends on, so this is purely a cost change. It takes the estimated box
arithmetic from about 42 minutes to about 28.

`scripts/check.sh`: PASS.

## Iteration 41 — T7b3g: emission, validated small first

`scripts/emit_dense_tail.py` writes three Lean files: the logarithm table, the
split tree, and the theorem. Rather than generate three megabytes and find out
at the end of an hour whether the format works, it has a `--small` mode that
does the same thing over a sub-region.

**The table holds only what the leaves ask for.** The search evaluates
`leafOK` at internal nodes too, to decide whether to split, but Lean never
does — it recurses there. So the emitter re-walks the *final* tree and records
logarithms at leaves only.

**Choosing the small region took two tries**, both informative. The first,
`ω ∈ [0.16, 0.17]`, hit the depth cap immediately — because the dense tail
only claims negativity for `ω ≤ 0.104`, and outside that the objective is
genuinely positive. A correct refusal, not a bug. The second closed in a
single node, too trivial to exercise anything.

The region now used is `α ∈ [0.01,0.2]`, `β ∈ [0.1,0.3]`, `ω ∈ [0.05,0.104]`,
chosen to touch both accumulator boundaries — `β = α/2` at `α = 0.2` and
`ω = β/2` at `β = 0.1` — so the clamped bound and the vacuous-leaf branch are
both on the path. It emits 2,831 nodes and 3,530 logarithms, 531 KB.

Generated files need their own `set_option maxRecDepth`: a 3,530-node literal
exceeds the default during elaboration, not during reduction.

Build of the small cover is running.

### Iteration 41, concluded — the small cover is proved

    Built SpinCodes.Cover.DenseTailSmall (831s)

    denseSmall     [propext, Classical.choice, Quot.sound]
    logsSmall_ok   does not depend on any axioms
    treeSmall_ok   [propext]

So the pipeline works end to end: a table Lean re-derives, a tree Lean walks,
a covering induction, and a leaf test — giving a kernel-checked strict
negative gap over a sub-region of the certificate's box, at the paper's own
constant `-7.68·10⁻⁸`. The claim is not vacuous: `(a,b,w) = (0.1, 0.2, 0.104)`
is feasible and inside the region.

`logsSmall_ok` depending on *no axioms* is worth noting — the table is
verified by pure computation, with no classical reasoning anywhere.

### What it costs, and where that leaves the covers

831 s for 2,831 nodes and 3,530 logarithms. Scaling by 9.5x in leaves and 7x
in logarithms puts the full cover near **85 minutes**.

That settles a question that had been open since iteration 33: the covers
cannot live in `scripts/check.sh`. The loop runs it every iteration, and even
the *small* cover would add fourteen minutes. So `scripts/check-cover.sh`
builds them on demand, and the honest consequence is recorded there and here:
**a cover's claim is certified only when that script has passed**, and each
passing run must be logged with its date and module.

The default check is unaffected — `SpinCodes.lean` does not import `Cover`,
and `scripts/check.sh` still passes in its usual time.

Full emission is running.

## Iteration 42 — the full certificate is emitted and building

    search closed: 27,025 nodes
    leaf-requested logarithms: 23,541
    DenseTailLogs.lean  2.78 MB
    DenseTailTree.lean  1.06 MB

23,541 rather than the 24,841 the search touched: the difference is the
internal nodes, where the search evaluates `leafOK` to decide whether to split
but Lean only recurses. Emitting from the final tree rather than from the
search's trace drops them.

`thr` and `terms` moved to `Cover/Params.lean` so both covers share them.

### The bridge was checked before the expensive part

`denseTail_real` restates the result over the paper's own intervals and
constant, and its proof is pure arithmetic — every endpoint is exact at
`scale = 10^30` (`10^30/125 = 8·10^27`, `10^30/500 = 2·10^27`,
`13·10^30/125 = 104·10^27`), so no second enclosure is involved.

Rather than discover a broken `rw` at the end of an hour and a half, that
proof was type-checked first in isolation, with `denseTail` assumed as a
hypothesis. It goes through. So the only thing the long build can still fail
on is a leaf.

The full build is running.

### Iteration 43 — the log table needed `maxHeartbeats 0`

The tree built — `Built SpinCodes.Cover.DenseTailTree (721s)` — but the table
failed after 82 s:

    error: (deterministic) timeout at «synthesize pending MVars»,
    maximum number of heartbeats (200000) has been reached

Two separate limits, and I had only raised one. The generated files carry
`set_option maxRecDepth`, which the *small* cover needed; the full table needs
`maxHeartbeats 0` as well. The small one never hit it because 3,530 entries
elaborate inside the default budget and 23,541 do not.

Worth noting where it failed: `synthesize pending MVars`, during elaboration
of the literal, not during the `decide`. The cost that overran was building
the 2.78 MB term, not checking it.

Both generated files and the emitter now set both options. The tree's `.olean`
survives — only the table module changed — so the rebuild skips the 721 s
already spent.

## Iteration 46 — T7c reconnaissance, while the dense tail builds

The dense-tail build is at 167 CPU-minutes with memory stable at 604 MB — slow
but not pathological, and abandoning it would discard three hours and require
93 minutes just to re-elaborate the table before reaching this point again. So
it runs on, and this tick did work that touches no Lean source.

### What T7c is actually up against

`golay_ba3_concave_majorant.json` holds the majorant certificate:

    accepted_boxes 101,940    deepest 55    interval_dps 70
    largest_accepted_upper  -2.2000492744786192e-09
    left_segments 15          weight interval [13/125, 112/125]

Two numbers make this harder than the dense tail, not easier:

* **101,940 boxes against 13,513** — seven and a half times as many.
* **a margin of `2.2e-9` against `7.69e-8`** — thirty-five times thinner, so
  each logarithm needs more series terms, not fewer.

The paper describes a *refined* majorant with 39 affine supports and 10,721
box checks, which is a different and much smaller object than this stored
31-support certificate (15 left, their reflections, and the central line). The
refinement adds left supports at slopes `0.09, 0.07, 0.03, 0.01` to the 15
here, giving 19 + 19 + 1 = 39.

I checked the paper's intercept recipe — an outward rational upper bound on
`ln(1+e^{-s}) - (ln 2)/2 + 10⁻⁵` — against the stored values: at slope `0.1`
it gives `0.297832` against `0.297841148`, at slope `0.2` `0.251575` against
`0.251585223`. Agreement to about `10⁻⁵`, which is the slack term itself, so
the recipe does generate them and the four new supports can be constructed
rather than guessed.

### What this implies

At the current architecture's cost, 101,940 boxes is not reachable, and even
10,721 at a `2.2e-9` margin needs more terms per logarithm than the dense
tail's ten. So T7c depends on D3's chunking rather than merely benefiting from
it, and the first thing to establish is which of the two box counts the
refined claim actually needs.

## Iteration 47 — the monolithic cover was killed at 16 GB

Both data modules finished — `DenseTailTree` and `DenseTailLogs`, 11,766 s
each — and then `DenseTail.lean`, which holds the two `decide`s, climbed to
**16,185 MB resident**. I stopped it rather than let it thrash the machine.

That number is the useful part of the run. It says the problem is not only the
246 MB table literal but also `tree_ok`: one `decide` over 13,513 leaves is
itself an unbounded object. The small cover's 1,416 leaves checked in 831 s at
ordinary memory, so the unit that works is about a thousand leaves.

`DECISIONS.md` D3 is updated accordingly: chunk the table *and* the tree.

**The table half is done.** `treeLogN` composes a list of chunks by trying each
in turn, and `treeLogN_oracle` derives `Oracle` from `checkAll`. Soundness did
not change shape at all — a lookup still ends at a node whose key was compared
equal inside a chunk that was verified, which is the same argument
`find_check` already made for one tree. It built first try.

**The tree half is next.** `checkTree` on a `.split` is definitionally the
conjunction of its children, so splitting at a fixed depth and proving each
subtree in its own module needs a rewrite to assemble, not a new argument. The
plan is to get that assembly working on the *small* cover first, where a wrong
guess costs minutes rather than hours — the same discipline that caught the
emission format earlier.

`scripts/check.sh`: PASS.

## Iteration 49 — the chunked cover works end to end

The small cover now builds as eight independent parts and assembles:

    PartSmall0..7   158-410 s each, bounded memory
    DenseTailSmall  27 s

    denseTailSmall   [propext, Classical.choice, Quot.sound]
    partSmall0_claim [propext, Classical.choice, Quot.sound]
    logsSmall0_ok    does not depend on any axioms

Against the monolithic form, which reached 16 GB before being killed.

### Two failures, both mine

**The first rebuild reported two parts failing** with
`failed to read file ... Factors.olean.private`. That was not a regression: I
had run `scripts/check.sh` concurrently, and its `lake build` raced the parts
build over Mathlib's `.olean` files. The same race had already produced a
spurious `CHECKS FAILED` a few minutes earlier. Re-running with nothing else
touching the build, all eight parts passed. The loop prompt now says never to
run `check.sh` alongside another `lake build`.

**The assembly then failed on unification.** `claim_split_fst` was applied with
only `m` given, leaving `A`, `B`, `W` as metavariables, and Lean cannot solve
`setHi ?A m ≡ aSmall0` — `setHi` is not injective as far as unification is
concerned. Passing the boxes explicitly fixes it, and the emitter already
knows them, since it threads them through `cut`. The generated term is uglier
but the defeq check it asks for is trivial: `setHi ⟨lo,hi⟩ m` reduces to
`⟨lo,m⟩`, which is exactly the part's own box.

### Why the claim-level split was the right choice

Splitting at the level of `checkTree` would have forced every part to share one
oracle, and therefore one table loaded into every module. Splitting at the
level of the conclusion lets each part carry its own table: `PartSmall0`
imports nothing from `PartSmall1`. The cost is duplicated logarithms — 5,427
across eight parts against 3,530 in a single table, 54% overhead — which is a
good trade for bounded memory and independent caching.

Regenerating the full cover in parts now.

`scripts/check.sh`: PASS.

## Iteration 50 — the dense-tail cover is proved

    44/44 parts built, root assembled in 50 s, DenseTailReal in 43 s
    wall clock 21:52 -> 22:55, about 63 minutes, memory bounded throughout

    denseTail_real  [propext, Classical.choice, Quot.sound]
    denseTail       [propext, Classical.choice, Quot.sound]
    logs0_ok        does not depend on any axioms
    tree0_ok        [propext]

`eq:structured-ba-dense-tail` is now a kernel-checked theorem:

    g(α) + π(α,β) + π(β,x) < -7.68·10⁻⁸

at every feasible point with `α ∈ [1/125, 1]`, `β ∈ [0, 1]`,
`x ∈ [1/500, 13/125]` — the paper's own region and its own constant. The
statement is not vacuous: `(0.05, 0.1, 0.06)` is feasible and inside it.

### What it took

    13,513 leaves   27,025 tree nodes   37,872 logarithms across 44 parts
    44 modules, averaging 307 leaves each

Against the monolithic attempt, which reached 16 GB and had to be killed after
three hours. The difference is entirely structural: nothing about the
mathematics changed, only what is asked of the elaborator at one time.

### The batching detail

Lake in this version has no job-count flag — `-j` is rejected outright — and
would default to all twenty cores. Parts peaked near 4 GB each, so twenty at
once would have exceeded the machine's 31.7 GB and reproduced the original
failure by a different route. `scripts/check-cover.sh` therefore builds the
parts in batches of six, which held total resident memory under 9 GB.

### What is trusted

Nothing in the emitted data. `logsN_ok` recomputes every table entry with
`flogQ`; `treeN_ok` re-runs the leaf test at every leaf; `claim_split_*`
supplies the covering argument from `le_or_gt`. The Python mirror chose the
partition and the witnesses and is checked against Lean on 22 cases, but even
if it were wrong the worst it could do is produce a certificate Lean rejects.

`scripts/check.sh`: PASS. `scripts/check-cover.sh full`: PASS (2026-09-23).

## Iteration 51 — T7c: what the majorant actually costs

### What the claim reduces to

`â_BA` is the lower envelope of affine supports, so `a_BA ≤ â_BA` needs
`a_BA ≤ ℓ_i` for every support. But on a given `x`-box only the *smallest*
line binds — dominate that one and the rest follow. So each segment
contributes one 3-D cover of

    g(a) + π(a,c) + π(c,x) - (s·x + intercept) < 0

over `a, c ∈ [0,1]` feasible and `x` in that segment's active interval.
`scripts/majorant.py` runs it, driven by the same `fixmirror` bound Lean
checks.

### The measurement

Segment 14 (slope `0.05`) passed 30,000 nodes without closing — more than the
*entire* dense-tail cover (27,025 nodes) for one segment out of nineteen.

That is not a bug in the setup. Computing the true residual on a coarse grid:

    segment 0   slope 1.666   residual  -1.9e-3
    segment 8   slope 1.0     residual  -2.7e-5
    segment 14  slope 0.05    residual  -2.0e-5

The support lines are near-tangent to `a_BA`, so the supremum is attained in
the interior rather than being comfortably negative. That is what the stored
certificate's depth of 55 was telling me.

### Why the dense tail was easier than it looked

Its `x ∈ [0.002, 0.104]` did not merely make the `x` axis short — feasibility
`c/2 ≤ x` also forced `c ≤ 2x ≤ 0.208`, collapsing the `c` axis to a fifth of
its range. Here `x ≈ 0.49` allows `c` up to `0.98`, so nearly all of `[0,1]²`
in `(a,c)` is in play. The margin is actually *larger* here (`2e-5` against
`7.7e-8`); it is the volume and the attained supremum that cost.

### The structural route the paper points at

Exchanging the suprema,

    sup_x [a_BA(x) - s·x] = sup_{a,c} { g(a) + π(a,c) + sup_x [π(c,x) - s·x] }

and the inner supremum is one-dimensional in `x`. Its value is the Legendre
transform of the accumulator exponent — which is almost certainly why the
paper's intercept recipe is `ln(1+e^{-s}) - (ln 2)/2 + 10⁻⁵`. That form is a
conjugate, not a fitted constant.

Taking that route would drop the cover from three dimensions to two and
replace the tangency in `x` with an exact maximisation. It is the difference
between a feasible task and an infeasible one, so T7c should be scoped around
it rather than around brute force.

A longer probe of segment 14 is running to put a number on the brute-force
cost before committing.

`scripts/check.sh`: PASS.

### Iteration 52 — segment 14 closes, and what that implies

    segment 14 (slope 0.05): closed, 64,107 nodes, 32,054 leaves,
                             max accepted upper -9.59e-09

So the brute-force route is not impossible — it terminates. It is just large:
one segment is 2.4x the entire dense-tail cover, and there are nineteen left
supports plus the central one.

Extrapolating from this segment alone would overstate the total, because the
cost tracks the residual and the residuals differ by two orders of magnitude
across segments (`-1.9e-3` at slope 1.666 against `-2.0e-5` at slope 0.05).
Probes of segments 0, 4, 8 and 11 are running to get the spread before putting
a number on it.

### On the `ln(1+e^{-s})` recipe

Worth recording what that form actually is, since I had hoped it was a
shortcut. For a random linear code of rate `R` the expected enumerator gives
exponent `ln(1+z) - ln 2 + R ln 2` at `z = e^{-s}`; at `R = 1/2` that is
exactly `ln(1+e^{-s}) - (ln 2)/2`. So the supports are the random-coding
(binomial) bound, and the `+10⁻⁵` is slack on top.

But the paper is explicit that this only *proposes* the line — "validity
follows from outward checks of `eq:structured-ba-variational` over its entire
active interval". So the random-coding form is a heuristic for choosing the
supports, not a proof of them, and the 3-D cover is the intended verification
after all. The exchange-of-suprema reduction remains mathematically available,
but it would be a *different* proof from the paper's, which is not what this
formalization is for.

## Iteration 53 — T7c blocked on cost; moving to T8

    segment 14   closed at 64,107 nodes / 32,054 leaves
    segment  0   >212,000 nodes, did not close in 30 CPU-minutes

My prediction that cost would track the midpoint residual was wrong, and
usefully so: segment 0's residual is a hundred times larger than segment 14's
and it is the more expensive of the two. A support *touches* `a_BA` somewhere
inside its active interval, and the residual goes to zero there no matter what
it is at the midpoint — so every segment pays for resolving a tangency.

Nineteen supports put the cover near a million leaves: about three days of
compute and on the order of 18 GB of `.olean`, against 455 MB for the dense
tail. That is not a reasonable thing to add to this project for one lemma.

Marked **BLOCKED on resources, not correctness**, with the numbers and the
unblocking route recorded as `DECISIONS.md` D4. The paper's argument is sound
and terminates; it is this setup's per-leaf cost that does not scale to it.

Next: T8, `lem:structured-route-domination` — the first genuinely
research-grade step, needing a Chernoff bound for a Poisson-binomial sum that
Mathlib does not have.

## Iteration 54 — T8a: one pigeonhole does two of the three jobs

`Structured/RouteDomination.lean`.

Reading `lem:structured-route-domination` closely, two of its three
ingredients are the same fact:

* the conditioning cost is `1/P[Bin(b,x) = w] ≤ b+1` when `w` is a mode;
* the method-of-types bound is `P[Bin(L,k/L) = k] ≥ 1/(L+1)`, again at a mode.

Both are the observation that a mass function on `n` points puts at least
`1/n` on its largest atom — the atoms sum to one, so the largest is at least
the average. `le_mode` states it for any `FinPMF`, in three lines, and
`inv_mode_le` gives the reciprocal form the conditioning cost wants.

Worth noticing because the paper presents them as two different arguments
("the inverse conditioning probability is at most `b+1`" and "the
method-of-types lower bound"), and the second is usually proved through
Stirling. Neither needs it.

`binPMF` carries `Bin(L,u)` as a `FinPMF` on `Fin (L+1)` so `le_mode` applies
directly; totality is `add_pow`. `types_lower_bound` is then immediate,
modulo the mode hypothesis, which is T8b.

Building the project's own `FinPMF` rather than Mathlib's measure-theoretic
probability keeps this consistent with `Prob.lean`, where the rest of the
development lives.

`scripts/check.sh`: PASS.

## Iteration 55 — T8b: the mode of `Bin(L, k/L)`

`P[Bin(L,k/L) = k] ≥ 1/(L+1)` now holds with no hypothesis left open.

### Clearing denominators first

The textbook proof runs a ratio test — `P(j+1)/P(j) ≥ 1` exactly when
`j+1 ≤ k(L+1)/L` — and then chains the resulting monotonicity in both
directions. That is a lot of bookkeeping in Lean.

Multiplying through by `L^L` instead leaves a statement with no division and
no reals at all:

    C(L,j) · k^j · (L-k)^(L-j)  ≤  C(L,k) · k^k · (L-k)^(L-k)

and that needs only two inequalities, both about factorials against powers:

    A!  ≤  j! · A^(A-j)                    (the factors j+1..A are each ≤ A)
    B! · B^d  ≤  (B+d)!                     (Mathlib's factorial_mul_pow_le_factorial)

`mode_half` is stated symmetrically in the two "sides" `A` and `B`, so the
case `j ≥ k` is the same lemma with `k` and `L-k` exchanged rather than a
second argument. Then `Nat.choose_mul_factorial_mul_factorial` converts the
factorial form to the binomial-coefficient form by multiplying up and
cancelling, and `binPMF_mode` divides back down by `L^L`.

### Two Lean notes

Mathlib's factorial notation is `scoped` in `Nat` **and** rejects a space
before the `!`. With neither `open Nat` nor the exact spacing the parser gives
"unexpected token; expected no space before", which does not obviously point
at the factorial. Writing `Nat.factorial n` sidesteps both.

`Nat.exists_eq_add_of_le` inside an induction leaves the original `≤`
hypothesis in the goal, so the induction hypothesis arrives as an implication
and needs `ih (by omega)` rather than `ih`.

`scripts/check.sh`: PASS.

## Iteration 56 — T8c1: the piece Mathlib does not have

    ∏ (1 + pᵢ·c)  ≤  (1 + q·c)^L,     q = (∑ pᵢ)/L,  c ≥ -1

At `c = e^t - 1` this says replacing the success probabilities by their
average can only increase the moment generating function, which is exactly
what lets the Chernoff bound for a Poisson binomial be read off the binomial
one. It is the only ingredient of `lem:structured-route-domination` with no
counterpart in Mathlib.

It is AM-GM with uniform weights. `Real.geom_mean_le_arith_mean_weighted`
gives `∏ zᵢ^(1/L) ≤ (1/L)∑ zᵢ` at `zᵢ = 1 + pᵢc`; the right side is `1 + qc`
by linearity, and raising both sides to the `L` finishes it. The only fiddly
part is moving between `rpow` and `pow` — `Real.finset_prod_rpow` to pull the
product inside, then `Real.rpow_natCast` and `Real.rpow_mul` to cancel
`(1/L)·L`.

The nonnegativity side condition needs `c ≥ -1` and `pᵢ ≤ 1`: for `c < 0`,
`pᵢc ≥ c ≥ -1`. Both hypotheses are used, neither is padding.

Checked on 20,000 random cases with no violation, and equality holds exactly
when the `pᵢ` agree — so the bound is tight rather than vacuously true.

`scripts/check.sh`: PASS.

## Iteration 57 — T8c2: the Chernoff bound, without ever writing `exp`

The task as scoped said "`P[S = k] ≤ exp(-L·D(k/L‖q))`, from T8c1 at
`c = e^t - 1`". Following that literally would have meant carrying an
exponential and a KL divergence through the optimisation and then unwinding
both again to meet T8b. Working out where the two halves have to meet first
showed they never need to be written that way.

### The observation

Optimising the free parameter in the Chernoff bound gives

    P[S = k]  ≤  q^k (1-q)^(L-k) · L^L / (k^k (L-k)^(L-k)),

and that *is* `exp(-L·D(k/L‖q))` — but as a ratio of powers, with no
transcendental function in it. T8b, meanwhile, says

    C(L,k) · k^k (L-k)^(L-k) / L^L  ≥  1/(L+1),

whose left side is the exact reciprocal of the Chernoff factor. So the two
cancel directly:

    P[S = k]  ≤  (L+1) · C(L,k) q^k (1-q)^(L-k),

which is the paper's "pointwise at most `L+1` times the iid Bernoulli-q law".
No `Real.exp`, no `Real.log`, no divergence — and no case split on whether
`k/L` is above or below `q`, because `z^k · P[S=k]` is a sub-sum of the
generating function for *every* `z > 0`, not only on one side of the mean.

### The four steps

1. `poissonBinom` — the law, indexed by the success *set* rather than by
   `Fin L → Bool`. That choice means one Mathlib lemma, `Finset.prod_add`,
   discharges both the total-mass obligation and the generating identity.
2. `poissonBinom_mgf` — `∑_T P(T) z^|T| = ∏ (1 - pᵢ + pᵢz)`: the same
   `Finset.prod_add` with `pᵢz` in place of `pᵢ`.
3. `poissonBinom_chernoff` — `z^k · P[S=k]` is a sub-sum of the left side
   (`Finset.sum_le_sum_of_subset_of_nonneg`), the right side is
   `∏ (1 + pᵢ(z-1))`, and T8c1 dominates it by `(1 + q(z-1))^L`.
4. `opt_algebra` / `poissonBinom_prob_le` — substitute
   `z = k(1-q)/(q(L-k))`, at which `1 + q(z-1) = (1-q)L/(L-k)`. The rest is
   `field_simp` once `(1-q)^L` and `(L-k)^L` are split as `^(L-k) · ^k`;
   `ring` cannot merge variable exponents, so the split has to be done by
   hand, but after it both sides match atom for atom.

### Validation

4,000 random Poisson binomials (`L` from 2 to 14, exact DP for the law): the
bound never fails, worst observed ratio `0.617`. The unscaled ratio
`P[S=k] / Bin(L,q)[k]` reaches `2.19`, so the `L+1` factor is load-bearing —
the same statement with `1` in its place would be false. Both hypotheses
`0 < k` and `k < L` are needed (the optimal `z` degenerates at either end).

Twice now `field_simp` has closed a goal outright and left the following
`ring` erroring with "no goals"; worth remembering that the error points at
the *next* tactic, not the one that over-applied.

`scripts/check.sh`: PASS.

## Iteration 58 — T8d: assembling route domination

`E_route[F] ≤ (b+1)^Q (L+1)^b E_Ber[F]`, for every nonnegative `F`, plus the
`exp(o(N))` corollary. New file `Structured/Domination.lean`.

### Proving the stronger statement is the easier one

The lemma quantifies over *all* nonnegative `F`, which is equivalent to a
*pointwise* bound `P ≤ c·Q` on the mass functions. Pointwise is the better
thing to carry: it composes (`Dominates.trans`, costs multiply) and it
multiplies over independent coordinates (`dominates_piPMF`), neither of which
is natural at the level of expectations. Expectations are recovered once, at
the end, by `Dominates.expect_le`.

So the file is a small calculus of domination, and the paper's two mechanisms
are two lemmas in it:

* `dominates_condition_le` — conditioning costs the reciprocal of the
  conditioning probability, which is `≤ b+1` when the event is "the sum hits
  a mode" (T8b);
* `dominates_of_fiber_uniform` — if both laws are constant inside each weight
  class, a bound on the *class* masses is a bound on the individual masses,
  because the class cardinality appears on both sides and cancels.

### Where the modelling boundary actually falls

The dominating side is provably fibre-uniform: `poissonBinom_const_apply`
shows the iid Ber(q) mass depends only on `|T|`, and `poissonBinom_const_prob`
shows its class mass is exactly `C(L,k) q^k (1-q)^(L-k)` — so the interface is
not vacuous, the dominating law really is the binomial.

The routed side is *not* fibre-uniform, and that is not a gap in the
formalization: it is the reason the paper needs the region shuffles at all
("the deterministic transpose alone spreads each row, but it does not erase
the row's coordinate pattern"). Fibre-uniformity of the routed law is
therefore a named hypothesis of `dominates_region`, which is exactly the right
place for the skeleton to stop.

### The endpoints, and a refactor that paid for itself

`poissonBinom_le_binom` (T8c2) needs `0 < k < L` because the optimal `z`
degenerates at the ends. Both ends turn out to be AM-GM with nothing else:
`P[S=0] = ∏(1-pᵢ) ≤ (1-q)^L` and `P[S=L] = ∏ pᵢ ≤ q^L`. The first is
`prod_one_add_mul_le` at `c = -1`; the second is not of that shape at all.

Rather than duplicate the AM-GM argument, `prod_le_pow_avg` —
`∏ zᵢ ≤ ((∑ zᵢ)/L)^L` — was factored out of T8c1, and `prod_one_add_mul_le`
re-derived from it as the specialisation at `zᵢ = 1 + pᵢc`. The pinned
statement of `prod_one_add_mul_le` is unchanged, so the pin confirmed the
refactor rather than hiding it. Net effect: one AM-GM application in the
project instead of two, and the endpoints come for free.

### `exp(o(N))`

`ln((b+1)^Q (L+1)^b) = Q ln(b+1) + b ln(L+1) ≤ Lb·(ln(b+1)/b + ln(L+1)/L)`,
and both terms vanish. Only `Q ≤ L` is used — the paper's extra hypothesis
that `Q/L` is bounded below is what makes `Lb` the right normaliser
elsewhere, not what makes this estimate true, so it is not a hypothesis here.
`Real.isLittleO_log_id_atTop` plus a squeeze against `2·ln(n+1)/(n+1)` gives
`ln(n+1)/n → 0`.

### Validation

The endpoint cases, which T8c2's 4,000-case run did not cover, were checked
on their own 4,000 cases: worst ratio `0.333`, never violated.

`scripts/check.sh`: PASS.

## Iteration 59 — T9a: Krawtchouk from scratch

T9 (`app:imt-finite-transfers`) is much the largest remaining item, so it was
split into five before starting: T9a Krawtchouk, T9b Fourier/Parseval, T9c the
integer fibre bounds, T9d the transfer entries, T9e the induction. This
iteration did T9a.

### What the appendix actually needs

Everything in the fibre-bound section flows from two identities,

    k_j = 2^{-19} (n_j + ∑_w b_w K_j(w)),
    V_j = 2^{-19} (n_j² + ∑_w b_w K_j(w)²),

and both reduce to one combinatorial fact: the character sum over a weight
layer, `∑_{|x|=j} (-1)^{|x∩v|}`, depends on `v` only through `|v|`. Without
that, the identities are not even well posed — the right-hand sides are
indexed by weights, the left-hand sides by vectors.

Mathlib has neither Krawtchouk polynomials nor MacWilliams, so this is built
from nothing. Representing vectors as subsets makes weight `Finset.card` and
the inner product `|v ∩ x|`, which keeps the whole development inside
`Finset`.

### The proof

Partition the weight-`j` layer by `h = |x ∩ v|`. The summand is constant
`(-1)^h` on each part, so only the part sizes matter, and the map
`x ↦ (x ∩ v, x \ v)` is a bijection from the part onto
`powersetCard h v × powersetCard (j-h) vᶜ` — giving `C(w,h)·C(n-w,j-h)`, which
is the Krawtchouk summand.

### Lean friction, not mathematical friction

Nearly all the iteration went on two mechanical things. `Finset.card_nbij'`
states its membership side conditions against the *set* coercion `↑s`, so the
hypotheses arrive as `a ∈ ↑s` and `rw [Finset.mem_product]` does not fire until
`Finset.mem_coe` is unfolded. And its four obligations are stated about the
un-beta-reduced applications `(fun x => …) x`, so goals mentioning
`#(x ∩ v, x \ v).2` do not match hypotheses about `#(x \ v)`; each needs a
`show` to beta-reduce before `omega` or `rw` will engage. Neither shows up in
the error message as what it is.

### Validation

Exhaustive: every `n ≤ 8`, every one of the `2^n` subsets `v`, and every `j`
from `0` to `n+1` (so the out-of-range case is covered too) — 4,608 instances,
zero mismatches. As an independent check of the *definition* rather than the
proof, the generating-function characterisation
`K_j(w) = [t^j](1-t)^w (1+t)^{n-w}` was verified to agree on all `n ≤ 6`.

`scripts/check.sh`: PASS.

## Iteration 60 — T9b: Parseval without the syndrome map

`|D|·k_j = ∑_{v∈D} K_j(|v|)` and `|D|·V_j = ∑_{v∈D} K_j(|v|)²`.

### The choice that shrank the task

`V_j` is defined in the paper as the sum of squares of the *syndrome fibre
sizes*, which sounds like it needs the 19×128 matrix, the syndrome space, and
a quotient. It does not. Two words share a syndrome exactly when their
symmetric difference is in the kernel, so

    ∑_s f_j(s)² = #{(x,y) : |x| = |y| = j, x ∆ y ∈ ker},

and the kernel is determined by the dual code. Both identities therefore need
nothing but a subgroup `D` of `(Finset (Fin n), ∆)`. No matrix, no syndrome
space, no quotient, and no linear algebra over `ZMod 2` — the development
stays inside `Finset`, where T9a already lives.

The same choice makes the two proofs nearly identical. Inversion swaps the
order of `∑_x ∑_{v∈D}`; Parseval swaps `∑_{(x,y)} ∑_{v∈D}` and then factors
`∑_x∑_y χ(v,x)χ(v,y)` as a square. The only extra ingredient is that `χ` is
multiplicative in its *second* argument as well as its first — one lemma,
since both come from `|·|` being additive mod 2 under `∆`.

### Orthogonality

The standard argument: if some `u ∈ D` has `χ(u,x) = -1`, reindex the sum by
`v ↦ u ∆ v` (a bijection of `D`, by closure and involutivity) to get
`S = χ(u,x)·S = -S`. `Finset.sum_nbij'` takes plain Finset membership, unlike
`Finset.card_nbij'` from the previous iteration, which takes the set
coercion — worth knowing, since the two look interchangeable.

`Finset`'s `∩` is not syntactically `⊓`, so Mathlib's
`inf_symmDiff_distrib_left/right` do not fire on it; two one-line `ext; simp;
tauto` lemmas were cheaper than converting.

### Validation

1,766 cases: random subgroups of `𝔽₂^n` for `n ≤ 9` (dimension up to 4), every
`j`, both identities checked against brute-force enumeration. Zero mismatches.
997 of the cases had a nontrivial code *and* a nonempty orthogonal layer, so
the agreement is not `0 = 0`.

`scripts/check.sh`: PASS.

## Iteration 61 — T9c1: three of the four fibre bounds

`a_j`; the Fourier bound `|D|·f_j ≤ ∑_{v∈D}|K_j(|v|)|`; and constant-weight
packing. The fourth, the variance bound, is split off as T9c2 because it is
the only one that needs *all* the fibre sizes at once rather than one fibre.

### Generalising inversion paid for itself immediately

The Fourier bound is the triangle inequality applied to inversion at a
*general* syndrome, so `card_orth_layer` from T9b (inversion at zero) was not
enough. Generalising it gives

    |D| · |fibre through x| = ∑_{v∈D} χ(v,x) · K_j(|v|),

and `card_orth_layer` then falls out as the `x = ∅` case, so the file got
shorter rather than longer. The pin confirmed the old statement was unchanged
by the refactor.

Defining the fibre as `{y : |y| = j, x ∆ y ∈ ker}` keeps the C-free setting of
T9b: "the fibre through x" needs no syndrome to name it.

### Packing, stated once

The appendix gives the packing bound with `h = min(j, 128-j)` and a remark
about complementing when needed. Rather than prove it twice, it is stated for
an arbitrary constant-weight set with a pairwise-distance guarantee. The fibre
gives one instance; the image of the fibre under complementation — same
pairwise distances, weight `n-j` — gives the other. The complemented corollary
is then eight lines.

One hypothesis was dropped after checking it was not needed: `r ≤ h`. With
`r > h` the bound is still true (and `omega` still closes the disjointness
step), because `2r < d ≤ |y ∆ y'| ≤ 2h < 2r` is already contradictory.

### Validation

~20,000 checks over random binary codes (`n ≤ 9`, dimension up to 4), every
weight, every fibre, and for packing every valid `r`. Zero violations. Each
bound is *tight* in hundreds of cases — 471, 1,216 and 836 respectively — so
none of the three is vacuously loose.

### Note

The `bash` heredoc failed to parse again on a file with mixed quoting; the
Write tool was the reliable path. Third occurrence, so: for Lean files with
doc comments, prefer Write over heredocs.

`scripts/check.sh`: PASS.

## Iteration 62 — T9c2: the variance bound, without the square root

### Stating it before the rounding

The appendix writes this bound as

    ⌊(a_j + ⌈√((M-1)(M(V_j - k_j²) - a_j²))⌉) / M⌋,

but the square root and the two roundings are not part of the argument — they
come from solving a quadratic and then rounding outward. What Cauchy-Schwarz
actually proves is

    (m+1)·t² + a_j²  ≤  2·a_j·t + m·(V_j - k_j²),

so that is what is stated. No `Nat.sqrt`, no floors, no ceilings, and nothing
to get wrong at the rounding step.

Written additively (`a = A + t`, `S = B + t²`, and `V_j = k_j² + ∑|Φ|²`) the
whole development stays in `ℕ` with no truncated subtraction anywhere, and the
core lemma reduces to exactly `A² ≤ m·B`, which is Mathlib's
`sq_sum_le_card_mul_sum_sq` on `F.erase i₀`.

Relaxing `m` upward is sound here, and worth noting: Cauchy-Schwarz gives
`(a-t)² ≤ |others|·(S-t²)`, and enlarging the count only weakens the
constraint. So `m+1` need only be an *upper bound* on the number of fibres,
which is what makes `M = 2^19 - 1` usable without proving the fibre count
exactly.

### The partition

This is the only one of the four bounds that needs all the fibre sizes at
once, so most of the file is the partition of the non-kernel layer into
fibres — still with no syndrome map, fibres being cosets named by any of their
members.

### Two timeouts worth remembering

`hrw ▸ h` on `Finset` goals blew the `whnf` heartbeat limit; rewriting the
hypothesis with `rw [hrw] at h` instead is predictable and fast.

More important: proving symmetric-difference identities by
`ext a; simp [Finset.mem_symmDiff]; tauto` also timed out once the identity
had four distinct variables. `(Finset α, ∆)` is a group, so
`symmDiff_symmDiff_symmDiff_comm` + `symmDiff_self` + `symmDiff_bot` prove
these in one `rw` each. Three such identities are now named lemmas.

### Validation

25,227 checks of the quadratic over random binary codes (`n ≤ 9`), every
weight, every fibre, and several values of `M ≥ #fibres`: zero violations,
3,843 tight. Both partition identities checked exactly on the same data. And
the appendix's own floor/square-root formula was evaluated on the same data
and also held everywhere — confirming the quadratic really is that bound
before rounding, rather than a different statement.

`scripts/check.sh`: PASS.

## Iteration 63 — T9d1: the transvection law

T9d is "the `T_j(z)` entries and the claim that they dominate the weighted
state measure". Reading the model first showed that the entries are not the
place to start: every one of them is built from a single probabilistic fact,

    L(M q) = ½ δ_q + ½ Unif(𝔽₂^s \ {0})    for q ≠ 0,

where `M = I + uvᵀ` with `u` uniform nonzero and `v` uniform in `u^⊥`. The
`½` on the lazy branch and the `1/(2M)` per target on the refresh branch are
exactly its two halves. So T9d was split: this iteration proves the law, and
T9d2 assembles the matrices on top of it.

### Stating it without `2^(s-1)`

Phrasing the law as probabilities would drag in `2^(s-1)`, which is awkward at
`s = 1` and needs division. Instead it is three exact counts over the sample
space of pairs `(u,v)`:

*  zero is unreachable from a nonzero `q`;
*  each other nonzero target is hit by exactly half of one `u^⊥`;
*  `q` itself is hit by *all* of `perp q` plus half of every other `perp`.

With `2·|u^⊥| = 2^s` these give `½ + 1/(2M)` and `1/(2M)`. No division, no
truncated subtraction, no edge case in `s`.

### One lemma does the work twice

Both `2·|u^⊥| = 2^s` and the fair-bit step are instances of: *a parity
functional splits a `∆`-closed family in half*, via the involution
`v ↦ v ∆ w` for any `w` the functional is odd on. The only real content is
finding that `w`: given `q ∉ {∅, u}` one needs `⟨u,w⟩ = 0` and `⟨q,w⟩ = 1`,
which is a two-case construction — take a single element of `q \ u` if there
is one, otherwise `q ⊊ u` and a pair `{a ∈ q, b ∈ u \ q}` works.

### A parity idiom worth keeping

Proving `Even |t ∩ (v ∆ w)| ↔ (Even |t ∩ v| ↔ Even |t ∩ w|)` by manipulating
`Even`'s existential witnesses failed — the witness is a metavariable when
`omega` runs, so it cannot close the goal. Rewriting all three with
`Nat.even_iff` first turns it into `%2` arithmetic and `omega` handles the
whole iff structure in one step.

### Validation

Exhaustive for `s = 1..7`: all four counting statements over every nonzero
`q`, zero violations. Separately, the *derived* law was checked as exact
rationals — `P[Mq = q] = ½ + 1/(2M)` and `P[Mq = q'] = 1/(2M)` for every other
nonzero `q'` — confirming the three counts really do encode the paper's
statement rather than something adjacent to it.

`scripts/check.sh`: PASS.

## Iteration 64 — T9d2: the law as an expectation

The transfer-matrix rows never use the transvection law as a distribution;
they use it as an expectation. So the useful form is

    E[g(M q)] = ½ g(q) + (1/(2M)) ∑_{w ≠ 0} g(w).

Two things about it are worth recording.

It is an **exact identity, not a bound**, and it holds for *every* `g` — no
nonnegativity, no monotonicity. That is a stronger and simpler statement than
the domination the matrix entries need, and it means the `½` and the `1/(2M)`
in the entries are not estimates at all: they are literally the two halves of
`½ δ_q + ½ Unif`, and any looseness in the matrix comes from elsewhere.

It also comes out with **no division**, as

    2 · ∑_p g(M_p q) = |u^⊥| · (M · g(q) + ∑_{w ≠ 0} g(w)),

because the three counts of T9d1 all carry a factor `|u^⊥|`. The normalised
version is a corollary; the integer version is what a later kernel
computation would want.

### Assembly note

The derivation is a fibrewise regrouping of the sample space by target, then
the two counting laws. Writing it as one long `rw` chain failed opaquely —
each step's output shape has to match the next pattern exactly. Naming the
three intermediate equalities (`hA`, `hB`, `hnz`) and finishing with a single
`ring` worked first time and is far easier to read.

`Finset.sum_fiberwise_of_maps_to` also needs its function argument explicitly
typed: passing `fun p => g (act p.1 p.2 q)` leaves `p`'s type unknown and the
projections `.1`/`.2` fail to elaborate, with an error that points at the
projections rather than at the missing annotation.

### Validation

1,235 instances: `s = 1..7`, every nonzero `q`, five random `g` per size with
values random *rationals including negatives*, checked exactly. Both the
integer and the normalised forms hold with no violations. Using negative `g`
is the point — it confirms the identity really is sign-free rather than a
bound that happens to be tight on nonnegative test functions.

`scripts/check.sh`: PASS.

## Iteration 65 — T9d3: the weighted-state abstraction

The seven coordinates `(Z, D, S_1..S_5)`, the domination relation, the
terminal bound, monotonicity, and the transfer row leaving a nonzero state.

### `D` is what makes the definition an existential

`Z` and the `S_i` are pointwise constraints — mass at zero, and mass below a
multiple of a shell's uniform measure. `D` is not: it "bounds arbitrary mass
on nonzero states", a *total* budget that may sit anywhere. So domination has
to quantify existentially over how that budget is spread, and everything else
stays pointwise. Getting this right matters, because a pointwise reading of
`D` would be a different (and false) statement.

### The shells stay abstract

Nothing in the architecture needs five shells, or the map `A`, or the specific
weights `(48,56,64,72,80)`. A `ShellSystem` is any family partitioning the
nonzero states, and the paper's weight classes are one instance. That keeps
the file independent of the 19×128 tables.

### The row is exact, and lossless

The row leaving a nonzero state is not an estimate:

    P(q, w) = ½·[w = q] + ∑_h (c_h/(2M))·Unif_h(w)

holds with equality, straight from T9d1/T9d2. The `½` is the lazy branch
sitting at `q` — a nonzero state, so it is charged to `D` — and the refresh
branch is flat at `1/(2M)`, which in shell coordinates is exactly `c_h/(2M)`.

Worth noting: the coordinate total on this row is `0 + ½ + ∑_h c_h/(2M)`,
and since the shells partition the `M` nonzero states this is `½ + ½ = 1`.
The row is a probability, so a total of exactly `1` means the abstraction
loses *nothing* on it — useful to know, since a lossy representation here
would silently weaken every later bound.

### Naming collision

`Dominates` was already taken by the route-domination calculus of T8d, which
is a different relation (pointwise on two mass functions, one scalar). The IMT
abstraction now lives in `Spin.Imt`. Two relations both fairly called
"dominates" is a sign these are genuinely separate mechanisms, not a
duplication to be merged.

### Validation

`s = 1..7`, every nonzero `q`, every target, with a random shell partition per
size: the pointwise law, the shell reconstruction, and the coordinate total
all exact. The total came out `1` on every instance.

`scripts/check.sh`: PASS.

## Iteration 66 — T9d4: the step kernel

Composing the transvection law with the syndrome shift `C x` gives the row of
`T_j(z)` leaving a state. The composition is **exact**:

    2 ∑_x ∑_p wgt(x) f(M_p q + C x)
      = |u^⊥| · ( M ∑_x wgt(x) f(q + C x)
                  + ∑_x wgt(x) ( ∑_w f(w) − f(C x) ) ).

The second bracket is the refresh branch. It is flat over *every* state with
the single state `C x` removed — not flat over the nonzero states — because
the refresh lands on `w + C x` for `w` uniform over the nonzero states, and
shifting moves the hole from `0` to `C x`.

### Three of the paper's entries fall out

*  **The zero row is exact.** `M_i` fixes zero, so that row carries no `½` at
   all: it is the input law pushed through `C`.
*  **The zero column** is `½ L_j(q;z) + (1/2M) E[z^{wt} 1(C x ≠ 0)]`, the
   indicator-of-zero case of the identity. The appendix remark that
   "termination also requires a nonzero syndrome" is exactly the `− f(C x)`
   term, so it is derived rather than imposed.
*  **The `min(m_{d,j}, ν_j)`** is the refresh term bounded two ways — by the
   total weight, and by the count of nonzero syndromes — needing only
   `0 ≤ wgt ≤ 1`, which is `z ≤ 1`.

### What the step does not need

`A` and `C` never have to be linear, or be anything but functions: `A` enters
only through the weight `z^{wt(x + Aq)}`, and `C` only as the shift. Inputs
are an arbitrary `wgt` over an arbitrary finite set, so there is no
normalisation and no probability measure — the identity is a statement about
sums. That keeps the file independent of the 19×128 tables and of the choice
of input law, and means the same lemma serves the fixed-weight `T_j`, the
occupation mixture `T_occ`, and the tilted `T_F`.

### Validation

796 nonzero-state instances over random `s ≤ 5`, `t ≤ 4`, random input sets,
random `C`, and random *signed rational* weights and test functions, plus 60
zero-row trials. All three identities exact, no violations. Signed values
again matter: they confirm these are identities, not bounds that happen to
hold on nonnegative data.

`scripts/check.sh`: PASS.

## Iteration 67 — T9e: the transfer induction

`E[z^{wt(Y)}] ≤ e_Z T^R 1`, plus mixtures for `T_occ` and the scalar envelope
`T_sc`. T9 is now complete.

### The induction is the short part

Once the one-step statement is in place the induction is ten lines: if the
step carries a measure dominated by `c` to one dominated by `c·T`, then
`step^[R] μ₀` is dominated by `(T.apply)^[R] c₀`, and `total_le` finishes.
Everything expensive was in getting the one-step statement right, which is
what the previous four iterations were.

The row-domination claim itself stays a hypothesis of `imt_moment`. That is
deliberate and is the scope boundary: it is where the verifier's numerical
envelopes (`m_{d,j}`, `ν_j`, `ℓ`) enter, and they are envelopes, not
architecture.

### Pinning the convention, not just the theorem

"Rows of a transfer matrix index the entering coordinate, so weighted measures
propagate as row vectors." A transposed `apply` would typecheck and prove an
equally plausible-looking theorem, so `Transfer.apply` and `Coords.eZ` are
pinned by `rfl`, and the numerical check verifies that iterating `apply` from
`e_Z` really equals `e_Z T^R` computed as a matrix power. That is the one
error this file could plausibly have contained.

### A `rw` trap worth recording

`rw [total_apply, Coords.total, ...]` unfolded the *wrong* `.total`: the goal
contained `T.rowZ.total` before `c.total`, and `rw` takes the first match.
The failure surfaces much later as an unhelpful `linarith failed`. Isolating
the intended rewrite in its own `have` fixed it, and is the safer habit
whenever a projection appears several times in a goal.

### Validation

200 random transfer matrices over `k ≤ 4` with rational entries, checking:
`apply` agrees with row-vector-times-matrix; `total_apply`; iterating from
`e_Z` agrees with the explicit matrix power `e_Z T^R`; and the scalar envelope
`e_Z T^R 1 ≤ λ^R` for `λ` the maximum row total. All clean.

`scripts/check.sh`: PASS.

## Iteration 68 — T10a: the KL chain-rule bound

`D(αx ‖ py) ≤ D(α‖p) + α D(x‖y)`, the likelihood cost of moving the iid
reference probability in `eq:imt-dense-exponent`. T10 split into T10a (this),
T10b (convexity and vertex bounds) and T10c (the Collatz bound).

### Why it is an inequality and not an identity

The paper calls it "the chain rule", which is exact — but only for the *joint*
experiment. A Bernoulli(α) active flag followed, when active, by a
Bernoulli(x) bit has joint divergence exactly `D(α‖p) + α D(x‖y)` against the
`(p,y)` reference, because the inactive branch contributes nothing. What the
argument actually needs is the divergence of the *observed bit*, which is a
function of that pair, so the step is chain rule followed by data processing.
Concretely that is one application of the two-term log-sum inequality, merging
"inactive" and "active-but-zero" into the single outcome `0`.

The whole analytic input is `log t ≤ t - 1`: writing `∑ aᵢ log(aᵢ/bᵢ) − A log(A/B)`
termwise and applying it to `aᵢB/(bᵢA)` gives `≥ A − B·(A/B) = 0`.

### The closed range

The paper covers `α ∈ [10⁻⁴, 1]`, and `α = 1` (full occupation) is attained,
so stopping at the open interval would have been a real gap rather than a
technicality. Stating log-sum with `0 ≤ aᵢ` rather than `0 < aᵢ` is what makes
the endpoints work: at `α = 1` the term `(1-α) log((1-α)/(1-p))` is `0 log 0`,
which is `0` both mathematically and in Lean. Only `x = 1` needed an explicit
case split, for the log-splitting step.

### Validation

200,000 random log-sum instances including zero numerators, and 200,000
chain-rule instances with `α` and `x` each forced to exactly `1` half the
time: no violations. About a quarter of the chain-rule cases are *exactly*
tight — those with `α = x = 1`, where both sides are `log(1/p) + log(1/y)` —
so the bound is attained and not merely true.

`scripts/check.sh`: PASS.

## Iteration 69 — T10b: convexity and the vertex check

What justifies the verifier "fixing a witness on each box and checking four
vertices".

### Why the change of variables is to `(α, αx)` and not `(α, x)`

In `(u, w) = (α, αx)` the three non-constant terms of the exponent are
`u·â_BA(w/u)`, `D(u‖p)` and `u·D(w/u‖y)`. The first is *linear* precisely
when `â_BA` is affine — which is what "on each affine segment" is buying, and
is the reason the segments appear at all. The third is the **perspective
transform** `u·g(w/u)` of a convex `g`, which is convex; Mathlib has no
perspective lemma, so it is proved here from Jensen with weights
`au₁/(au₁+bu₂)`.

### A test failure that was really a missing hypothesis

The first numerical run of the vertex principle reported five violations. They
were not violations: they were boxes not contained in the cone
`{0 < w/u < 1}`, where the exponent is not even real-valued — `D(x‖y)` at
`x > 1` takes a logarithm of a negative number. Restricting to boxes inside
the cone, 3,655 boxes with 40 interior samples each gave no violations.

This is worth recording because the constraint is real and easy to miss: in
`(α, αx)` coordinates the feasible region is a **cone, not a rectangle**, so
a box `[α₁,α₂]×[w₁,w₂]` is admissible only when `w₂ < α₁`. The `ConvexOn`
hypothesis of `le_max_vertices` is stated on the box itself, so it encodes
exactly this requirement rather than silently assuming it.

### Lean notes

`λ` is a reserved token in Lean 4, so `hλ₁` as a hypothesis name is a parse
error — and the error points at the *next* line, not the name.

`Set.mem_prod.mpr ⟨hx, hc⟩` leaves the pair as a metavariable and fails to
elaborate; `Set.mk_mem_prod hx hc` states the pair explicitly and works.

### Validation

60,000 convex-combination tests: the full exponent convex in `(u,w)`, the
perspective term alone, and `D(·‖p)` in its first argument — all clean; plus
the vertex principle as described above.

`scripts/check.sh`: PASS.

## Iteration 70 — T10c: the Collatz bound

`T w ≤ λ w` for a positive column `w` bounds `e_Z T^R 1`. T10 is complete.

### One identity does all the work

The whole argument is adjointness, `(c T) · w = c · (T w)`. With it, each
application of the row action becomes one factor of `λ`:

    (e_Z T^R) · w  ≤  λ·(e_Z T^{R-1}) · w  ≤ … ≤  λ^R · (e_Z · w)  =  λ^R w_Z,

and `c.total · w_min ≤ c · w` converts the pairing back to the all-ones
column. Proving the adjointness is the only real work — it is the same
double-sum exchange as `total_apply`, which is unsurprising, since `total` is
the pairing against the all-ones column.

### Two small sharpenings over the paper's phrasing

The paper writes the factor as `λ^R · (max w / min w)`. Starting from `e_Z`
only `w_Z` is ever seen, so the bound here is `λ^R · w_Z / w_min`, which is at
least as good. It is also stated multiplicatively,
`(e_Z T^R 1) · w_min ≤ λ^R · w_Z`, so no division appears at all.

`0 ≤ λ` is *derived* rather than assumed: `T w` is nonnegative when `T` and
`w` are, and `w_Z ≥ w_min > 0`, so the eigen-inequality forces it. Worth
doing, because a hypothesis that looks harmless can hide a case where the
bound is vacuous.

### Is it redundant with the scalar envelope?

No — checked. On 300 random matrices, the Collatz bound is *strictly* better
than `λ_scalar^R` (the maximum row total) in 20 of them, even though `w` here
is random rather than an approximate eigenvector. With a tuned `w` it would
win far more often; the point of the check was only to confirm the two bounds
are genuinely different, so that `T10c` is not restating `T9e`.

### Validation

300 random transfer matrices with rational entries: adjointness exact, the
Collatz bound never violated.

`scripts/check.sh`: PASS.

## Iteration 71 — T11a: empty-epoch estimates

T11 (`app:imt-fixed`) is the largest remaining item and was flagged as the
highest risk, so it was read end to end and split first: T11a empty epochs
(this), T11b the lift/projection and coarse-grained kernel, T11c the simplex
integral `K_a(θ)`, T11d the weighted-norm bound.

### The matrix exponential is not needed

The appendix states the empty-epoch mixing through `P = (I+Π)/2` and
`P^g = 2^{-g}I + (1-2^{-g})Π`. But `Π` sends every live state to the same
stationary mean, so the only thing that ever varies is a scalar, and the whole
analysis collapses to

    m_{g+1} = (m_g + μ)/2.

That is a two-line induction, and both bounds the appendix needs fall out of
it plus two finite geometric sums: `|E S_g − μg| ≤ 2t` because `∑ 2^{-i} ≤ 2`,
and `Var(S_g) ≤ 3t²g` because each covariance row is damped by
`∑_{d≥1} 2^{-d} ≤ 1`, so the off-diagonal costs at most twice the diagonal.
No probability space, no operator algebra, no idempotent bookkeeping.

The variance statement is formalised as a pure finite-sum inequality with the
covariance *shape* as the hypothesis, which is slightly stronger than needed
(the inner sum runs over all `d < g` rather than `d < g-i`), and stronger is
free here since the extra terms are nonnegative.

### Checking the model, not just the algebra

The risk with collapsing to a scalar recurrence is that the collapse itself is
wrong. So the first numerical check builds an actual state space, an actual
`Π` that resamples uniformly over the nonzero states, applies `P = (I+Π)/2`
repeatedly, and compares against the closed form — exactly, over rationals.
It agrees. That validates the modelling step rather than the induction.

### Sharpness

`|E S_g − μg| ≤ 2t` is attained: the worst observed ratio to the bound is
exactly `1.0`, as `g → ∞` with `|m_0 − μ| = t`. So the `2` is the right
constant and not slack.

### Validation

200 state-space simulations for the modelling step; 20,000 recurrence and
mean-bound instances over rationals; 5,000 variance instances; 20,000
`G_Q`-antitone instances. No violations.

`scripts/check.sh`: PASS.

## Iteration 72 — T11b: the two-state coarse-graining

`R J = I₂`, the entries of `J D R`, and the Lipschitz step.

### The approximation statement needs no further estimate

The point worth extracting: `J D R` has its **live block constant** at `d/M`,
zero-to-zero equal to `1`, and both cross entries `0`. So
`eq:imt-empty-limit`, which bounds `|W₀^g(q,q') − e^{-θμg/L}/M|` uniformly
over live `q,q'`, *is* the entrywise distance to `J D R` — there is nothing to
prove between "close to the constant `d/M`" and "close to `J D R`". That is
what `entry_error_transfer` records, and the content is entirely in the entry
formula.

### Which direction is lossy

`R J = I₂` holds always. `J R`, the other composite, is *not* the identity on
the fine space — checked numerically, and it came out the identity only in the
degenerate `M = 1` cases (`s = 1`), never for `M > 1`. That asymmetry is the
reason the appendix can only claim an *approximation* of `W₀^g` by `J D_γ R`
rather than an identity: the coarse space genuinely forgets which live state
it is in, and the empty-epoch mixing is what makes that forgetting harmless.

### A `fin_cases` trap

`fin_cases i` on `i : Fin 2` leaves the goal mentioning `(fun i => i) ⟨0, ⋯⟩`
rather than `0`, so subsequent `rw`s with `0`/`1` literals do not match, and
the error is reported as a missing pattern rather than as a normalisation
problem. Case-splitting on the *if-conditions* instead — `by_cases hi : i = 0`
— keeps `i` abstract, makes the hypotheses usable by `simp`, and avoids the
issue entirely.

### Validation

`s = 1..7`: `R J = I₂` and the `J D R` entry formula exact over rationals with
random `d`; 200,000 random instances of the Lipschitz bound. No violations.

`scripts/check.sh`: PASS.

## Iteration 73 — T11c1: the continuum kernel `K_a`

The `a = 1` and `a = 2` integrands, the three scalar integrals, and `K_1^+`.

### Checking the paper, not a transcription of it

The `K_a^+` formulas are where a sign or a transposition would hide, so the
integrands are stated against **Mathlib's own matrix product** rather than
against a hand-expanded version. Working the `a = 2` case out by hand first
gave `(r d₁, d₁d₂; r d₀d₁, r d₀d₂ + d₀d₁d₂)`, which integrates to the paper's
`[[rg, e], [re, re + a₀]]` — so the closed forms are right, and that was
confirmed independently by quadrature over the triangle (agreement ~2·10⁻⁴ on
a 900×900 midpoint grid, consistent with the quadrature error alone).

### The fix that unstuck the matrix identities

Stating the `(2,2)` entry as `e^{-γ}` cost four failed rounds: `simp` kept
renormalising the exponent arguments into a form none of my prepared
identities matched, and the residual goal was *literally* one of the
hypotheses yet neither `assumption` nor a targeted rewrite fired.

Stating the entries as the **raw products the matrix multiplication actually
produces** — `e^{-γu₀}·e^{-γu₁}·e^{-γu₂}` rather than `e^{-γ}` — made every
remaining goal a pure ring identity in the exponential atoms, and it went
through immediately. The exponent arithmetic then lives in its own two-line
lemmas (`exp_prod_one`, `exp_prod_two`), where the simplex constraint is used
explicitly. This is the better factoring anyway: the matrix identity is now
independent of the constraint `∑ uᵢ = 1`, and `exp_prod_two` is the only place
that constraint is needed.

### Integrals without substitution or parts

All three scalar integrals are one `integral_eq_sub_of_hasDerivAt` each, with
an antiderivative written down directly — no change of variables and no
integration by parts. `∫₀¹ e^{-γ(1-u)}` is *not* done by substituting `v=1-u`;
its antiderivative `e^{-γ(1-u)}/γ` is just as easy to exhibit.

### Validation

20,000 random instances of each integrand identity; the three integrals
against numerical quadrature (max error 1·10⁻¹⁴); and the paper's `K_2^+`
closed form against 2-simplex quadrature at two parameter settings.

`scripts/check.sh`: PASS.

## Iteration 74 — T11c2: `K_2^+` over the 2-simplex

All four entries, matching the paper's `[[rg, e], [re, re + a₀]]`. T11c is
complete.

### No Fubini, no measure theory

The simplex integral is written as a **nested `intervalIntegral` with a
variable inner endpoint**,

    2 ∫_{u₀=0}^{1} ∫_{u₁=0}^{1-u₀} … du₁ du₀,   u₂ = 1 - u₀ - u₁,

which is a well-defined term needing no product measure and no Fubini
theorem. Each level is then one antiderivative: the inner integral is
evaluated for fixed `u₀`, and the result integrated again. Four
variable-endpoint integral lemmas carry the whole computation.

That the parametrisation and the `a! = 2` factor are the right reading of
`eq:imt-continuum` was checked numerically before formalising: nested
quadrature against the paper's closed forms agrees to `7·10⁻¹⁵`.

### The entries are not symmetric in the way one might guess

`(0,1)` and `(1,0)` both come out proportional to `e`, but for different
reasons: `(0,1)` because `e^{-γu₁}e^{-γu₂} = e^{-γ(1-u₀)}` makes the inner
integrand constant in `u₁`, and `(1,0)` because the inner integral of
`e^{-γ(u₀+u₁)}` telescopes against `e^{-γ}`. Only `(0,0)` produces `g` rather
than `e`, and it is the one entry whose inner integrand does not involve `u₀`
at all.

### Two normalisation traps

`field_simp` rewrote some occurrences of `-gam * u0` to `-(gam * u0)` but not
others, leaving a goal that `ring` could not close because the two forms are
different atoms. `simp only [neg_mul]` first, on both the goal and the
hypothesis, makes them agree.

`congr 2` on a goal of the form `(exp A - exp B)/γ = (exp C - exp B)/γ` does
not leave `A = C`; rewriting with `sub_sub_cancel` to fix `1 - (1 - u₀) = u₀`
before the arithmetic is both shorter and predictable.

`scripts/check.sh`: PASS.

## Iteration 75 — T11d1: the weighted row norm

The paper says outright what proves `eq:imt-product-norm`: "Submultiplicativity,
`‖D_σ‖_v ≤ 1`, and `‖(I+uP_+)D_σ‖_v = m_v(u)`." Those three, plus the
composition along the interleaved product, are what this iteration does.

### A predicate instead of a supremum

The weighted row norm is `max_i (1/v_i) ∑_j |A_ij| v_j`. Defining that maximum
would bring in a supremum over `Fin 2` and a division by `v₁`, neither of which
is ever used. Carrying the bound as a predicate instead,

    RowNormLe v₁ c A  :  ∑_j |A_ij| v_j ≤ c · v_i   for each row,

is exactly the form every step consumes, makes submultiplicativity a two-row
calculation, and keeps `v₁` out of every denominator except inside `m_v` where
the paper puts it.

Seen this way, `m_v(u) = max{1 + uσv₁, ur/v₁ + σ(1+u)}` is not an arbitrary
maximum: the two arguments are literally the two row requirements for
`(I+uP_+)D_σ`. Checked numerically with the paper's own constants
(`σ = 127/250`, `v₁ = 3/1600`, `r = 1/(2¹⁹-1)`), the bound is an **exact
equality** in all 20,000 cases — so the paper's `=` is right and `m_v` is the
norm, not merely an upper bound.

### What is not used

The paper notes "this reasoning requires nonnegative coefficients, not
stochasticity of `P_+`", and that is visible in the hypotheses: only
`0 ≤ u, r, σ, v₁` and `σ ≤ 1` appear, and they are used only to drop absolute
values. Nothing anywhere needs a row to sum to one.

One hypothesis was removed after checking it was unnecessary: `RowNormLe v₁ 1 I`
holds for *every* `v₁`, including negative, since the two obligations reduce to
`1 ≤ 1` and `v₁ ≤ v₁`.

### Validation

20,000 random rational matrices for submultiplicativity; 20,000 instances of
the two factor norms at the paper's constants; 3,000 interleaved products of
random length. No violations.

`scripts/check.sh`: PASS.

## Iteration 76 — T11d2: the Beta identity

`E[(1-Y+Y/σ)^{-(Q+1)}] = σ^l` for `Y ~ Beta(l, Q+1-l)`. T11 is complete.

### The Beta value is never needed

The identity looks like it requires the Beta function, but the same
`∫₀¹ u^a(1-u)^b du` sits on both sides and cancels. Stating it as

    ∫₀¹ y^a(1-y)^b / (1-y+y/σ)^{a+b+2} dy = σ^{a+1} ∫₀¹ u^a(1-u)^b du

(with `a = l-1`, `b = Q-l`, so no truncated subtraction appears) removes all
Beta and Gamma machinery from the proof — what is left is one change of
variables. Dividing through by the Beta value recovers the paper's
expectation form.

Under `y = uσ/D`, `D = 1-u(1-σ)`, the substitution is *exact*: the three
factors contribute `D^{-a}`, `D^{-b}` and `D^{a+b+2}`, the Jacobian
contributes `D^{-2}`, and the powers sum to zero. The transformed integrand is
literally `σ^{a+1}u^a(1-u)^b`, with nothing left over — which is what makes
the answer a clean power of `σ`.

### A `field_simp` failure mode worth naming

`field_simp` on an expression containing `(1/D)^n` produces terms like
`(D⁻¹)^a · D^a`, which `ring` cannot cancel because it does not know
`x⁻¹^a · x^a = 1` without `x ≠ 0` — and the resulting goal is a page long and
uninformative. Eliminating the reciprocal *first*, by rewriting
`(1/D)^n = (D^n)⁻¹` and then `X / Y⁻¹ = X · Y`, leaves only ordinary
denominators, after which `field_simp; ring` closes it immediately.

### A validation scare that was not one

The first numerical check reported a maximum relative error of `1.1·10⁻³`,
which for an exact identity is alarming. Splitting by `σ` showed it is
entirely quadrature: `7.7·10⁻¹⁶` for `σ ∈ [0.5,1]`, `1.6·10⁻¹⁵` for
`σ ∈ [0.2,0.5]`, degrading only as `σ → 0` where `(1-y+y/σ)^{-(a+b+2)}` becomes
sharply peaked at the origin. The paper's `σ = 127/250 ≈ 0.508` is in the
well-conditioned range. Worth recording because "the identity is wrong" and
"the quadrature is wrong" look identical until you stratify.

`scripts/check.sh`: PASS.

## Iteration 77 — T12: the scope call, and the end of the queue

T12 was flagged "may be out of scope — flag before starting". It is two
separate claims from the theorem statement, and they get different answers.

### The rate half is already done

"Every realized code has rate `1/2`" and, through the requested-length
wrapper, `1/2 - o(1)`. `Wrapper.lean` proves the wrapper is injective and
weight-preserving; `Schedule.lean`'s `wrapper_rate` proves
`N_{m(n)}/n → 1`. Together those are the `1/2 - o(1)`. What is *not* modelled
is the native rate being exactly `1/2`, which is a dimension count of the
encoder; this development never fixes the encoder's dimension because the
distance argument never uses it. That is a small, clearly located gap, not a
hidden one.

### The linear-work half is declined

`O(N)` bit operations needs a **cost model for bit operations**, and neither
Lean nor Mathlib has one. If I define one, three things follow:

*  the theorem becomes a statement about *my* model, not about the paper's
   claim;
*  the mathematics is trivial — `N/128` steps times a constant — so the whole
   content sits in whether the model faithfully represents bit operations,
   which is exactly the part a Lean proof cannot establish;
*  a model that charges `1` per step makes the claim true by construction.

That last point is decisive. This project has hit vacuous measurements twice
already (the `Option.isSome` timing that the kernel satisfied without
evaluating the payload), and the pin discipline exists to catch precisely this
shape of error. Producing a green `example` in `Pin.lean` asserting `O(N)`
work would look like certification and would certify nothing. Writing three
lines of trivial arithmetic dressed as the result would be worse than leaving
it out, because it would appear in the pin file as if the claim were covered.

Nothing in `distance_whp` — the theorem's mathematical content, and the whole
subject of this formalization — depends on either claim.

### Queue empty

`scripts/check.sh`: PASS. Axiom closure of the headline theorem is still
exactly `[propext, Classical.choice, Quot.sound]`, with zero `sorry`, zero
`admit`, zero project axioms, and no `native_decide`.

## 2026-09-24 — T7c re-measurement (out of loop; user-directed)

Not a loop iteration: the queue is empty and this was a direct request to
resolve T7c's blockage. No Lean files touched, so no `check.sh` run.

**Method change that made the difference.** Node counts say nothing about how
much of a search remains. Added `scripts/majorant_progress.py`, which tracks
*resolved volume* — the summed volume of closed and pruned boxes. For a DFS
that is exact progress, and elapsed/fraction is a real ETA. Every earlier
estimate in D4 came from extrapolating node counts and was wrong.

**Result: segment 0's search does not terminate.** Resolved volume stalls at
32% while the ETA climbs monotonically (0.16h → 0.37h over 60k nodes).
Extrapolation gives ~`1e16` boxes. Two hypotheses tested and killed:

* feasibility-tightened start — big throughput gain, no convergence gain;
* `MAJ_MAXDEPTH=120` — **no effect at all**; observed depth pins at 56 under
  both caps, so the cap was never binding and its closeness to the
  certificate's `deepest_box = 55` was coincidence.

Also ruled out earlier: tightness (segment 0 has the *most* slack) and `Fix`
precision (enclosure converges to ~`8e-13`).

**What is left.** The per-box bound. ~`1e16` boxes versus 101,940 for the
published cover is eleven orders of magnitude, which no split heuristic
reaches — the naive interval extension loses to dependency (`c` appears in two
`π` terms, `a` in two places). A centred/mean-value form or monotonicity in
`x` is the untested lead.

**Artifact finding.** The certificate is 3,702 bytes and stores only summary
statistics — slopes, intercepts, box counts, `status: proved`. The boxes are
not in it and neither is the search procedure, so it **cannot be re-checked
from its own contents** and there is no cover to transcribe. This is the gap
the Lean replay exists to close.

**Corrections made.** `DECISIONS.md` D4 carried three errors, all inflating
cost in the same direction (nineteen supports vs fifteen; a unit-confused
comparison of the certificate's whole-claim 101,940 boxes against a
single-segment 32,054; and a near-tangency explanation contradicted by the
measured margins). Corrected in place with the evidence, and `TASKS.md` T7c
restated: the blockage is search, not compute.

New, uncommitted: `scripts/majorant_progress.py`, `scripts/majorant_split.py`,
`scripts/majorant_floor.py`.


## 2026-09-25 — T7c closed with the original three-dimensional argument

User-directed continuation of the stalled concave-majorant formalization.

The original verifier's `Box.taylor_upper` was missing from the Lean search.
Added `Numeric/MajorantDefs.lean` and `Numeric/MajorantEval.lean`: a centered
bound for the complete objective, keeping derivative cancellations before
multiplication by the coordinate widths. Soundness reuses the existing
one-dimensional mean-value expansions and fixed-point interval lemmas.

Added `Numeric/BAClosed.lean` to identify the entropy and division-free
accumulator exponents on the **closed** feasibility region, and strengthened
the new majorant evaluator interface to include those boundaries. Existing
pinned statements and old evaluator interfaces were not changed.
`Numeric/BACentral.lean` proves the central support via entropy concavity,
the generating-function witness u=1, and a checked rational rounding of log 2.

The current paper uses `d11/OUTER_REFINED.json`, not the older artifact targeted
by the stalled probe. Generated all 19 nonconstant left-support covers with
outward-rounded endpoints. Search: 108,045 nodes, 54,032 leaves, maximum depth
56, approximately 22.01 seconds. The original absolute-width split rule and
memoized exact log evaluations are used. No floating-point result is trusted.

Verification already completed:

- all 387 independent numerical modules passed kernel checking;
- all 19 segment assemblies passed;
- exact support/reflection/active-interval proofs in `Majorant/RefinedData.lean` passed;
- `Majorant/Refined.lean` compiled, proving `baExponent_le_refined` on
  [13/125,112/125], with all feasibility boundaries included;
- `lake build SpinCodes`: passed (8989 jobs, mostly cached);
- `bash scripts/check.sh`: ALL CHECKS PASS, original pins unchanged,
  no forbidden proof escapes, headline closure the three standard axioms.

The numerical replay reports, source hashes, and full witnesses are retained
in `scripts/majorant_data/`. `scripts/check-majorant.sh` replays the numerical
modules, assemblies, final theorem, new pins, and axiom audit. Final new-pin
and axiom-check results are appended below when complete.

D5 corrects the earlier unsupported nontermination and missing-source claims.
`CLOSURE_AUDIT.md` records that the concrete distance theorem still requires
T3a, the actual family instance, and the remaining asymptotic/occupation
connections. Closing this majorant does not silently discharge those inputs.

**Final audit — 2026-09-25 05:14 America/Los_Angeles: PASS.**
`Majorant/Pin.lean` compiled. `baExponent_le_refined`, the combined evaluator,
the closed-boundary identity, and the analytic central bound close over
`[propext, Classical.choice, Quot.sound]`. `S0.P0.logs_ok` and `S0.P0.tree_ok`
depend on no axioms. The new final theorem and exact support list are pinned.
The large replay phases were run separately with bounded concurrency; the
complete reproduction sequence is now `scripts/check-majorant.sh`.

## 2026-09-25 — T3a: complete polynomial program verified

The full concrete distance theorem remains the goal. This iteration closes
the coefficient-identity and normalization gap for the sparse program, but
does not mark T3a or the overall theorem complete.

Added the actual fixed-weight matrix formulas, their binomial mixture,
column-action identities, and nonnegativity and monotonicity lemmas in
`Occupation.lean`. `SparseWitness.lean` proves the exact affine witness's
weighted average and uniform coordinate lower bound. `SparseModel.lean`
defines the numerical max/min matrix and identifies its zero-coordinate
action with the program.

The pure polynomial program follows the original verifier's finite branch
choices. Lean now checks 1,032 identities across 129 weight modules, all 693
maximum comparisons, all 119 contribution sums, and all seven final residual
identities against the original integer rows. The recovered positive
denominators and initial factors of `x` are included in the proof.

Direct convolution was too expensive for the large sums. `PolyPacked.lean`
proves exact integer encoding sound: the radix exceeds the sum of absolute
coefficients of the difference, so an encoded zero cannot hide a carry.
The packed checker verifies all scaling denominators and uses no compiled
evaluation. Explicit power tables avoid repeated expansion of 128th powers.

`SparseMaximumAll.lean` proves the real maximum inequalities for all 129
input weights. `SparseContributionSound.lean` and `SparseResidualsSound.lean`
connect all arithmetic pieces. `SparseProgram.lean` verifies exact coverage
of weights 0 through 128 and proves `programResidual_neg` for every `Fin 7`
coordinate on `0 < x ≤ 1`, with contraction `1 - 96*(x/10000)`.

Verification completed:

- 129/129 weight modules, 129/129 dominance modules, 119/119 sum modules;
- seven residual identities and all real assembly files;
- `SparsePin.lean` and `SparseBridge/Pin.lean`; original pins unchanged;
- `lake build SpinCodes`: PASS, 9,005 jobs;
- `bash scripts/check.sh`: ALL CHECKS PASS;
- final theorem and packed-checker soundness: the three standard axioms;
- sampled packed numerical checks: no axioms; other sampled identities use
  at most `propext`.

The first weight run had 27 missing-import failures when a concurrent Lake
build temporarily removed dependency objects. All 27 passed after the build
finished. Both logs are retained; `kernel_weights_complete.json` combines
the successful results. The replay sequence now builds dependencies first,
and the numerical runner records and checks dependency hashes.

Reports and source hashes are in `scripts/sparse_data/final_verification.json`.
Replay with `scripts/check-sparse.sh`. `SPARSE_BRIDGE.md` describes the exact
remaining work: selected cancellation bounds, evaluation of the six live
actions, and the binomial-mixture assembly. Then continue with the concrete
family and the outstanding asymptotic/occupation connections in
`CLOSURE_AUDIT.md`. No commits were made.

## 2026-09-25: T3a closed for the numerical occupation matrix

`SparseCancellation.lean` proves positive denominators and the selected
cancellation/fresh-state bounds from the checked finite indices. It also
identifies the integer subtraction in the polynomial cap with the natural
subtraction in the real shell formula using the checked kernel bound.

`SparseColumnBounds.lean` evaluates all six live polynomial actions and proves
their fixed-weight column inequalities. The zero-live-mass branch is matched
exactly with `kernel = choose(128,j)`, including both endpoint weights.

`SparseContraction.lean` assembles the binomial mixture and all seven program
residual inequalities. `Spin.Imt.Occupation.Sparse.sparse_collatz` now states
the paper's inequality directly at `(β,z)=((4/5)α,1-(8/5)α)` throughout
`0 < α ≤ 1/10000`. No polynomial identity, sign, or selection inequality is
left as a hypothesis of this theorem.

`SparseNonneg.lean` and `SparseIteration.lean` establish matrix nonnegativity
and `e_Z T^R 1 ≤ 2048 (1-96α)^R` for every natural `R`. The supplemental pins
check both the contraction and iteration statements; the original pins are
unchanged.

Verification: all five new semantic modules and the supplemental pin passed
a fresh sequential kernel replay. Both final theorems depend only on
`propext`, `Classical.choice`, and `Quot.sound`. `lake build SpinCodes` passed
(9,005 jobs), and the actual `bash scripts/check.sh` reports `ALL CHECKS PASS`.
The existing numerical certificate sources are unchanged; their prior
verified module coverage is retained. The new semantic replay records source
hashes and individual logs in `scripts/sparse_data/semantic_verification.json`.

T3a is DONE. T1a is OPEN: next define the paper's exact maps, identify their
spectra and fiber data, and instantiate the one-step domination law, then the
concrete family and its remaining selection/occupation estimates. The full
distance goal remains open. No commits were made.

## 2026-09-25: concrete map structure proved; full spectrum replay started

`PackedMapDefs.lean` gives the least-significant-bit-first XOR evaluator.
`PackedMap.lean` proves its XOR law, composition and inverse rules, output
bounds, one-hot evaluation, support bijection, and Hamming-weight/cardinality
identity. The selected maps use the exact 19 rows in the paper's map table.

`ConcreteMapData.lean` checks 13 finite claims: dimensions and row bounds,
explicit inverse basis identities, nonzero expansion columns, and feedback
column weights/distinctness. `ConcreteMaps.lean` lifts the checked inverse
identities to injectivity of A and C transpose and surjectivity of C for all
inputs, and transfers these results to the finite-support representation.
It also proves weight five, nonzero value, and distinctness of the actual
feedback images `C (inputBasis i)`. New supplemental pins preserve these
statements; the original pins remain unchanged.

All nine structural modules passed a fresh sequential replay. The semantic
results use only the three standard axioms; the sampled finite checks use
only `propext`. The actual `bash scripts/check.sh` reports `ALL CHECKS PASS`.
`scripts/verify_map_checkpoint.py` passes and records the scope and source
hashes in `scripts/map_data/checkpoint.json`.

The spectrum generator emits 512 blocks of 1,024 states for both A and C
transpose. Direct checks passed but repeated too much work. The final
version shares a kernel-checked low-ten-bit image table, uses one XOR with a
checked high image, and proves its faster weight recurrence equivalent to
the original. `PackedMap.eval_split` and `MapSpectrum.block_split_weights`
connect this computation to the exact evaluator. Generic histogram and
range-partition lemmas ensure every state is counted once.

The basis replay passed in about 10 seconds; the nonzero-high-part sample
block 257 passed in about 61 seconds. The full three-worker replay is now
running under `scripts/check-map-spectrum.py --jobs 3`. Preserve and poll its
existing session; do not start a duplicate or rebuild its numerical imports.
The full report is `scripts/map_data/kernel_blocks.json`. Only PASS proves a
completed replay. `map_spectrum_assemble.py` refuses to emit the total and
cardinality bridge until all 512 current source hashes have passed.

Next: finish the spectrum replay/assembly, prove the C/C-transpose pairing,
derive the frozen fiber and cancellation data, and instantiate one-step
domination. The concrete family and asymptotic estimates are still required
for the full distance theorem. T1a remains in progress; no commits were made.

## 2026-09-25 — Concrete character pairing and fiber bounds

The actual `Cset` and `CtransposeSet` now satisfy the character-pairing
identity. The proof derives transposition on basis vectors from the packed
table and extends it by the established XOR laws. Consequently the abstract
orthogonal code is the actual kernel, and its fibers are the actual syndrome
fibers. Fourier inversion and Parseval are instantiated with no assumed map
or rank properties. A grouping theorem connects their support-weight counts
to the packed counts produced by the spectrum checker.

The actual feedback kernel has minimum distance at least four: its words
have even weight because every column has weight five, and distinct columns
exclude weight two. This supplies both packing bounds. The variance theorem
is instantiated over all nonzero syndromes, including empty fibers; the
number of terms is exactly `2^19-1`. Integer certificate rules now justify
each of the paper's four fiber caps, including the quadratic test at `cap+1`.
The fresh `check-map-fibers.py` replay passed all eight modules, including
the supplemental pins, and all 16 axiom audits. The proof modules use only
standard axioms; sampled numerical checks use only `propext`. Source and
output-object hashes and logs are recorded in `fiber_verification.json`.

The spectrum replay was deliberately restarted after profiling Lean's
`decide +kernel` mode. It avoids the duplicated elaborator computation and
uses only `propext` in the sampled numerical checks. The earlier partial
report is preserved as `kernel_blocks_plain_decide.json`; no old source
hash is counted as a pass for a regenerated module. The new replay records
output-object hashes and supports hash-checked resume. The full replay
remains active; preserve it and its dependencies across turns.

An explicitly untrusted arithmetic preview reproduces all 129 frozen kernel
counts from the proposed spectrum and finds witnesses for all 129 caps:
125 Fourier bounds, two packing bounds, and two endpoint complement bounds.
These witnesses are saved in `fiber_candidates.json`; they still need Lean
checking and the completed concrete spectrum identification.

Next: finish and assemble the current spectrum replay, then kernel-check
the proposed Fourier sums and cap witnesses. Use the new weight-count bridge
to derive the frozen kernel counts and caps for the actual map. Low-input
cancellation counts, shells, one-step domination, the concrete family, and
the asymptotic estimates remain. The original pins are unchanged; no commits
were made. The full distance theorem remains open.

The earlier structural checkpoint and its invariant-check log remain valid
for their recorded scope. This continuation uses direct Lean checks for the
new modules; it does not rerun Lake while the numerical workers depend on
the current object files. Repeat the full project check after the workers
and spectrum assembly finish.

## 2026-09-28 — concrete counts and memory-bounded assembly

The spectrum replay finished with 512/512 blocks passing; the histogram
and semantic assemblies also compiled. The Fourier replay finished with
129/129 modules passing. No numerical rejection occurred in the six spectrum
blocks retried after timeouts.

The first combined count assembly exhausted practical memory. Abstract
finite elimination and separate per-weight cap lemmas isolated a kernel
memory failure at the weight-128 complement bound. Moving its arithmetic
rewrite into a generic lemma before specializing the weight resolved it.
The revised `FiberNumericsAll` compiled in 39.73 seconds with one worker and
a 10,000 MB Lean memory limit. `ConcreteCounts` compiled in 28.59 seconds.
The final eleven axiom audits passed using only standard axioms. The report
is `scripts/map_data/concrete_counts_verification.json`; numerical dependencies
were reused from their completed, hash-recorded replays.

`ConcreteCounts` now identifies the actual transpose spectrum, kernel and
pair counts, nonzero-syndrome fiber caps, and five actual expansion-map shell
sizes with the frozen numerical data. This discharges the former spectrum
premise. The corrected live/nonnegative induction and its separate connection
audit also passed. Original statement pins were not changed.

The user authorized Peach for additional compute. Its isolated workspace is
`/tmp/spin-lean-peach/project`, with the matching Lean 4.34.0 toolchain beside
it. The transferred dependency manifest records 691 locally checked project
modules. Remote checks reuse those modules; they are not an independent full
numerical replay. Remote setup/check results are recorded separately.

Next: prove the low-input cancellation facts and the actual encoder's
one-step `LiveDominates` bound. Then instantiate the concrete random family
and discharge the selection/occupation asymptotics. The unconditional
distance theorem remains open. The default project build and `scripts/check.sh`
passed after the final count audit (`counts_invariant_check.log`). The new
count modules remain supplemental direct checks rather than default imports.

The Peach assembly check subsequently passed with the same source hashes:
6.00 seconds for `FiberNumericsAll`, 4.00 seconds for `ConcreteCounts`, and
3.75 seconds for the eleven-closure audit, with peak resident memory below
7 GiB. The local report and logs are `peach_counts_verification.json` and
the associated `peach_*.log` files. Windows and Linux produced identical
object hashes for both final modules. No verification job remains running.

## 2026-09-28 — concrete one-step law for the general input weights

`ConcreteTransfer.fixedStep_liveDominates_generic` now proves the actual
fixed-weight transfer law for every `j : Fin 129` except 1 and 2, for
`0 ≤ z ≤ 1`. It applies to arbitrary measures satisfying `LiveDominates`
and arbitrary nonnegative coordinate budgets. It uses the existing numerical
matrix and has no assumed row-domination or one-step premise.

The new modules identify actual emitted moments with the hypergeometric
formula; prove pointwise and shell-average cancellation bounds from the
concrete fiber data; decompose the actual transvection/input transition;
prove the exact zero row and refresh zero-column minimum; handle the
all-kernel shell-preserving branch; and lift the resulting row bounds
through a nonnegative kernel. The zero row is proved at all 129 weights.
`fixedStep_liveDominates_of_cancellation` isolates the cancellation bounds
needed to extend the result to the two special weights.

The fresh sequential Peach replay passed eleven semantic/pin modules and
seventeen axiom audits. Only standard axioms occur; the small check that
special patterns appear only at weights 1 and 2 uses no added axiom.
The report is `scripts/map_data/concrete_step_verification.json`, with
source/object hashes and per-module logs. Existing numerical dependencies
were reused. `scripts/peach-lean.py` transfers changed inputs through the
host-key-pinned SSH connection, compiles the requested module, and fetches
its proof objects. No simultaneous benchmarks were run.

Next: enumerate the 128 singleton inputs and 8,128 unordered pairs using
the packed maps. Check their feedback syndromes, expansion weights, and
emitted weights against `weight1` and `weight2` in `SparseModelData.lean`.
The diffuse bound needs the syndrome-grouped output-weight patterns; the
shell bound needs the output-weight histogram grouped by the expansion
weight of the feedback syndrome. A finite enumeration certificate plus
semantic reindexing should discharge both remaining cases. Then form the
Bernoulli mixture and apply the sparse iteration result. The concrete
family and the final selection/occupation asymptotics remain open.

The original statement pins and the paper were not modified. The 127/129
count is a count of input weights in the new one-step theorem, not a
completion percentage for the full distance theorem. No commits were made.

Final verification for this continuation passed: the fetched proof objects
passed the local supplemental pin/axiom check, and the default build plus
`scripts/check.sh` reported `ALL CHECKS PASS`. All new source/object hashes
still match the semantic replay. See `step_final_verification.json`.
No verification job remains running.

## 2026-09-28 — low cancellation tables and actual occupation transfer

Closed both exceptional fixed-input weights. The new generator groups all
128 singleton inputs and 8,128 pairs by their concrete feedback syndromes.
Lean checks each packed input, syndrome, expansion weight, and emitted
weight. It proves complete layer coverage from cardinality and absence of
duplicates, checks distinct syndrome keys, and checks all five shell sums.
The pair layer has 7,909 occupied syndromes. All 63 finite blocks passed.

`LowCancellationBridge`, `LowCancellationMoments`, and
`LowCancellationShells` connect these certificates to actual fiber moments.
`ConcreteLowCancellation` discharges both special numerical envelopes.
`fixedStep_liveDominates` in `ConcreteTransferAll.lean` now covers all 129 weights.
The earlier generic theorem and original statement pins were preserved.

Kernel reduction of the library merge sort stopped at an opaque
accessibility proof. `LowCancellationSort` instead uses explicit list
recursors and a bounded structural sort. A theorem proves multiset
preservation; each certificate separately checks the needed sortedness.
The generator's first aggregate proof also needed conjunction reassociation.
After correcting that proof assembly, the weight-2 aggregate passed in
171.62 seconds including transfer, under the 20,000 MB Lean limit.

`ConcreteBernoulli` defines the direct transition under independent input
bits and proves equality with the binomial mixture. `ConcreteOccupation`
then proves the actual one-step domination, the finite-iterate matrix bound,
and `2048 (1-96α)^R` at the sparse parameters. No cancellation-table or
one-step premise remains in those concrete statements.

The sequential certificate replay and assembly checked 76 modules in total.
All 16 supplemental axiom audits passed with only the standard axioms.
Reports and hashes are in `scripts/map_data/low_cancellation_verification.json`;
the local follow-up is recorded in `low_final_verification.json`.
The earlier spectrum, Fourier/count, and sparse-polynomial dependencies
were reused. No simultaneous benchmarks were run.

The full distance theorem remains open. Next, define the concrete finite
encoder experiment and identify its emitted-weight moment with the kernel
iterate. Discharge the routing-law model assumptions, instantiate the shared
outer/route/inner family and actual threshold, and apply the selection and
occupation estimates to its conditional first moments. See
`LOW_CANCELLATION.md` and `CLOSURE_AUDIT.md`.

Final verification passed locally: the supplemental import/axiom audit took
208.23 seconds, and the default build plus `scripts/check.sh` took 101.67
seconds and reported `ALL CHECKS PASS`. All 76 source/object hashes match
the successful replay, and the original distance pins retain their hash.
No verification job remains running. No commits or paper changes were made.

## 2026-09-28 — concrete routed IMT moment, parallel proof assembly

The user authorized three parallel subagents. Their separate file ownership
covered the finite encoder, actual route law, and permutation/row laws.
The root supplied finite-law composition, reshaping, region-major
serialization, final integration, consolidated verification, and documentation.
The explicit goal was to close this finite routed moment bridge; it did not
claim to close the full asymptotic distance theorem.

`ConcreteEncoder` defines the actual 128-bit/19-bit recurrence using the
checked `Aset` and `Cset`, independent uniform valid transvection pairs,
empty initial state, and no terminal tail. `ConcreteEncoderMoment` proves
that its iid Bernoulli generating function equals the iterated direct
Bernoulli kernel, and hence obeys the numerical matrix and sparse bounds.

`ConcreteShufflePermutation` identifies the layer law with the pushforward
of actual uniform permutations. `ConcreteRoutePermutation` supplies the
independent row/region seed experiment and discharges the identification
premises in the route inequality. The cost is exactly
`(b+1)^Q (L+1)^b`, with cost one for empty rows.

`ConcreteReshape` proves that every fixed position bijection preserves iid
bits. `ConcreteSerialization` pins the region-major index relation.
`ConcreteRoutedEncoder.stream_moment_bound` combines the route and encoder
over an explicit product sample space. It needs only positive dimensions,
divisibility by 128, interior total input weight, and `0 ≤ z ≤ 1`.
Independent semantic review confirmed the recurrence, serialization, and
independence structure against the paper.

The sparse specialization retains the full route factor and assumes actual
input density `(4/5)α`. It does not close sparse occupation for arbitrary
outer profiles. That requires the fair-row counting comparison and
marked-position conditioning; the generic `(L+1)^b` loss is too costly.

The new consolidated script is `scripts/check-routed-bridge.py`. Its report
`scripts/map_data/routed_bridge_verification.json` records all 19 replayed
modules, 26 standard-axiom audits, the local statement check, and the default
project invariants. It reuses earlier numerical dependencies and checks the
previous 76-module checkpoint hashes. Only `status: PASS` denotes completion.
No simultaneous benchmarks, commits, or paper changes were made.

Next recommended parallel milestone: concrete outer-family/selection
assembly, sparse counting-and-conditioning, and fixed-occupation kernels.
Then combine their first moments with the positive-occupation argument.
See `ROUTED_MOMENT.md`, `STATUS.md`, and `CLOSURE_AUDIT.md`.

Final verification: **PASS**. The 19-module sequential Peach replay took
189.98 seconds, the Windows statement/26-axiom check took 282.73 seconds,
and the default build plus invariant script took 140.00 seconds and reported
`ALL CHECKS PASS`. Source/object hashes match, the original distance pin
retains its hash, and the checked paper sources are unchanged. The finite
routed-moment goal is complete. No verification job remains running.


## 2026-09-28: continuing parallel full-proof goal

The user requested an ongoing delegation goal. Three workers continue in
parallel while the coordinator integrates and checks; the runtime allows
four total agents. Work is replenished after each bounded result.

The stable finite checkpoint passed: 44 modules replayed on Peach in
624.84 seconds; 90 standard-only axiom audits; Windows combined pin in
266.66 seconds; default invariants in 99.17 seconds. The preceding 95
module records and original pins match, and paper files are unchanged.
The report is `scripts/map_data/parallel_finite_verification.json`.

Subsequent native sparse assembly passed seven modules and 19 audits.
Its actual binary-message sum transports exactly to the outer tuple sum;
the selected b² cost is explicit; the sparse exponential rate is linked
to `concreteFamily.EZ`. The two outer asymptotic facts (tail and shell
ratio envelope) are still hypotheses. The canonical family totalizes
conditioning by using the whole space only when raw Good has zero mass;
the actual bad-distance event is unaffected. Eventual agreement follows
from a vanishing raw selection failure probability, not an assumption
that raw Good is positive at every finite index.

Workers also closed the parity/boundary entropy bridge to the certified
refined majorant, exact impulse composition, and the coarse-product error.
The dense replay contains 1,023 boxes, but no numerical dense-box proof
has yet been completed in Lean. The rational occupation checker is sound;
a fixed-point evaluator is being built for efficient kernel verification.
The coordinator is proving the scalar and Fourier transfer alternatives.

Next recommended work: certify the scalar outer defect and tail, finish
actual placement/Riemann limits for fixed occupation, and certify all dense
witness families. Then discharge the remaining asymptotic interfaces. The
full proof goal remains ACTIVE; this finite checkpoint is not closure.

## 2026-09-28 — actual distance semantics and mixed dense integration

The outer tail, native selection failure and growing sparse regime are now
unconditional. Actual outer/inner XOR laws, emitted word weight and injection,
linear image dimension N/2, rate 1/2, and minimum-distance event equality all
pass their dedicated kernel checks and standard-only audits.

Global geometry covers all 1,023 dense source boxes. The common PointRate
interface supplies selected profile/layer bounds for actual seeds; root
proved native prefactors have a vanishing remainder and connected DenseRates
to the exact actual minimum-distance theorem. Numeric replay is in progress
for all three families. More than 100 scalar local boxes have passed; one
transport/cache write failure was retried successfully and awaits manifest
reconciliation after the batch. The first exact Fourier witness passes with
radius approximately 2.87e-32. Indexed wrappers separately verify geometry
association and actual certified outer-support membership.

Fixed occupation now has uniform actual simplex convergence, robust margin
transfer through arbitrary region products, and row-profile coefficient
extraction. Actual continuum contraction bounds and their finite numeric
margins remain open; Q=1 explicit integration is underway.

Next recommended work: finish and combine all indexed dense certificates,
and discharge Q=1/Q=2 and larger-Q fixed continuum bounds, followed by actual
fixed EZ limits and final unconditional theorem replay. Goal remains ACTIVE.

The actual Q=1 native first moment is now closed unconditionally:
ConcreteNativeFixedOneLimit.native_one_EZ_tendsto. Its eventual bound is
1000 exp(-b/50). Exact region-major serialization, fair routed probability,
selected outer tuple counting and native schedule are all linked. Dedicated
report native_fixed_one_verification.json records four modules and fourteen
standard-only audits. Q=2 native integration assigned to route_law; imt_encoder
now owns fixed subset insertion and shuffle_laws owns larger-Q simplex bounds.
Scalar indexed wrappers resumed from hashes with three workers after a
controlled coordinator restart; no simultaneous compilation of a module.

The full subset insertion identity and all 62 larger-fixed numerical endpoint
inequalities now pass. Exact actual routing for every fixed row profile passes
in ConcreteFixedProfileRouted; its continuum norm premise remains explicit.
fixed_routing_verification.json records six modules and eleven standard-only
audits. Canonical order-statistic integration now passes via exact scaled
simplex volumes and Fubini, supplying the last spacing-law analytic ingredient.
All 283 scalar local boxes pass, including the recorded successful B028 retry;
566 standard-only audits. Indexed geometry wrappers are finishing. The three
workers now assemble actual Laplace, matrix norm, and native profile limits.
Next: close those large-Q bounds, finish indexed dense replay, instantiate the
actual minimum-distance theorem and check its full pin. Goal remains ACTIVE.

All fixed occupations CLOSED unconditionally. ConcreteNativeFixedLargeFinal
proves EZ(m,Q)->0 for every Q>=3 with eventual bound (1600/3)exp(-Qb/2000).
ConcreteNativeFixedAll combines Q=1,2,and>=3 and its fixed-range theorem PASS.
All283 scalar indexed geometry wrappers PASS; all105 Fourier matrix witnesses
PASS. Fourier indexed boxes run in four disjoint lanes; occupation matrices
in four workers, with separate local and indexed wrappers. Native semantic
review found no construction/rate/distance mismatch; paper complexity is a
separate unformalized clause. Final theorem and pin sources are prepared but
unverified until dense mixed certificate finishes. Root final checker session
92570 waits for that aggregate; then compiles both sources and records their
full project dependency graph and standard-axiom audits. Goal remains ACTIVE.
Next recommended work: finish indexed occupation/Fourier replay, check mixed
aggregate, resolve any final assembly errors, run invariant checks and finish
reproducibility/semantic reviews before declaring closure.

Default scripts/check.sh PASSED after current proof additions: default build,
original statement pin, no holes/new project axioms/compiled-trust certificates,
headline axiom closure and numerical certificate axiom audit. Original pin SHA256
remains447d959bec1591dd678a9696262ceb7a3df818b59cdda38086c2feb259f7ceca.
Scalar aggregate DenseOccupationScalarCertified PASS. Independent native and
dense semantic reviews found no distance/rate mismatch. Reproduction dry-run
caught two unchecked occupation wrapper substitutions (B152 geometry/B372 name);
generator fixed, affected wrappers regenerated, all433identitychecks PASS.
Occupation pipeline now matrices4(session25477), localboxes2(48758), indexed2
(13227); Fourier4disjointlanes. Mixed aggregate watcher35873 and root final
checker92570 remain active. Agent imt has checked generic success/existence
corollaries and will attach unconditional wrappers after final theorem PASS.
Reproducibility docs/script are ready; a fullsource replay has NOT been run.

Fourier family and aggregate COMPLETE: 105 witnesses, 307 boxes, 824 paired
source/object records, 1543 standard-only audits. Scalar aggregate also PASS.
Occupation pipeline safely drained and resumed with matrices5 (22634),
local3 (31043), indexed3 (5404), and two shared remote semantic permits.
Matrix processes measured about10.08GiB peak each. Ready boxes are scheduled
in witness order. Aggregate31234 and mixed47389 wait on audited PASS/hash
records. Root final theorem waiter15735 and consequence waiter94013 are live;
old handles in preceding entries are historical. Original helpers/stable
checked proof sources remain unchanged. Thirteen outer semantic modules
freshly rechecked in isolated Windows outputs, all byte-identical to live
objects, all standard-only axioms. Pending final theorem namespace preflight
fixed missing ConcreteOuter opens. Goal remains ACTIVE; next work is the
occupation certificate, final theorem/pins/corollaries, and final evidence
and invariant checks. No full-source replay claim is made.

Next goal continuation: previous turn classified PROGRESS. All five occupation
controllers and both final waiters were confirmed live by Windows process
inspection; root waiter15735 also polled successfully. Matrix witnesses
advanced to87/133, local220/433, indexed217/433. Recent matrix throughput is
about1.1 witnesses/minute; actual admitted-run minimum headroom36.60GiB.
Final evidence checker now enforces14 gates, exact module/audit identities,
exception-safe INCOMPLETE publication, and final default invariants. Read-only
run checked2515 current pairs with no hash conflicts or nonstandard axioms;
only the expected occupation/aggregate/final gates remain incomplete.
Root compared assembly import discovery with the independent nested-comment
parser: both find exactly5215 project modules, identical sets. Result in
native_import_graph_crosscheck.json. Final invariant runner prepared and
independently reviewed, not executed before its actual theorem prerequisites.
NATIVE_RESULT.md now documents exact statement/imports and explicitly PENDING
status. A fixed-occupation semantic review is assigned independently.
Twelve old foundations are being recompiled in isolated outputs. Fresh source
checks pass, but some bytes differ from live objects; no false pairing or live
replacement has occurred. Root found the old files include package metadata
spincodes and are96 bytes larger, consistent with original Lake --setup versus
direct Lean invocation. Owner is testing that explanation and comparing full
declarations before claiming source/object identity. Goal remains ACTIVE.

Continuation classified PROGRESS plus verified waits: occupation advanced to
98/133 witnesses and274/433 indexed boxes; controller processes and root
waiter15735 confirmed live. All12 old foundation sources and12 standard-only
audits PASS in isolation; live objects preserved. Declaration/package-metadata
comparison is ongoing and does not yet establish byte identity. Fixed-native
semantic review completed with no blocking mismatch, recording conservative
tilts/margins and the valid finite4095-occupation sum. Complete occupation
evidence watcher was hardened to exact1566 module/2132 audit identity sets,
immutable pre/post hashes, and FAIL publication on exceptions; root inspected
the revision. Consolidated completion now requires15 gates. Route agent waits
on live final controller for independent full-statement inspection; imt owns
artifact comparison and corollaries; shuffle owns remaining certificates and
aggregates. Next: finish numeric replay, actual assembly, corollaries, final
invariants, and inspect all completion evidence. Goal remains ACTIVE.

Next continuation PROGRESS:105/133 occupation witnesses and300/433 indexed
boxes PASS, no numeric failures. Root final waiter15735 and all occupation
controllers confirmed live. Framework live/fresh comparison PASS for23 named
declarations: full printed types, bodies, and axiom lists identical (348253
bytes; shared log SHA256fe3e431ce141ee820a9fcf69943ef9af3c94f12f7698df5d568e56aea4e3d8f3).
Local package-metadata rebuild stopped at the3GiB available-RAM guard after
unrelated user Chrome memory grew; live artifacts unchanged, no proof failure.
Shuffle verified Peach minimum measured headroom36.60GiB and current87GiB;
root authorized one isolated -j1 -M4000 diagnostic with20GiB headroom guard.
It uses a separate tree, verifies reused import hashes, and must reproduce
Framework exactly before extending to the remaining11 foundations. Main
numeric concurrency unchanged. Goal ACTIVE; next remains occupation aggregate,
unconditional theorem/corollaries, invariant check, and completion audit.

Further continuation PROGRESS:114/133 occupation witnesses and344/433 indexed
boxes PASS. Root main waiter15735 and occupation controller processes verified
live. Package-aware diagnostic exceeded4GB compiler cap with ample host RAM;
single8GB retry succeeded under20GiB headroom guard. Ten foundation artifacts
matched exactly after preserving package metadata; remaining two reproduced
exactly using unchanged source bytes with original filename metadata via
Lean --stdin. All12 source checks and their axiom audits PASS; root independently
rehashed current sources, live objects, fresh objects, and recorded identities.
No binary editing or live artifact replacement occurred. Reports:
encoder_foundations_package_peach_verification.json and
encoder_foundations_filename_peach_verification.json. An optional11-module
semantic batch is assigned while matrix work remains; it must stop new
submissions when actual base assembly starts and is not an added completion
gate. Goal remains ACTIVE, with the same actual-theorem completion criteria.
Continuation PROGRESS:120/133 occupation matrix witnesses PASS; local/indexed
queues have consumed currently available witnesses, throughB389 in their
witness-sorted order. Both corrected wrappers B152/B372 are individually
kernel-checked, eliminating the known generator defects. All6 occupation
controllers/watchers and root final waiter are live. Eleven more old semantic
modules reproduced exactly with11 standard-only audits; total36 repaired
semantic modules now have direct matching current source/object evidence.
One last bounded optional arithmetic-helper/definition batch is assigned
under the same8GB/20GiB limits, stopping new submissions when base assembly
starts. No heavy historical numerical replay is added. Remaining matrix work
estimated12–14min; final theorem/corollary/invariant gates are still pending.
Goal remains ACTIVE. Next: finish those gates and inspect actual statements.

2026-09-28 closing assembly: all 133 occupation witnesses and all 433 local
and indexed boxes now PASS. Occupation aggregate and mixed DenseRates PASS;
all 1,023 global rectangles are certified. The actual unconditional native
theorem and full-statement pin compiled on Peach, with four standard-only
recursive audits and a 5,215-module source/object dependency snapshot.
All four success/existence corollary modules compiled, with seven
standard-only audits. Independent final-statement review found no hidden
premise or construction mismatch. Comment-aware scan of 5,299 project Lean
sources found zero prohibited constructs. Final default invariant check is
running; consolidated closure follows its PASS.

The timestamp provenance gate accepts 68 known timestamp inversions only
with direct compilation evidence: 63 existing paired low-cancellation
records and five freshly reproduced exact objects. No live object or mtime
was changed for this gate. Fifty modules have now been reproduced exactly
in isolated outputs. Final evidence coverage distinguishes 4,320 individual
compiler-pair records, 795 historical source-only records, 34 initial
snapshot-only records, and 66 final/transfer snapshot-only records. No
pending files or pair mismatches. A clean full-project source replay was
not run. Paper and original pin remain preserved against the checkpoint;
no commits. Next: final invariants, consolidated PASS, and goal completion.

2026-09-28 FINAL CLOSURE PASS: default build/invariants completed with ALL
CHECKS PASS and 5,219 actual theorem/corollary pairs preserved. Consolidated
audit passes all 17 gates and 5,220 recorded source/object pairs, with zero
issues. Root independently revalidated the 5,299-source trust-scan inventory
and scanner hashes. Actual native minimum-distance probability tends to
zero at floor(0.11N), every realization has rate1/2, and unconditional
success/existence corollaries pass. Current status, closure audit, parallel
assignment record, and result guide updated. Goal complete within agreed
native distance/rate scope. Recommended next: separate clean source replay
and packaging for review; no further mathematical lemma is pending.

2026-09-28 documentation and public-interface cleanup: added SpinCodes.Native
and NativePin without editing existing proof sources or generated certificates.
Both modules compiled on Peach; the concise pin exposes four complete types
and four standard-only recursive axiom audits. native_public_entry_verification.json
records the new sources/objects and links to the preserved closure evidence.
A separate scan covers all 5,301 current Lean sources with zero prohibited
constructs. The public import replay dry run discovers 5,218 modules; no full
replay was executed. Revalidated original closure: all 17 gates PASS.

Replaced the stale skeleton README; added PROOF_GUIDE, CONTRIBUTING, and the
scripts index; clarified reproduction scope/commands and marked ten earlier
documents as historical. Validated 76 local links across 17 documents. Disabled
the obsolete existence/old-log occupation controller before any imports or
I/O; its original implementation remains for historical inspection, and a
sample invocation exits1 with the retirement notice. Existing theorem sources,
paper, original pin, and certificate objects remain unchanged. No commits.
Next: present the native result in the paper with precise formalization scope;
retain independent full source replay as a separate verification follow-up.

2026-09-28 approved paper revision applied: five TeX sources updated with
formal distance/rate claim, appendix correspondence and conservative bounds,
conditioning clarification, artifact wording, and AI disclosure. The theorem
statement is unchanged; all 5,301 Lean source hashes still match the cleanup
inventory. Public84-page and anonymous82-page PDFs compile with resolved
references and no overfull boxes. Edited pages passed root and independent
visual QA; anonymous metadata/repository-link checks pass. Original proof
closure reports retain pre-revision paper hashes. The separate manuscript
record is scripts/map_data/paper_revision_verification.json. No commits.
Next: final reading of the revised discussion and versioned artifact packaging;
full project-source replay remains a separate follow-up.
