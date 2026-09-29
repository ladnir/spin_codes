# Formalizing the Structured SPIN distance theorem in Lean 4

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

Target: `thm:structured-spin-scalable` (Scalable Structured SPIN), stated in
`research/paper/structured_proof.tex`, with certificates in
`research/paper/structured_appendix.tex` and
`research/paper/structured_imt_appendix.tex`.

No `sorry`, no project `axiom`, no compiled evaluation. Every result depends
only on Mathlib's `propext`, `Classical.choice`, `Quot.sound`; the numerical
certificates are checked by the Lean kernel (`sparseRows_cert` and
`Golay.hist_eq` do not even use choice). `scripts/check.sh` enforces all of
this, and `SpinCodes/Pin.lean` makes statement drift a build error.

## What was proved at this checkpoint

### Framework and architecture

| File | Content | Paper |
|---|---|---|
| `SpinCodes/Prob.lean` | finite probability: expectation, Markov, union bound, conditioning, product laws, `condition_prod_left`, uniform law | the measure-comparison bookkeeping |
| `SpinCodes/Framework.lean` | `expect_Z_eq`, `expect_Z_le_class_sum`, `prob_bad_le_class_sum` | `thm:class-first-moment` |
| `SpinCodes/Framework.lean` | `ZQ`, `Z_eq_sum_ZQ`, `prob_bad_le_cond_sum` | occupation refinement + the conditioning remark |
| `SpinCodes/Distance.lean` | `exists_nonzero_light`, `one_le_Z_iff` | "the bad event occurs exactly when `Z_d ≥ 1`" |
| `SpinCodes/Wrapper.lean` | `zeroExtend_preserves` | `eq:structured-requested-length-wrapper` |
| `SpinCodes/Selection.lean` | `prob_not_good_le` (quantitative) | `lem:structured-one-sample-selection` |
| `SpinCodes/Structured/Regimes.lean` | `fixed_regime`, `sparse_regime`, `dense_regime`, `total_regime_sum` | "Completing the Argument" |
| `SpinCodes/Structured/Certificate.lean` | `ScalableCertificate`, `distance_whp` | `thm:structured-spin-scalable` |
| `SpinCodes/Structured/Instantiation.lean` | `Family`, `ofFamily`, `ofNativeFamily` | the concrete construction |

### The schedule

`SpinCodes/Structured/Schedule.lean` — `eq:structured-native-schedule`.
`exists_sched` proves a positive multiple of 24 satisfying
`b ≥ (39/4) log₂(L_m b)` exists, so the family is non-empty and the theorem is
not vacuous. `bsched_mono` + `bsched_succ_le` give the `o(1)` relative gap
without needing `b_m = (39/4) log₂ N_m + O(1)`, and `wrapper_rate` gives
`N_{m(n)}/n → 1`.

### The outer majorant

`SpinCodes/Structured/Majorant.lean` — `â_BA` as the lower envelope of
rational affine supports, which is how the certificate represents it.
Concavity and `eq:structured-ba-jensen` follow from the representation, so they
hold whatever the 39 supports turn out to be.

### The exact BA spectrum — complete

`lem:structured-exact-ba-spectrum` is fully proved, Golay-instantiated, with
every ingredient computed or derived rather than assumed.

| File | Content |
|---|---|
| `Accumulator.lean` | `accCount_eq_accT`: `T_b(a,c)` counts exactly the words claimed. Splitting on the *last* letter reduces the closed form to one Pascal step per parity |
| `AccTuple.lean` | the same count over `Fin b → Bool`, plus `card_filter_ofFn_eq_countP` |
| `Interleaver.lean` | `card_perm_accWt_eq`: `P_b(a,c)` is the transition law of a *uniformly interleaved* accumulator — the paper's asserted "the second permutation makes the word uniform on its Hamming slice" |
| `Composition.lean` | `prob_two_stage`: two stages compose by Chapman–Kolmogorov. `sum_Pt` checks the law is a distribution |
| `Enumerator.lean` | `card_tuples_weight`: `#{k-tuples of weight a} = [u^a] W(u)^k` |
| `Golay.lean` | `hist_eq`: the weight distribution of all 4096 codewords, **computed by the kernel** |
| `GolayEnum.lean` | `enumerator_eq`: `G(u) = 1 + 759u⁸ + 2576u¹² + 759u¹⁶ + u²⁴`; `min_distance` |
| `BASpectrum.lean`, `BAGolay.lean` | `Ā_b(w) = Σ_{a≥1} G_b(a) Σ_c P_b(a,c) P_b(c,w)` |

### Numerical certificates

`SpinCodes/Structured/PolyCert.lean` + `SparseData/Row0..6` + `SparseCert.lean`
— all seven degree-257 residuals of `eq:imt-sparse-collatz` are kernel-checked.
The coefficients were regenerated from `certify_imt_sparse.py` and matched
against the `coefficients_sha256` recorded in `SPARSE_EXACT.json`, so the Lean
data is the certified polynomial.

`SpinCodes/Numeric/LogBounds.lean` — rigorous rational enclosures of
`Real.log`, the basis for the box covers. See `DECISIONS.md` D1.

## Work remaining at this checkpoint

The historical queue is in [TASKS.md](TASKS.md); `LOOP_LOG.md` records its iterations.

| Task | Content | Note |
|---|---|---|
| T3a | Bridge `T_occ(β,z)v ≤ (1-96α)v` to the seven certified polynomials | needs `T_occ`, so belongs with T9 |
| T6 | `lem:structured-ba-tails`, sparse half | pure inequalities, no boxes |
| T7 | BA dense tail + concave majorant box covers | needs LogBounds |
| T8 | `lem:structured-route-domination` | **Chernoff for a Poisson-binomial sum is not in Mathlib** — the first research-grade step |
| T9 | `app:imt-finite-transfers`: Krawtchouk/Parseval over `𝔽₂¹⁹`, fiber bounds, transfer induction | |
| T10 | Positive-occupation box cover | needs LogBounds |
| T11 | Fixed occupation `Q < 4096`: simplex matrix integrals, continuum limit | hardest; highest risk |
| T12 | `O(N)` ordinary and transposed encoder work | may be out of scope for a distance formalization |
| T1a | Concrete randomness parametrization | deferred until the outer/inner maps exist |

## Build

```bash
lake exe cache get && lake build
bash scripts/check.sh
```

Pinned to Lean `v4.34.0` / Mathlib `v4.34.0`. First build downloads the
Mathlib cache (~5 GB). The Golay module runs a ~2.5-minute kernel computation
on first build and is cached thereafter.
