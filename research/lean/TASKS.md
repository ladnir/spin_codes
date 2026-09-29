# Work queue — Structured SPIN formalization

> **Historical snapshot — archived 2026-09-28.** Status statements and task recommendations below describe an earlier stage.
> See [STATUS.md](STATUS.md) and [NATIVE_RESULT.md](NATIVE_RESULT.md) for the closed native distance/rate result and verification scope.

This historical queue used one task per loop iteration. Large tasks could
span iterations; later tasks waited until the current task was `DONE` or `BLOCKED`.

## Historical loop protocol

The protocol at this checkpoint was:

1. Read this file. Take the first task not `DONE`/`BLOCKED`.
2. Do the work. Prefer removing a `ScalableCertificate` field over adding
   library surface.
3. Run `bash scripts/check.sh`. It must pass. If it cannot be made to pass,
   revert the iteration's edits and record why.
4. Update this file (status + one-line note) and `LOOP_LOG.md`.

### Invariants — never violate to make progress

- `lake build` clean, `SpinCodes/Pin.lean` included.
- Zero `sorry`, zero `admit`, zero project `axiom`.
- Axiom closure exactly `[propext, Classical.choice, Quot.sound]`.
- **`SpinCodes/Pin.lean` is not edited to make the build pass.** If a task
  genuinely requires a statement change, stop the loop and report it. A
  changed pin is a claim about the paper and needs review.
- No hypothesis added to a lemma merely to close a goal. If a lemma needs a
  hypothesis the paper does not state, that is a finding: record it and stop.

### Stop conditions

- Queue empty.
- Pin change required.
- Same task fails 3 iterations running (thrash).
- A certificate field is found false, or unprovable as stated, under the
  paper's hypotheses. **This is the valuable outcome** — report it in full:
  the field, the obstruction, and whether it looks like a formalization
  artifact or a real gap in the paper.

## Queue at this checkpoint

| # | Task | Status | Note |
|---|---|---|---|
| T1 | Concrete per-`m` family. | **DONE** | `Structured/Instantiation.lean`. `probBad_le`, `probBad_nonneg`, `EZ_nonneg` now proved via `ofFamily`. Field kept (not deleted) so the pin is untouched; randomness spaces left abstract — see note below. |
| T1a | Instantiate the paper's actual maps and shared outer/routing/transvection randomness, occupation count, and distance threshold. | **IN PROGRESS — 2026-09-25** | The exact maps, ranks, XOR laws, feedback columns, transpose pairing, concrete Fourier/Parseval identities, kernel distance, and four fiber bounds now compile. The 512-block spectrum replay is running; numerical spectrum/kernel/cap identification, low-input cancellation data, one-step law, and the actual family remain. See `CONCRETE_MAPS.md` and the scoped verification reports. |
| T2 | Jensen + native schedule. | **DONE** | `Structured/Majorant.lean`, `Structured/Schedule.lean`. Five schedule fields discharged via `ofNativeFamily`. |
| T2a | Zero-extension wrapper. | **DONE** | `Wrapper.lean`. Injectivity + weight preservation proved. |
| T2b | Schedule gap + wrapper rate. | **DONE** | `Schedule.lean`: `Nsched_ratio_tendsto_one`, `mOf`, `wrapper_rate`. Proved via monotonicity + step ≤ 24, so no `O(1)` characterization of `b_m` was needed. |
| T3 | Sparse Collatz certificate — numerical core. | **DONE** | `PolyCertDefs/PolyCert/SparseData/SparseCert`. All 7 degree-257 rows kernel-checked; closure `[propext]`. Data regenerated from `certify_imt_sparse.py` and SHA-matched to `SPARSE_EXACT.json`. |
| T3a | Bridge `T_occ(β,z)v(α) ≤ (1-96α)v(α)` to the seven certified polynomials — i.e. prove the polynomials really majorize the residuals. Needs `T_occ`, so it belongs with T9. | **DONE — 2026-09-25** | `Sparse.sparse_collatz` proves the numerical-matrix inequality for `0 < α ≤ 1/10000`; all cancellation choices, six live rows, endpoint branches, and 129-weight mixture are connected. `sparse_iterate_bound` gives the prefactor-2048 bound for every `R`. Fresh semantic replay, pins, standard-axiom audit, and invariant check pass. Concrete map identification remains T1a. |
| T4 | Exact accumulator transition `T_b(a,c)`. | **DONE** | `Structured/Accumulator.lean`. `accCount_eq_accT`: the closed form counts exactly the words claimed. |
| T4b | Uniform interleaver + tuple bridge. | **DONE** | `AccTuple.lean`, `Interleaver.lean`. `card_perm_accWt_eq`: for fixed `u` of weight `a` and uniform `τ`, `Pr[accWt(u∘τ)=c] = T_b(a,c)/C(b,a)`. This is the paper's asserted "the second permutation makes the word uniform on its Hamming slice". |
| T4c | Two-stage composition. | **DONE** | `Composition.lean`. `prob_two_stage`: Chapman-Kolmogorov for the BA-3 inner pair. `sum_Pt` confirms the transition law is a probability distribution. |
| T4f | The exact BA spectrum. | **DONE** | `BASpectrum.lean`, `BAGolay.lean`. `lem:structured-exact-ba-spectrum` complete, Golay-instantiated, every ingredient proved. |
| T4d | Direct-sum enumerator. | **DONE** | `Enumerator.lean`. `card_tuples_weight`: `#{k-tuples of total weight a} = [u^a] W(u)^k`, stated over an arbitrary finite block type. |
| T4e | Golay weight enumerator. | **DONE** | `Golay.lean`, `GolayEnum.lean`. `G(u)` computed by kernel evaluation of all 4096 codewords, not assumed. Minimum distance 8 falls out. |
| T5 | Interval-arithmetic decision. | **DONE** | `DECISIONS.md` D1: build, do not import. `Numeric/LogBounds.lean` prototypes and validates the enclosure. |
| T6 | BA sparse tail, combinatorial core. | **DONE** | `BATails.lean`. `eq:structured-ba-one-acc-tail` up to the middle equality, division-free in `ℕ`. |
| T6b | `eq:structured-ba-one-acc-tail`, even case. | **DONE** | `BATails.lean`. `sum_accT_even_le_pow`: `p_{2ℓ} ≤ (13/28)^ℓ`, division-free in `ℕ`. |
| T6d | Odd-weight tail case. | **DONE** | `BATails.lean`. `sum_accT_odd_le_pow`. `eq:structured-ba-one-acc-tail` now holds for every input weight. |
| T6c | Negative-binomial bound. | **DONE** | `NegBinom.lean`. `moment_sum_le`: `Σ_{h≥ℓ} C(h-1,ℓ-1)z^h ≤ (z/(1-z))^ℓ` for partial sums — no convergence argument needed. |
| T6e | Sparse moment ratio. | **DONE** | `BATails.lean`. `choose_ratio_le`: the last step of `eq:structured-ba-sparse-moment`, `ζ`-free. |
| T6f | Split-and-bound structure of `S_b`. | **DONE** | `Regimes.lean`. `sum_split_geom` + `pow_mul_rpow_sqrt_tendsto`. |
| T6h | The constant `C_*η³ < 0.657`. | **DONE** | `Numeric/BAConstants.lean`. `ζ`, `κ` enclosed; `Cstar_eta_cubed_lt` proved. Margin is `4.1e-4` — verified at 50 digits first. |
| T6j | The monotonicity check. | **DONE** | `Numeric/BAConstants.lean`. `logDeriv_lt`. First use of the `LogBounds` enclosure; argument reduction by 8, four series terms. |
| T6k | `b·S_b → 0`. | **DONE** | `Regimes.lean`. `sparse_sum_mul_tendsto`, plus `tendsto_sqrt_atTop` extracted for reuse. |
| T6g | Upper tail: the complement identity. | **DONE** | `Accumulator.lean`. `accWtL_flipHead` / `accWtL_ge_iff`: `Acc(V+e₁) = ¬Acc(V)`, so the upper tail *is* a lower tail. |
| T6i | Upper-tail counting step. | **DONE** | `UpperTail.lean`. `card_upper_tail_le`: upper-tail words inject into the two neighbouring lower-tail classes. Division-free. |
| T6l | Upper tail vs. the neighbouring lower tails. | **DONE** | `UpperTail.lean`. `upper_tail_le_sums`, `choose_pred_ratio`, `choose_succ_ratio`, `sum_accT_eq_card`. |
| T6m | Upper tail at odd weight, numerically. | **DONE** | `UpperTail.lean`. `upper_tail_odd_bound`: `U·28^{m+1} ≤ 41·13^m·b·C(b,2m+1)`. Edge-case free. |
| T6n | Upper tail at even weight. | **DONE** | `UpperTail.lean`. `upper_tail_even_bound`. With T6m, **the sparse half of `lem:structured-ba-tails` is complete for every input weight.** |
| T7a | Integer fixed-point arithmetic with directed rounding (`ℤ` numerators over a fixed denominator) plus soundness against `ℝ`. | **DONE** | `Numeric/FixedDefs.lean` (defs, import-free) + `Numeric/Fixed.lean` (soundness). `Fix.Mem` pinned by `Iff.rfl` so it cannot be weakened to a vacuous predicate. Measured: **~21 ms per box** at 8 operations. |
| T7a2 | Fixed-point `log`: the `RatLog` series and power-of-two reduction rebuilt on `Fix`, with soundness chained to `Real.log`. | **DONE** | `Numeric/FixedLog.lean`. `flog_mem` pinned with its window hypotheses. Horner evaluation, window `[3/4, 3/2]`. Measured **~280 ms per log** at `n = 48`; `n` is the caller's knob and `n = 36` already gives `1.5e-11`. |
| T7b1 | Interval reciprocal, division, and `log` of an enclosure — the primitives `π` and `g` need beyond the ring operations. | **DONE** | `inv` / `div` / `flogI`, soundness and pins. `div = mul a (inv b)` reuses the corner analysis; `flogI` is two `flog` calls by monotonicity. |
| T7b2 | The real-side definitions `hEnt`, `piBA`, `gObj`, `gBA`, the identity `piBA = piEval`, and the infimum-witness bound. | **DONE** | `Numeric/BAExponent.lean`. `piBA_eq_piEval` proved; `gBA_le` needs `0 ≤ a ≤ 1` because the family is genuinely unbounded below above 1. `GolayG_eq_sum` derives `759/2576/759` from the verified code rather than retyping them. |
| T7b2b | The `Fix` evaluators for `hEnt`, `piEval` and `gObj`, with soundness by composition from T7b1. | **DONE** | `Numeric/BAEvalDefs.lean` + `Numeric/BAEval.lean`. Shifts are searched for and the window is *checked*, with the result in `Option` — so no soundness lemma needs a hypothesis about the search. |
| T7b2c | An `xlogx` interval extension, with `fhEnt` and `fpiEval` rebuilt on it. | **DONE** | Upper end by the maximum principle (`Real.convexOn_mul_log`), lower end by the tangent `x·(log m + 1) - m` at the box midpoint — no `1/e` constant needed. The degenerate `{0}` box (the boundary `c = a/2`) is its own branch. |
| T7b3 | **Scoped, not started.** The cover replay. Measurements in iteration 33 settled the design; see `LOOP_LOG.md`. Split below. | — | |
| T7b3a | The odd (`atanh`) series, replacing the `log(1-x)` series in the box path. | **DONE** | `oddSeries` / `flogA` / new `log2`, window `[2/3, 3/2]`. Measured **42 ms a log at m=9, against 245 ms** — 5.8x — and the enclosures came out tighter, not looser. |
| T7b3b | Tabulation: a lookup tree for the distinct log arguments, verified once, with the evaluators generalized over a log oracle. | **DONE** | `LogTree` / `Oracle` / `treeLog_oracle`. Soundness needs **no ordering invariant** — `find` answers only where it compared the key equal. Measured **2.3 ms a lookup against 42 ms a computation**. |
| T7b3c1 | The partial derivatives `Dpa`, `Dpc` of `piEval`, as `HasDerivAt`. | **DONE** | `BAExponent.lean`. Match the verifier's `p_gradient` to 1e-15. |
| T7b3c2 | The centered bound: two-step mean-value assembly, plus its `Fix` evaluator. | **DONE** | `piEval_centered` (real) and `fpiCenteredL` / `fpiCenteredL_mem` (`Fix`). On a box where the naive bound gives **+0.0057** the centered one gives **-0.0066**, against a true max of `-0.0071`. |
| T7b3d | The certificate as a **split tree**, with the covering property as structural induction. | **DONE** | `BoxTree` / `checkTree` / `checkTree_sound`, parametric in the leaf test and the property. Plus `meet`, `fpiCenteredG` (guard checked, not assumed) and `fpiBestL`, which takes whichever bounds are available. |
| T7b3e | The leaf predicate and the bridge from a checked tree to the paper's claim. | **DONE** | `infeasible` / `leafOK` / `DenseClaim` / `leafOK_sound` / `denseTail_of_checkTree`. Both leaf kinds tested: a feasible box encloses `[-0.1477, -0.1183]` around the true `[-0.1447, -0.1232]` and passes; an infeasible box passes vacuously. |
| T7b3f | The generator. | **DONE** | `scripts/fixmirror.py` (cross-checked against Lean, 22 cases, 0 mismatches) + `scripts/emit_cover.py`. |
| T7b3g | Emit the cover as Lean data and instantiate the bridge. | **DONE** | 44 self-contained parts + `DenseTail.lean` + hand-written `DenseTailReal.lean`. |
| T7b3h | Check the full cover and record the run. | **DONE — 2026-09-23** | **44/44 parts + root, ~63 min wall at 6-way batching, bounded memory.** `denseTail_real` closes over `[propext, Classical.choice, Quot.sound]`; every `logsN_ok` over **no axioms at all**. Run with `bash scripts/check-cover.sh full`. |
| T7c | Concave majorant cover (`eq:ba-spectrum-majorant`). | **DONE — 2026-09-25** | `Majorant/Refined.lean`: the full variational bound for the current paper's exact 39 supports, on the closed window and feasibility region. Whole-objective Taylor bound recovered from the original verifier; 54,032 leaves, 387/387 kernel-checked numerical modules, 19/19 segment assemblies, final theorem compiled. Central support proved analytically and right supports by reflection. Existing pins unchanged; new pins in `Majorant/Pin.lean`. See `DECISIONS.md` D5 and `CLOSURE_AUDIT.md`. |
| T8 | `lem:structured-route-domination`. Split below. | — | |
| T8a | The pigeonhole behind two of the three ingredients: a mass function puts at least `1/n` on its largest atom; plus the binomial mass function as a `FinPMF`. | **DONE** | `Structured/RouteDomination.lean`. `le_mode`, `inv_mode_le`, `binPMF`, `types_lower_bound` (mode hypothesis still open). |
| T8b | Discharge the mode hypothesis at `u = k/L`. | **DONE** | `binPMF_mode_nat` / `binPMF_mode` / `types_lower_bound_at`. Clearing denominators by `L^L` turns it into a natural-number inequality needing two factorial-versus-power bounds and one symmetry — no ratio test, no Stirling. |
| T8c1 | Poisson-binomial MGF domination `∏(1 + pᵢc) ≤ (1 + qc)^L`. | **DONE** | `prod_one_add_mul_le`. Uniform-weight AM-GM (`Real.geom_mean_le_arith_mean_weighted`); equality exactly when the `pᵢ` agree. |
| T8c2 | The Chernoff bound itself, and the pointwise domination it yields. | **DONE** | `poissonBinom` / `poissonBinom_mgf` / `poissonBinom_chernoff` / `poissonBinom_prob_le` / `poissonBinom_le_binom`. Stated without `exp` or `log`: at the optimal `z` the bound is `q^k(1-q)^(L-k)L^L/(k^k(L-k)^(L-k))`, whose factor cancels exactly against T8b, giving `P[S=k] ≤ (L+1)·Bin(L,q)[k]` with no entropy algebra at all. |
| T8d | Assemble `E_route[F] ≤ (b+1)^Q (L+1)^b E_Ber[F]`, and the `exp(o(N))` corollary. | **DONE** | `Structured/Domination.lean`. `Dominates` calculus (`expect_le`, `trans`), `dominates_of_fiber_uniform` (the shuffle step), `piPMF` + `dominates_piPMF` (costs multiply), `dominates_condition_le` (the row step), `poissonBinom_le_binom_all` (endpoints, which need AM-GM directly — prompted factoring `prod_le_pow_avg` out of T8c1), `poissonBinom_const_prob` (the dominating law really is the binomial), `dominates_region`, `route_domination`, `route_factor_isLittleO`. |
| T9 | `app:imt-finite-transfers`. Split below. | — | |
| T9a | Binary Krawtchouk polynomials and the weight-layer character sum `∑_{|x|=j} (-1)^{|x∩v|} = K_j(|v|)`. | **DONE** | `Structured/Krawtchouk.lean`. Not in Mathlib (neither Krawtchouk nor MacWilliams). Proved by splitting the layer at `|x∩v|`; the fibre bijection `x ↦ (x∩v, x
)` gives `C(w,h)·C(n-w,j-h)`. |
| T9b | Character orthogonality and the Fourier/Parseval identities. | **DONE** | `Structured/Parseval.lean`. `chi`, `sum_chi` (orthogonality on `(D,∆)`), `card_orth_layer` (inversion), `card_orth_pairs` (Parseval). Both proved with **no syndrome map**: `V_j` in coset form, since two words share a syndrome iff their symmetric difference lies in the kernel. So the whole development needs only a subgroup of `(Finset (Fin n), ∆)` — no matrix, no quotient. |
| T9c1 | Three of the four integer fibre bounds: `a_j`, the `∑|K_j|` Fourier bound, and constant-weight packing. | **DONE** | `Structured/FiberBounds.lean`. Needed generalising inversion to an arbitrary syndrome (`card_fiber`), which also made `card_orth_layer` its `x = ∅` case. Packing is stated for an abstract constant-weight set, so the complemented form (`h = n-j`) is a second corollary rather than a second proof. |
| T9c2 | The fourth bound: the variance/Cauchy-Schwarz bound. | **DONE** | `Structured/VarianceBound.lean`. Stated as the *quadratic* `(m+1)t² + a² ≤ 2at + m·∑|Φ|²`, which is the appendix formula before it is solved for `t` — so no square root and no rounding. Needed the fibre partition: `sum_card_fibers` (`∑|Φ| = a_j`) and `card_orth_pairs_split` (`V_j = k_j² + ∑|Φ|²`). |
| T9d1 | The IMT transvection law `L(Mq) = ½δ_q + ½Unif(nonzero)` for `q ≠ 0`. | **DONE** | `Structured/Transvection.lean`. Stated as three exact counts over the sample space, so no `2^(s-1)` and no edge case in `s`. Turns on one lemma: a parity functional splits a `∆`-closed family in half. |
| T9d2 | The transvection law as an expectation: `E[g(Mq)] = ½ g(q) + (1/(2M)) ∑_{w≠0} g(w)`. | **DONE** | `sum_act_eq` (integer form) / `expect_act` (normalised) / `card_transPairs`. An *exact* identity for every `g`, no sign condition — this is the form every row of `T_j(z)` consumes. |
| T9d3 | The weighted-state abstraction and the transfer row out of a nonzero state. | **DONE** | `Structured/Transfer.lean`. `ShellSystem` (shells abstract, any partition of the nonzero states), `Coords`, `Dominates`, `total_le`, `mono`, and `dominates_actLaw` — the row is *exact*, and its coordinate total is exactly `1`, so the abstraction is lossless there. Also `prob_act`, the law in pointwise form. |
| T9d4 | The step kernel: the transvection law composed with the syndrome shift; the zero row and zero column. | **DONE** | `Structured/Step.lean`. `step_sum` is *exact*; `step_zero_column` and `step_from_zero` are corollaries, and `sum_refresh_le_min` derives the paper’s `min(m_{d,j}, ν_j)`. Needs neither `A` nor `C` to be linear — `A` enters only through the weight, `C` only as the shift. Inputs are an arbitrary weighted finite set, so no normalisation. |
| T9e | The transfer induction `eq:imt-moment`, plus mixtures (`T_occ`) and the scalar envelope (`T_sc`). | **DONE** | `Structured/Induction.lean`. `Transfer`/`apply` (row-vector convention, pinned by `rfl`), `dominates_eZ`, `imt_moment`, `imt_moment_eZ`, `dominates_smul`/`dominates_sum` for the binomial mixture, `total_apply` and `total_iterate_le` for `T_sc`. The row-domination claim stays a hypothesis — that is where the verifier’s envelopes enter. |
| T10a | The KL chain-rule bound `D(αx‖py) ≤ D(α‖p) + α D(x‖y)` behind `eq:imt-dense-exponent`. | **DONE** | `Structured/KLChain.lean`. `logSum2` (two-term log-sum, zero numerators allowed) and `binKL_chain`. Holds on the paper’s *closed* range — `α, x` up to `1` — with the degenerate terms vanishing under `0 log 0 = 0`. |
| T10b | Convexity in `(α, αx)` and the vertex-maximum principle. | **DONE** | `Structured/Convexity.lean`. `le_max_vertices` (four corners suffice), `convexOn_perspective` (the perspective transform — not in Mathlib, and what makes `α·D(x‖y)` convex in `(α,αx)`), `convexOn_binKL_left`, and the slice lemmas. The box must lie in the cone `{0 < w/u < 1}`; outside it the exponent is not even real-valued. |
| T10c | The Collatz bound `T w ≤ λ w` ⇒ `e_Z T^R 1 ≤ λ^R·w_Z/w_min`. | **DONE** | `Structured/Collatz.lean`. `Coords.dot`, `Transfer.applyCol`, `dot_apply` (adjointness `(cT)·w = c·(Tw)`, the whole argument), `collatz_total`. Stated multiplicatively — no division — and `0 ≤ λ` is derived, not assumed. T10 is now complete. |
| T11a | Empty-epoch estimates: the mixing recurrence, `|E S_g - μg| ≤ 2t`, `Var(S_g) ≤ 3t²g`, and `G_Q ≤ G_3`. | **DONE** | `Structured/EmptyEpoch.lean`. `P = (I+Π)/2` collapses to the scalar recurrence `m_{g+1} = (m_g+μ)/2`, so no probability space and no matrix exponential is needed. The `2t` constant is sharp. |
| T11b | The lift/projection `RJ = I₂` and the coarse-grained empty kernel `J D_γ R`. | **DONE** | `Structured/Coarse.lean`. `proj_lift`, `lift_coarseD_proj` (the live block is *constant* at `d/M`, so `eq:imt-empty-limit` is exactly the entrywise distance to `J D R`), and `abs_expNeg_sub_le`, the Lipschitz step. |
| T11c1 | The `K_a` integrands, the scalar integrals, and `K_1^+`. | **DONE** | `Structured/Simplex.lean`. Integrands checked against Mathlib’s *own* matrix product (`integrand_one`, `integrand_two`), stated as the raw products it yields with the exponent arithmetic split off (`exp_prod_one`/`exp_prod_two`); three scalar integrals by exhibited antiderivatives, no substitution or parts; `K_1^+` entry by entry. |
| T11c2 | `K_2^+` by integration over the 2-simplex. | **DONE** | All four entries, matching the paper’s `[[rg, e], [re, re+a₀]]`. Written as a nested `intervalIntegral` with a variable inner endpoint (`u₂ = 1-u₀-u₁`), so no Fubini and no measure theory — each level is one antiderivative. T11c complete. |
| T11d1 | The weighted row norm: submultiplicativity, `‖D_σ‖_v ≤ 1`, `‖(I+uP_+)D_σ‖_v = m_v(u)`, and the product bound. | **DONE** | `Structured/WeightedNorm.lean`. Carried as a predicate `RowNormLe`, so no supremum and no division by `v₁`. `m_v(u)` is the norm *exactly* (checked). Nonnegativity of `u, r, σ, v₁` is all that is used — never stochasticity, as the paper notes. |
| T11d2 | The Beta identity `E[(1-Y+Y/σ)^{-(Q+1)}] = σ^l`. | **DONE** | `Structured/BetaSub.lean`. The Beta value **cancels**, so no Beta/Gamma machinery is needed — only the change of variables `y = uσ/(1-u(1-σ))`, under which every power of the denominator cancels and the integrand becomes *exactly* `σ^{a+1}u^a(1-u)^b`. T11 complete. |
| T12 | Rate and linear work: `O(N)` bit operations for the ordinary and transposed encoders. | **OUT OF SCOPE (work half) / DONE (rate half)** | The theorem’s *rate* claim is already covered: `Wrapper.lean` (injectivity, weight preservation) and `Schedule.lean`’s `wrapper_rate` give rate `1/2 - o(1)`; native rate exactly `1/2` is a dimension count of the construction, which this development deliberately does not model because nothing in the distance argument uses it. The *`O(N)` bit operations* claim is **declined**: it needs a cost model for bit operations, which Lean/Mathlib does not have. Any model defined here would be one I wrote, the mathematical step (`N/128` steps × constant) is trivial, and all the content sits in whether the model faithfully represents bit operations — which a Lean proof cannot establish. A model charging `1` per step makes the claim true by construction, i.e. exactly the vacuity failure mode the pin discipline exists to catch. `distance_whp` does not depend on this claim. |

## Not committed

The historical loop did not run `git commit`. Changes accumulated in the
working tree, and `LOOP_LOG.md` recorded each iteration.


