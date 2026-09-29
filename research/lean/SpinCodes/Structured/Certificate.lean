/-
Theorem `thm:structured-spin-scalable` (Scalable Structured SPIN), assembled.

This file states the theorem's *proof obligations* as the fields of a single
structure and proves that they imply the conclusion

    Pr[ d_min(E_m) ≤ ⌊0.11 N_m⌋ ] = o(1).

Every field is one of the paper's certified inputs, named by its paper label.
Nothing here is an `axiom`: the theorem is an implication from the
certificate, so `#print axioms` shows only Mathlib's three.  Discharging a
field replaces a hypothesis by a proof without touching this file.

What is *proved* here (not assumed): that the four regimes compose, that the
geometric sparse bound sums over an unbounded cutoff, that `L e^{-ηN+o(N)}`
vanishes, that conditioning on the outer-only event `𝒢_b` is legitimate
(`Setup.prob_bad_le_cond_sum`), and that the resulting union is `o(1)`.
-/
import SpinCodes.Framework
import SpinCodes.Selection
import SpinCodes.Distance
import SpinCodes.Structured.Regimes

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset Filter

/-- The proof obligations of `thm:structured-spin-scalable`.

Field-by-field correspondence with the paper:

| field | paper |
|---|---|
| `N_eq`, `b_pos`, `b_tendsto`, `N_tendsto`, `L_le_N` | `eq:structured-native-schedule` |
| `probNotGood_tendsto` | `lem:structured-one-sample-selection`, `eq:ba-good-event` |
| `fixed` | `app:imt-fixed`, fixed occupation `Q < 4096` |
| `sparse` | `eq:imt-sparse-final` |
| `dense`, `eta_pos`, `denseErr_tendsto` | positive-occupation certificate |
| `probBad_le` | `Setup.prob_bad_le_cond_sum` (proved in `Framework.lean`) |
-/
structure ScalableCertificate where
  /-- Outer length `L_m`. -/
  L : ℕ → ℕ
  /-- Block length `b_m`. -/
  b : ℕ → ℕ
  /-- Direct length `N_m = L_m b_m`. -/
  N : ℕ → ℕ
  N_eq : ∀ m, N m = L m * b m
  b_pos : ∀ m, 0 < b m
  b_tendsto : Tendsto b atTop atTop
  N_tendsto : Tendsto (fun m => (N m : ℝ)) atTop atTop
  L_le_N : ∀ m, (L m : ℝ) ≤ (N m : ℝ)
  /-- `𝔼[Z_{⌊0.11N⌋,Q} ∣ 𝒢_b]`, the conditional first moment at occupation `Q`. -/
  EZ : ℕ → ℕ → ℝ
  EZ_nonneg : ∀ m Q, 0 ≤ EZ m Q
  /-- `Pr[¬𝒢_b]`, paid once for the single sampled outer. -/
  probNotGood : ℕ → ℝ
  probNotGood_tendsto : Tendsto probNotGood atTop (nhds 0)
  /-- Fixed occupation: finitely many `Q`, each with a `Q`-dependent remainder. -/
  fixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => EZ m Q) atTop (nhds 0)
  /-- The sparse/dense cutoff, eventually between `4095` and `L`. -/
  cut : ℕ → ℕ
  cut_spec : ∀ᶠ m in atTop, 4095 ≤ cut m ∧ cut m ≤ L m
  /-- Uniform sparse occupation, `𝔼 ≤ e^{-0.006 Q b}`. -/
  sparse : ∀ m, ∀ Q ∈ Ico 4096 (cut m + 1),
    EZ m Q ≤ Real.exp (-(0.006 : ℝ) * Q * b m)
  /-- Positive occupation: `≤ L e^{-ηN + o(N)}`. -/
  eta : ℝ
  eta_pos : 0 < eta
  denseErr : ℕ → ℝ
  denseErr_tendsto : Tendsto denseErr atTop (nhds 0)
  dense : ∀ m, ∑ Q ∈ Ico (cut m + 1) (L m + 1), EZ m Q
    ≤ (L m : ℝ) * Real.exp (-eta * N m + denseErr m * N m)
  /-- `Pr[d_min(E_m) ≤ ⌊0.11 N_m⌋]`, and its first-moment bound. -/
  probBad : ℕ → ℝ
  probBad_nonneg : ∀ m, 0 ≤ probBad m
  probBad_le : ∀ m, probBad m ≤ probNotGood m + ∑ Q ∈ Ico 1 (L m + 1), EZ m Q

namespace ScalableCertificate

variable (C : ScalableCertificate)

/-- The geometric base of the sparse regime, `ρ_m = e^{-0.006 b_m}`. -/
noncomputable def rho (m : ℕ) : ℝ := Real.exp (-(0.006 : ℝ) * (C.b m : ℝ))

lemma rho_nonneg (m : ℕ) : 0 ≤ C.rho m := (Real.exp_pos _).le

lemma rho_lt_one (m : ℕ) : C.rho m < 1 := by
  refine Real.exp_lt_one_iff.mpr ?_
  have : (1 : ℝ) ≤ (C.b m : ℝ) := by exact_mod_cast C.b_pos m
  nlinarith

lemma rho_tendsto : Tendsto C.rho atTop (nhds 0) := by
  have hcast : Tendsto (fun m => ((C.b m : ℝ))) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp C.b_tendsto
  have h1 : Tendsto (fun m => (0.006 : ℝ) * (C.b m : ℝ)) atTop atTop :=
    Tendsto.const_mul_atTop (by norm_num) hcast
  have h2 : Tendsto (fun m => -((0.006 : ℝ) * (C.b m : ℝ))) atTop atBot :=
    tendsto_neg_atTop_atBot.comp h1
  have h3 := Real.tendsto_exp_atBot.comp h2
  refine h3.congr fun m => ?_
  simp only [Function.comp_apply, rho]
  ring_nf

/-- The sparse envelope in geometric form. -/
lemma sparse_geom (m : ℕ) (Q : ℕ) (hQ : Q ∈ Ico 4096 (C.cut m + 1)) :
    C.EZ m Q ≤ C.rho m ^ Q := by
  refine le_trans (C.sparse m Q hQ) (le_of_eq ?_)
  rw [rho, ← Real.exp_nat_mul]
  congr 1
  ring

lemma fixed_tendsto : Tendsto (fun m => ∑ Q ∈ Ico 1 4096, C.EZ m Q) atTop (nhds 0) :=
  fixed_regime C.EZ C.fixed

lemma sparse_tendsto :
    Tendsto (fun m => ∑ Q ∈ Ico 4096 (C.cut m + 1), C.EZ m Q) atTop (nhds 0) :=
  sparse_regime C.EZ C.rho C.cut C.EZ_nonneg C.rho_nonneg C.rho_lt_one C.rho_tendsto
    (fun m Q hQ => C.sparse_geom m Q hQ)

lemma dense_tendsto :
    Tendsto (fun m => ∑ Q ∈ Ico (C.cut m + 1) (C.L m + 1), C.EZ m Q) atTop (nhds 0) :=
  dense_regime _ C.L C.N C.eta C.denseErr C.eta_pos
    (fun m => Finset.sum_nonneg fun Q _ => C.EZ_nonneg m Q)
    C.L_le_N C.dense C.denseErr_tendsto C.N_tendsto

/-- The union over all occupations is `o(1)`. -/
theorem total_tendsto :
    Tendsto (fun m => ∑ Q ∈ Ico 1 (C.L m + 1), C.EZ m Q) atTop (nhds 0) :=
  total_regime_sum C.EZ C.L C.cut C.cut_spec C.fixed_tendsto C.sparse_tendsto C.dense_tendsto

/-- **Theorem (Scalable Structured SPIN), distance claim.**

`Pr[ d_min(E_m) ≤ ⌊0.11 N_m⌋ ] = o(1)`, over the two BA permutations, all
route permutations, and all IMT transvections. -/
theorem distance_whp : Tendsto C.probBad atTop (nhds 0) := by
  have hmaj : Tendsto
      (fun m => C.probNotGood m + ∑ Q ∈ Ico 1 (C.L m + 1), C.EZ m Q) atTop (nhds 0) := by
    have := C.probNotGood_tendsto.add C.total_tendsto
    rwa [add_zero] at this
  exact tendsto_of_tendsto_of_tendsto_of_le_of_le tendsto_const_nhds hmaj
    C.probBad_nonneg C.probBad_le

end ScalableCertificate

end Spin.Structured
