/-
The numeric constants of the BA sparse tail.

    c_o = 13/28,   ζ = √c_o,   κ = 4ζ/(1-ζ),   η = 1/1000,
    C_* = (4095 e/24)·(4κ/(1-8η))^4.

The paper asserts `C_* η³ < 0.657`.  That margin is thin — the true value is
`0.6565638…`, so the claim holds by about `4.4·10⁻⁴`, a relative margin of
`0.07%` — which makes it exactly the kind of constant worth checking rather
than trusting.  It does check out.

Everything here is rational arithmetic plus one bound on `e`; no logarithms
are involved, so `LogBounds` is not needed.  (The companion monotonicity claim
*does* need a logarithm and is tracked separately.)
-/
import Mathlib
import SpinCodes.Numeric.LogBounds

set_option linter.unusedSectionVars false

namespace Spin.Numeric

open Real

/-- `ζ = √(13/28)`. -/
noncomputable def zeta : ℝ := Real.sqrt (13 / 28)

lemma zeta_sq : zeta ^ 2 = 13 / 28 := Real.sq_sqrt (by norm_num)

lemma zeta_pos : 0 < zeta := Real.sqrt_pos.mpr (by norm_num)

lemma zeta_lt_one : zeta < 1 := by
  nlinarith [zeta_sq, zeta_pos]

/-- `0.6813851 < ζ < 0.6813852`. -/
lemma zeta_lower : (6813851 : ℝ) / 10000000 < zeta := by
  nlinarith [zeta_sq, zeta_pos]

lemma zeta_upper : zeta < (6813852 : ℝ) / 10000000 := by
  nlinarith [zeta_sq, zeta_pos]

/-- `κ = 4ζ/(1-ζ)`. -/
noncomputable def kappa : ℝ := 4 * zeta / (1 - zeta)

lemma one_sub_zeta_pos : 0 < 1 - zeta := by
  have := zeta_lt_one
  linarith

lemma kappa_pos : 0 < kappa := by
  unfold kappa
  have := zeta_pos
  have := one_sub_zeta_pos
  positivity

/-- `κ < 8.5544`. -/
lemma kappa_upper : kappa < (85544 : ℝ) / 10000 := by
  unfold kappa
  rw [div_lt_iff₀ one_sub_zeta_pos]
  nlinarith [zeta_upper, zeta_lower, one_sub_zeta_pos]

/-- `8.5543 < κ`. -/
lemma kappa_lower : (85543 : ℝ) / 10000 < kappa := by
  unfold kappa
  rw [lt_div_iff₀ one_sub_zeta_pos]
  nlinarith [zeta_upper, zeta_lower, one_sub_zeta_pos]

/-- `η = 1/1000`. -/
noncomputable def eta : ℝ := 1 / 1000

/-- `C_* = (4095 e/24)·(4κ/(1-8η))^4`. -/
noncomputable def Cstar : ℝ := (4095 * Real.exp 1 / 24) * (4 * kappa / (1 - 8 * eta)) ^ 4

/-- **`C_* η³ < 0.657`.**

The true value is `0.6565638…`.  The rational chain used here (via
`κ < 8.5544` and `4κ/(1-8η) < 34.4936`) proves the slightly weaker
`< 0.6565855`, which is still below `0.657` with margin `4.1·10⁻⁴`. -/
theorem Cstar_eta_cubed_lt : Cstar * eta ^ 3 < 657 / 1000 := by
  have he : Real.exp 1 < 2.7182818286 := Real.exp_one_lt_d9
  have hk : kappa < (85544 : ℝ) / 10000 := kappa_upper
  have hkpos : 0 < kappa := kappa_pos
  have hratio : 4 * kappa / (1 - 8 * eta) < (344936 : ℝ) / 10000 := by
    unfold eta
    rw [div_lt_iff₀ (by norm_num)]
    nlinarith [hk]
  have hrpos : 0 < 4 * kappa / (1 - 8 * eta) := by
    unfold eta
    positivity
  have hpow : (4 * kappa / (1 - 8 * eta)) ^ 4 < ((344936 : ℝ) / 10000) ^ 4 := by
    exact pow_lt_pow_left₀ hratio hrpos.le (by norm_num)
  have hfac : 4095 * Real.exp 1 / 24 < (4095 : ℝ) * 2.7182818286 / 24 := by
    have : (0 : ℝ) < 4095 / 24 := by norm_num
    nlinarith [he]
  have hfacpos : (0 : ℝ) < 4095 * Real.exp 1 / 24 := by
    have := Real.exp_pos (1 : ℝ)
    positivity
  unfold Cstar eta
  calc (4095 * Real.exp 1 / 24) * (4 * kappa / (1 - 8 * (1 / 1000))) ^ 4 * (1 / 1000) ^ 3
      < ((4095 : ℝ) * 2.7182818286 / 24) * ((344936 : ℝ) / 10000) ^ 4 * (1 / 1000) ^ 3 := by
        have h1 : (4 * kappa / (1 - 8 * ((1 : ℝ) / 1000))) ^ 4 < ((344936 : ℝ) / 10000) ^ 4 := by
          simpa [eta] using hpow
        have h2 : (0 : ℝ) < ((344936 : ℝ) / 10000) ^ 4 := by norm_num
        nlinarith [hfac, hfacpos, h1, h2]
    _ < 657 / 1000 := by norm_num


/-! ## The monotonicity check

`eq:structured-ba-sparse-moment` needs `(κℓ/(b-2ℓ+1))^ℓ` to decrease in `ℓ`
for `ℓ/b ≤ 12η`.  The paper checks this through the logarithmic derivative

    ln(κx/(1-2x)) + 1 + 2x/(1-2x) < -1.22   at x = 12η = 3/250.

The true value is `-1.2275267…`, so the margin is `0.0075`.

The argument of the logarithm is `≈ 0.1052`, far from `1`, where the log series
converges slowly.  Reducing by a factor of `8` moves it to `≈ 0.8424`, and then
four terms suffice for this margin. -/

/-- The argument of the logarithm at `x = 12η`. -/
noncomputable def logArg : ℝ := kappa * (3 / 250) / (1 - 2 * (3 / 250))

lemma logArg_pos : 0 < logArg := by
  unfold logArg
  have := kappa_pos
  positivity

/-- `κ < 8.5544` gives `logArg ≤ 0.1053`. -/
lemma logArg_le : logArg ≤ (1053 : ℝ) / 10000 := by
  unfold logArg
  rw [div_le_iff₀ (by norm_num)]
  nlinarith [kappa_upper]

/-- Four terms of the reduced series. -/
lemma series_lower : (17147 : ℝ) / 100000 ≤ Spin.Numeric.logSeries 4 ((1576 : ℝ) / 10000) := by
  unfold Spin.Numeric.logSeries
  simp only [Finset.sum_range_succ, Finset.sum_range_zero]
  norm_num

lemma rem_upper : Spin.Numeric.logRem 4 ((1576 : ℝ) / 10000) ≤ (12 : ℝ) / 100000 := by
  unfold Spin.Numeric.logRem
  rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1576 / 10000)]
  norm_num

/-- `ln(0.1053) ≤ -2.2507`, by reducing the argument with `ln 8 = 3 ln 2`. -/
lemma log_1053_le : Real.log ((1053 : ℝ) / 10000) ≤ -((22507 : ℝ) / 10000) := by
  have hsplit : Real.log ((8424 : ℝ) / 10000)
      = Real.log 8 + Real.log ((1053 : ℝ) / 10000) := by
    rw [show ((8424 : ℝ) / 10000) = 8 * (1053 / 10000) by norm_num]
    exact Real.log_mul (by norm_num) (by norm_num)
  have hlog8 : Real.log (8 : ℝ) = 3 * Real.log 2 := by
    rw [show (8 : ℝ) = 2 ^ (3 : ℕ) by norm_num, Real.log_pow]
    push_cast
    ring
  have hser := (Spin.Numeric.log_mem_Icc (q := (8424 : ℝ) / 10000)
    (by norm_num) (by norm_num) 4).2
  rw [show (1 : ℝ) - 8424 / 10000 = 1576 / 10000 by norm_num] at hser
  have hl2 : (0.6931471803 : ℝ) < Real.log 2 := Real.log_two_gt_d9
  have h1 := series_lower
  have h2 := rem_upper
  rw [hsplit, hlog8] at hser
  linarith

/-- **The monotonicity check.**

`ln(κx/(1-2x)) + 1 + 2x/(1-2x) < -1.22` at `x = 12η`. -/
theorem logDeriv_lt :
    Real.log logArg + 1 + 2 * (3 / 250) / (1 - 2 * ((3 : ℝ) / 250)) < -((122 : ℝ) / 100) := by
  have hle : Real.log logArg ≤ Real.log ((1053 : ℝ) / 10000) :=
    Real.log_le_log logArg_pos logArg_le
  have hb := log_1053_le
  norm_num
  linarith

end Spin.Numeric
