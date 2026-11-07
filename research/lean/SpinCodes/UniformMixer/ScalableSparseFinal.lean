import SpinCodes.UniformMixer.ScalableSparseSelected
import SpinCodes.Structured.EventualRegimes

/-! Summation of the full growing sparse occupation interval. -/

noncomputable section
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily
  Finset Filter

theorem scalableGoodLayer_sparse_sum_tendsto :
    Tendsto (fun m : ℕ =>
      ∑ Q ∈ Ico 4096 (nativeCut m + 1), scalableGoodLayer m Q)
      atTop (nhds 0) := by
  let ρ : ℕ → ℝ :=
    fun m => Real.exp (-(0.006 : ℝ) * (bsched m : ℝ))
  have hρ0 (m : ℕ) : 0 ≤ ρ m := (Real.exp_pos _).le
  have hρ1 (m : ℕ) : ρ m < 1 := by
    apply Real.exp_lt_one_iff.mpr
    have hb : (0 : ℝ) < bsched m := by
      exact_mod_cast bsched_pos m
    nlinarith
  have hρ : Tendsto ρ atTop (nhds 0) := by
    have hb := (tendsto_natCast_atTop_atTop (R := ℝ)).comp
      bsched_tendsto
    have hm := Tendsto.const_mul_atTop
      (by norm_num : (0 : ℝ) < 0.006) hb
    have hh := Real.tendsto_exp_atBot.comp
      (tendsto_neg_atTop_atBot.comp hm)
    apply hh.congr
    intro m
    change Real.exp (-(0.006 * (bsched m : ℝ))) =
      Real.exp (-(0.006 : ℝ) * (bsched m : ℝ))
    rw [neg_mul]
  apply sparse_regime_eventually scalableGoodLayer ρ nativeCut
    scalableGoodLayer_nonneg hρ0 hρ1 hρ
  filter_upwards [scalableGoodLayer_sparse_eventually] with m hm
  intro Q hQ
  refine (hm Q hQ).trans_eq ?_
  rw [show ρ m = Real.exp (-(0.006 : ℝ) *
    (bsched m : ℝ)) from rfl, ← Real.exp_nat_mul]
  congr 1
  ring

#print axioms scalableGoodLayer_sparse_sum_tendsto
end Spin.UniformMixer
