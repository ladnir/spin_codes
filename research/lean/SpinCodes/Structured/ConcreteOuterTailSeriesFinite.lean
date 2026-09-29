import SpinCodes.Structured.ConcreteOuterTailSeries

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset Filter

lemma sparseScaledSeries_eq (b : ℕ) :
    (∑' n, sparseScaledTerm b n) = (b:ℝ)*sparseBoundSeries b := by
  have hout : ∀ n ∉ range (b/1000), sparseScaledTerm b n=0 := by
    intro n hn
    have hh : ¬1000*(n+1)≤b := by
      have : ¬n<b/1000 := by simpa only [mem_range] using hn
      omega
    simp only [sparseScaledTerm, if_neg hh]
  rw [tsum_eq_sum hout, sparseBoundSeries, sum_Ico_eq_sum_range]
  simp only [Nat.add_sub_cancel]
  rw [mul_sum]
  apply sum_congr rfl
  intro n hn
  have hh : 1000*(n+1)≤b := by have := mem_range.mp hn; omega
  simp only [sparseScaledTerm, if_pos hh, Nat.add_comm 1 n]

theorem sparseBoundSeries_scaled_tendsto :
    Tendsto (fun b : ℕ => (b:ℝ)*sparseBoundSeries b) atTop (nhds 0) := by
  simpa only [sparseScaledSeries_eq] using sparseScaledSeries_tendsto

lemma sparseBoundSeries_nonneg (b : ℕ) : 0 ≤ sparseBoundSeries b := by
  apply sum_nonneg
  intro j hj
  exact pow_nonneg (mul_nonneg Cstar_nonneg (by positivity)) _

theorem sparseBoundSeries_tendsto : Tendsto sparseBoundSeries atTop (nhds 0) := by
  have hh := sparseBoundSeries_scaled_tendsto.div_atTop
    (tendsto_natCast_atTop_atTop : Tendsto (fun b : ℕ => (b:ℝ)) atTop atTop)
  apply hh.congr'
  filter_upwards [eventually_ge_atTop 1] with b hb
  have hb0 : (b:ℝ)≠0 := by exact_mod_cast (show b≠0 by omega)
  exact mul_div_cancel_left₀ _ hb0

end Spin.Structured.ConcreteOuter
