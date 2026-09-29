import SpinCodes.Structured.ConcreteOuterTailSparse

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric Finset Filter

lemma sparseLowerTail_nonneg (k : ℕ) : 0 ≤ sparseLowerTail k := by
  apply sum_nonneg
  intro x hx
  exact sum_nonneg (fun w _ => FinPMF.prob_nonneg _ _)

lemma sparseUpperTail_nonneg (k : ℕ) : 0 ≤ sparseUpperTail k := by
  apply sum_nonneg
  intro x hx
  exact sum_nonneg (fun w _ => FinPMF.prob_nonneg _ _)

theorem sparseLowerTail_tendsto : Tendsto sparseLowerTail atTop (nhds 0) := by
  have hk : Tendsto (fun k : ℕ => k*24) atTop atTop := tendsto_atTop_mono (fun k => by change k≤k*24; omega) tendsto_id
  have hh := (sparseBoundSeries_tendsto.comp hk).const_mul zeta⁻¹
  simp only [mul_zero] at hh
  apply squeeze_zero' (Eventually.of_forall sparseLowerTail_nonneg) _ hh
  filter_upwards [eventually_ge_atTop 1] with k hk
  exact sparseLowerTail_le hk

theorem sparseUpperTail_tendsto : Tendsto sparseUpperTail atTop (nhds 0) := by
  have hk : Tendsto (fun k : ℕ => k*24) atTop atTop := tendsto_atTop_mono (fun k => by change k≤k*24; omega) tendsto_id
  have hh := (sparseBoundSeries_scaled_tendsto.comp hk).const_mul (zeta⁻¹^2+1)
  simp only [mul_zero, Function.comp_def, ←mul_assoc] at hh
  apply squeeze_zero' (Eventually.of_forall sparseUpperTail_nonneg) _ hh
  filter_upwards [eventually_ge_atTop 1] with k hk
  exact sparseUpperTail_le hk

/-- The actual sparse-message expectation outside the window vanishes at both
ends; the upper-tail prefactor b has already been paid. -/
theorem sparse_outer_tails_tendsto :
    Tendsto (fun k => sparseLowerTail k+sparseUpperTail k) atTop (nhds 0) := by
  simpa using sparseLowerTail_tendsto.add sparseUpperTail_tendsto

end Spin.Structured.ConcreteOuter

