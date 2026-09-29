import SpinCodes.Structured.SparseConditioningMass
import SpinCodes.Structured.ConcreteMarkedReferenceMoment

noncomputable section
namespace Spin.Structured.ConcreteMarked

/-- The paper's conditioning and IMT factors for the actual fair-row reference. -/
theorem fairExperiment_sparse_probability_exp {L b R : ℕ} (S : Finset (Fin L))
    (hS : S.Nonempty) (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000)
    (hα : (S.card:ℝ)/L = α) (d : ℕ) :
    (fairExperimentLaw S b R).prob
      (fun ω => ConcreteRoutedEncoder.weight e ω.1 ω.2 ≤ d) ≤
        ((8*Real.sqrt (S.card:ℝ))^b *
          Real.exp (((L:ℝ)*b)*Spin.binKL α ((8/5)*α))) *
          (2048*(1-96*α)^R) / ((1-(8/5)*α)^d) := by
  have hp0 : 0 < (8/5:ℝ)*α := by positivity
  have hp1 : (8/5:ℝ)*α < 1 := by linarith
  have hc := markMass_inverse_pow_le (Finset.card_pos.mpr hS) (Routing.card_le_width S)
    hp0 hp1 b
  rw [hα] at hc
  refine (fairExperiment_sparse_probability S e hα0 hα1 d).trans ?_
  apply div_le_div_of_nonneg_right _ (pow_nonneg (by linarith) _)
  apply mul_le_mul_of_nonneg_right hc
  apply mul_nonneg (by norm_num)
  exact pow_nonneg (by linarith) R

end Spin.Structured.ConcreteMarked
