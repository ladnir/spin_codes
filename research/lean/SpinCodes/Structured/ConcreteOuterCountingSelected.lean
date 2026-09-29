import SpinCodes.Structured.ConcreteOuterCountingProbability

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset ConcreteMarked

/-- The sparse count in the framework's nonzero-message convention. -/
theorem nonzero_occupation_sparse_failure {L k R : ℕ} (seed : Seed k) (Q : ℕ)
    (e : (Fin (k * 24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) {B α : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (∑ x ∈ (Spin.nonzeroMsgs (Fin L → LocalMessage k)).filter (fun x => occupation x = Q),
      failureProbability e (rowSupports seed x) d) ≤
      (L.choose Q : ℝ) * (((2 : ℝ)^(k * 24) * B)^Q *
        ((((1 / markMass L Q ((8/5)*α)) ^ (k * 24)) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d))) := by
  apply le_trans _ (occupation_sparse_failure seed Q e d hB hα0 hα1)
  apply Finset.sum_le_sum_of_subset_of_nonneg
  · intro x hx
    exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, (Finset.mem_filter.mp hx).2⟩
  · intro x _ _
    exact FinPMF.prob_nonneg _ _

/-- The selected shared outer incurs its b² spectrum penalty once per active message row. -/
theorem good_nonzero_occupation_sparse_failure {L k R : ℕ} (seed : Seed k) (Q : ℕ)
    (e : (Fin (k * 24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) (W : Finset ℕ) {B α : ℝ}
    (hB0 : 0 ≤ B) (hgood : Spin.Good (seedLaw k) spectrum (k * 24) W seed)
    (hB : ∀ w ∈ W, Spin.Abar (seedLaw k) spectrum w ≤ B * ((k * 24).choose w : ℝ))
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (∑ x ∈ (Spin.nonzeroMsgs (Fin L → LocalMessage k)).filter (fun x => occupation x = Q),
      failureProbability e (rowSupports seed x) d) ≤
      (L.choose Q : ℝ) * (((2 : ℝ)^(k * 24) * (((k * 24 : ℕ) : ℝ)^2 * B))^Q *
        ((((1 / markMass L Q ((8/5)*α)) ^ (k * 24)) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d))) :=
  nonzero_occupation_sparse_failure seed Q e d (good_spectrum_bound seed W hB0 hgood hB) hα0 hα1

end Spin.Structured.ConcreteOuter
