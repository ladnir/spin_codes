import SpinCodes.Structured.ConcreteMarkedFairRows
import SpinCodes.FiniteMoment

noncomputable section
namespace Spin.Structured.ConcreteMarked
open ConcreteRoute ConcreteEncoder

theorem fairRows_route_expect {L : ℕ} (S : Finset (Fin L)) (b : ℕ)
    (F : (Fin b → Finset (Fin L)) → ℝ) :
    (fairRows S b).expect (fun rows => (seedLaw L b).expect (fun seed => F (routeEval rows seed))) =
      (regionsLaw S b).expect F := by
  rw [← fairRows_route, FinPMF.expect_bind]
  simp_rw [← permutation_law_eq, FinPMF.expect_map]

/-- Fair active-row values, actual route permutations, and transvections are independent. -/
def fairExperimentLaw {L : ℕ} (S : Finset (Fin L)) (b R : ℕ) :=
  (fairRows S b).prod (ConcreteRoutedEncoder.experimentLaw L b R)

theorem fairExperiment_sparse_moment {L b R : ℕ} (S : Finset (Fin L))
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (fairExperimentLaw S b R).expect
      (fun ω => (1 - (8/5)*α) ^ ConcreteRoutedEncoder.weight e ω.1 ω.2) ≤
        ((1 / markMass L S.card ((8/5)*α)) ^ b) *
          (2048 * (1 - 96*α) ^ R) := by
  rw [fairExperimentLaw, FinPMF.expect_prod]
  simp_rw [ConcreteRoutedEncoder.moment_eq]
  rw [fairRows_route_expect S b
    (fun regions => inputMoment (1-(8/5)*α) (reshape e regions) ∅)]
  exact sparse_moment S e hα0 hα1

theorem fairExperiment_sparse_probability {L b R : ℕ} (S : Finset (Fin L))
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) (d : ℕ) :
    (fairExperimentLaw S b R).prob
      (fun ω => ConcreteRoutedEncoder.weight e ω.1 ω.2 ≤ d) ≤
        (((1 / markMass L S.card ((8/5)*α)) ^ b) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d) := by
  apply FinPMF.prob_weight_le_of_moment _ _ d (by linarith : 0 < 1-(8/5:ℝ)*α)
    (by linarith : 1-(8/5:ℝ)*α ≤ 1)
  exact fairExperiment_sparse_moment S e hα0 hα1

end Spin.Structured.ConcreteMarked

