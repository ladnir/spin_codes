import SpinCodes.Structured.ConcreteOuterCountingSelected

namespace Spin.Structured.ConcreteOuter
open Finset ConcreteMarked

example {L k R : ℕ} (seed : Seed k) (Q : ℕ)
    (e : (Fin (k * 24) × Fin L) ≃ (Fin R × Fin 128)) (d : ℕ) {B α : ℝ}
    (hB : ∀ w ≤ k * 24, (spectrum seed w : ℝ) ≤ B * ((k * 24).choose w : ℝ))
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (∑ x ∈ (Spin.nonzeroMsgs (Fin L → LocalMessage k)).filter (fun x => occupation x = Q),
      (ConcreteRoutedEncoder.experimentLaw L (k * 24) R).prob
        (fun ω => ConcreteRoutedEncoder.weight e (rowSupports seed x) ω ≤ d)) ≤
      (L.choose Q : ℝ) * (((2 : ℝ)^(k * 24) * B)^Q *
        ((((1 / markMass L Q ((8/5)*α)) ^ (k * 24)) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d))) :=
  nonzero_occupation_sparse_failure seed Q e d hB hα0 hα1

#print axioms Spin.Structured.ConcreteOuter.failureProbability_eq_route
#print axioms Spin.Structured.ConcreteOuter.failureProbability_fair_eq
#print axioms Spin.Structured.ConcreteOuter.activeMessages_failure_le_fair
#print axioms Spin.Structured.ConcreteOuter.activeMessages_sparse_failure
#print axioms Spin.Structured.ConcreteOuter.sum_occupation_eq_activeSets
#print axioms Spin.Structured.ConcreteOuter.occupation_sparse_failure
#print axioms Spin.Structured.ConcreteOuter.nonzero_occupation_sparse_failure
#print axioms Spin.Structured.ConcreteOuter.good_nonzero_occupation_sparse_failure

end Spin.Structured.ConcreteOuter
