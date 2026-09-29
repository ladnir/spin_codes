import SpinCodes.Structured.ConcreteMarkedExperiment
import SpinCodes.Structured.ConcreteMarkedReferenceMoment

open Spin Spin.Structured Spin.Structured.ConcreteMarked

example {L : ℕ} (S : Finset (Fin L)) (b : ℕ) :
    (fairRows S b).bind ConcreteRoute.law = regionsLaw S b := fairRows_route S b

example {L : ℕ} (S : Finset (Fin L)) (p : ℝ) (hp0 : 0 ≤ p) (hp1 : p ≤ 1)
    (h : 0 < (markBitLaw L p hp0 hp1).prob (fun x => x.1.card = S.card)) :
    ((markBitLaw L p hp0 hp1).condition (fun x => x.1.card = S.card) h).map
      (fun x => x.1 ∩ x.2) = regionLaw S := markBitLaw_conditioned S p hp0 hp1 h

example {L b R : ℕ} (S : Finset (Fin L))
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128))
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) (d : ℕ) :
    (fairExperimentLaw S b R).prob
      (fun ω => ConcreteRoutedEncoder.weight e ω.1 ω.2 ≤ d) ≤
        (((1 / markMass L S.card ((8/5)*α)) ^ b) *
          (2048 * (1 - 96*α) ^ R)) / ((1 - (8/5)*α) ^ d) :=
  fairExperiment_sparse_probability S e hα0 hα1 d

#print axioms iid_thinning
#print axioms conditioned_marks_eq_shuffleLaw
#print axioms regionLaw_dominates
#print axioms regionsLaw_dominates
#print axioms shuffled_fairOnMarks
#print axioms fairRows_route
#print axioms markBitLaw_unconditioned
#print axioms markBitLaw_conditioned
#print axioms sparse_moment
#print axioms fairExperiment_sparse_moment
#print axioms fairExperiment_sparse_probability
