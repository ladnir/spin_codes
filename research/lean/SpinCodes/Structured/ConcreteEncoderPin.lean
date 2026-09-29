import SpinCodes.Structured.ConcreteEncoderMoment

open Spin.Structured.ConcreteEncoder Spin.Structured.ConcreteMaps

example (P : Spin.FinPMF (Finset (Fin 128))) (β z : ℝ)
    (hP : ∀ x, P.p x = β ^ x.card * (1 - β) ^ (128 - x.card)) (R : Nat) :
    (Spin.piPMF (fun _ : Fin R => P)).expect (fun xs => inputMoment z xs ∅) =
      ∑ q, ((bernoulliStep β z)^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q :=
  averagedMoment_eq_iterate P β z hP R

example (P : Spin.FinPMF (Finset (Fin 128))) {α : ℝ}
    (hP : ∀ x, P.p x = ((4 / 5) * α) ^ x.card *
      (1 - (4 / 5) * α) ^ (128 - x.card))
    (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : Nat) :
    (experimentLaw P R).expect
        (fun ω => (1 - (8 / 5) * α) ^ outputWeight ω.1 ω.2 ∅) ≤
      2048 * (1 - 96 * α) ^ R :=
  sparse_encoder_experiment P hP h0 h1 R

#print axioms Spin.Structured.ConcreteEncoder.inputMoment_nonneg
#print axioms Spin.Structured.ConcreteEncoder.transvection_expect
#print axioms Spin.Structured.ConcreteEncoder.nextState_expect
#print axioms Spin.Structured.ConcreteEncoder.averagedMoment_eq_iterate
#print axioms Spin.Structured.ConcreteEncoder.encoder_moment_bound
#print axioms Spin.Structured.ConcreteEncoder.sparse_encoder_moment
#print axioms Spin.Structured.ConcreteEncoder.sparse_encoder_experiment
