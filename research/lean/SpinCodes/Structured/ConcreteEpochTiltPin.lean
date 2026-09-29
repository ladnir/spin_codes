import SpinCodes.Structured.ConcreteEpochTilt

open Spin.Structured.ConcreteEncoder

example (g : Nat) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect
      (fun ts => ((outputWeight (fun _ => ∅) ts q : ℝ) - epochMean * g) ^ 2) ≤
        3 * 128 ^ 2 * g :=
  emptyCenteredSecond_bound g hq

example {θ L : ℝ} (hθ : 0 ≤ θ) (hL : 0 < L) (g : Nat)
    {q r : Finset (Fin 19)} (hq : q ≠ ∅) (hr : r ≠ ∅) :
    |emptyKernel (Real.exp (-(θ / L))) g q r -
      Real.exp (-(θ / L) * (epochMean * g)) / 524287| ≤
      θ / L * (128 * Real.sqrt (3 * g)) + (1 / 2 : ℝ) ^ g :=
  emptyKernel_tilt_scaled hθ hL g hq hr

#print axioms Spin.Structured.ConcreteEncoder.emptyCenteredSecond_succ
#print axioms Spin.Structured.ConcreteEncoder.emptyCenteredSecond_bound
#print axioms Spin.Structured.ConcreteEncoder.emptyVariance_bound
#print axioms Spin.Structured.ConcreteEncoder.emptyAbsDeviation_bound
#print axioms Spin.Structured.ConcreteEncoder.exp_neg_lipschitz
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_tilt_deviation
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_tilt_live
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_tilt_scaled
