import SpinCodes.Structured.ConcreteEpochMean
import SpinCodes.Structured.ConcreteEpochKernel
import SpinCodes.Structured.ConcreteImpulse

open Spin.Structured.ConcreteEncoder

example (g : Nat) (f : Finset (Fin 19) → ℝ) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect (fun ts => f (emptyState ts q)) =
      liveMean f + (1 / 2 : ℝ) ^ g * (f q - liveMean f) :=
  emptyExpect_live g f hq

example (g : Nat) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    |(Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect
        (fun ts => (outputWeight (fun _ => ∅) ts q : ℝ)) -
      (128 * (262144 / 524287 : ℝ)) * g| ≤ 256 :=
  emptyOutputMean_bound g hq

example : uniformImpulseCoarse = !![(0 : ℝ), 1; 1 / 524287, 1 - 1 / 524287] :=
  uniformImpulseCoarse_eq

#print axioms Spin.Structured.ConcreteEncoder.transvection_nonzero
#print axioms Spin.Structured.ConcreteEncoder.emptyState_nonzero
#print axioms Spin.Structured.ConcreteEncoder.emptyExpect_live
#print axioms Spin.Structured.ConcreteEncoder.emitted_liveMean
#print axioms Spin.Structured.ConcreteEncoder.emptyOutputMean_bound
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_one_live
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_one_error
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_live_zero
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_zero
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_total
#print axioms Spin.Structured.ConcreteEncoder.impulse_from_zero
#print axioms Spin.Structured.ConcreteEncoder.impulse_from_live
#print axioms Spin.Structured.ConcreteEncoder.uniformImpulseCoarse_eq
