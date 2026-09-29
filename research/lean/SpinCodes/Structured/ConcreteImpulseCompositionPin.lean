import SpinCodes.Structured.ConcreteImpulseComposition

open Spin.Structured.ConcreteEncoder

example {m n : Nat} (z : ℝ) (xs : Fin m → Input) (ys : Fin n → Input) :
    endpointKernel z (Fin.append xs ys) = endpointKernel z xs * endpointKernel z ys :=
  endpointKernel_append z xs ys

example {a : Nat} (z : ℝ) (gaps : Fin a → Nat) (last : Nat) :
    matrixExpect (Spin.piPMF (fun _ : Fin a => coordinateLaw))
      (fun coords => endpointKernel z (impulseInputs gaps coords last).get) =
        uniformImpulseProduct z gaps last :=
  endpointKernel_uniform_impulses z gaps last

#print axioms Spin.Structured.ConcreteEncoder.terminalState_append
#print axioms Spin.Structured.ConcreteEncoder.outputWeight_append
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_eq_path
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_append
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_empty
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_total
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_impulseInputs
#print axioms Spin.Structured.ConcreteEncoder.uniformImpulseKernel_apply
#print axioms Spin.Structured.ConcreteEncoder.impulseProduct_average
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_uniform_impulses
