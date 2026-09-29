import SpinCodes.Structured.ConcreteImpulseCoarseProductUniform

open Spin.Structured.ConcreteEncoder

example {a : Nat} {θ L : ℝ} (hθ : 0 ≤ θ) (hL : 0 < L)
    (gaps : Fin a → Nat) (last : Nat) (q r : State) :
    |(matrixExpect (Spin.piPMF (fun _ : Fin a => coordinateLaw))
        (fun coords => endpointKernel (Real.exp (-(θ / L))) (impulseInputs gaps coords last).get)) q r -
      (liveLift * coarseImpulseProduct (θ / L) gaps last * liveProjection) q r| ≤
        impulseErrorBudget (θ / L) gaps last :=
  endpointKernel_uniform_impulses_coarse_error hθ hL gaps last q r

#print axioms Spin.Structured.ConcreteEncoder.liveProjection_lift
#print axioms Spin.Structured.ConcreteEncoder.uniformImpulseKernel_one_coarse
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_lifted_error
#print axioms Spin.Structured.FiniteKernel.product_rowError
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_substochastic
#print axioms Spin.Structured.ConcreteEncoder.emptyKernel_rowError
#print axioms Spin.Structured.ConcreteEncoder.roundKernel_exp_rowError
#print axioms Spin.Structured.ConcreteEncoder.uniformImpulseKernel_exp_rowError
#print axioms Spin.Structured.ConcreteEncoder.approximateImpulseProduct_eq_lift
#print axioms Spin.Structured.ConcreteEncoder.impulseErrorBudget_sum
#print axioms Spin.Structured.ConcreteEncoder.uniformImpulseProduct_substochastic
#print axioms Spin.Structured.ConcreteEncoder.uniformImpulseProduct_rowError
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_uniform_impulses_coarse_error
#print axioms Spin.Structured.ConcreteEncoder.impulseErrorBudget_gap_bound
#print axioms Spin.Structured.ConcreteEncoder.endpointKernel_uniform_impulses_good_gaps
