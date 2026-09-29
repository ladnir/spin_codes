import SpinCodes.Structured.ConcreteOuterEnvelope

namespace Spin.Structured.ConcreteOuter

example {n k : ℕ} (hk : k ≤ n) :
    Real.exp ((n : ℝ) * Spin.Numeric.hEnt ((k : ℝ)/n)) / ((n : ℝ)+1) ≤
      (n.choose k : ℝ) ∧
    (n.choose k : ℝ) ≤ Real.exp ((n : ℝ) * Spin.Numeric.hEnt ((k : ℝ)/n)) :=
  ⟨exp_entropy_div_le_choose hk, choose_le_exp_entropy hk⟩

example {k : ℕ} (hk : 0 < k) (w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w : ℝ)) ≤
      Real.exp (finiteEnvelope k w + 4 * Real.log (((k * 24 : ℕ) : ℝ) + 1)) :=
  expected_spectrum_explicit_log_bound hk w

#print axioms Spin.Structured.ConcreteOuter.type_power_eq_exp
#print axioms Spin.Structured.ConcreteOuter.choose_le_exp_entropy
#print axioms Spin.Structured.ConcreteOuter.exp_entropy_div_le_choose
#print axioms Spin.Structured.ConcreteOuter.golayCoefficient_le_exp_gBA
#print axioms Spin.Structured.ConcreteOuter.transition_le_exp_entropy
#print axioms Spin.Structured.ConcreteOuter.expected_spectrum_finite_bound
#print axioms Spin.Structured.ConcreteOuter.expected_spectrum_explicit_log_bound

end Spin.Structured.ConcreteOuter
