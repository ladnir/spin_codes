import SpinCodes.Structured.ConcreteOccupation

/-! Supplemental pins: all weights, actual Bernoulli input, and the explicit
sparse finite-round bound. The original distance theorem pins are unchanged. -/
namespace Spin.Structured.ConcreteMaps
open Spin.Imt

example {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (j : Fin 129)
    (c : Coords 5) (μ : Finset (Fin 19) → ℝ)
    (hc : c.Nonneg) (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).apply c) (fixedStep j z μ) :=
  fixedStep_liveDominates hz hz1 j c μ hc hμ

example {α : ℝ} (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : Nat) :
    ∑ q, ((bernoulliStep ((4 / 5) * α) (1 - (8 / 5) * α))^[R]
      (fun q => if q = ∅ then (1 : ℝ) else 0)) q ≤ 2048 * (1 - 96 * α) ^ R :=
  sparse_bernoulli_moment h0 h1 R

end Spin.Structured.ConcreteMaps

#print axioms Spin.Structured.LowCancellation.Data.groups1_valid
#print axioms Spin.Structured.LowCancellation.Data.groups2_valid
#print axioms Spin.Structured.LowCancellation.Data.groups1_inputs
#print axioms Spin.Structured.LowCancellation.Data.groups2_inputs
#print axioms Spin.Structured.LowCancellation.inputs_layer
#print axioms Spin.Structured.ConcreteMaps.cancellationMoment_le_patterns_one
#print axioms Spin.Structured.ConcreteMaps.cancellationMoment_le_patterns_two
#print axioms Spin.Structured.ConcreteMaps.shellCancellationMoment_one
#print axioms Spin.Structured.ConcreteMaps.shellCancellationMoment_two
#print axioms Spin.Structured.ConcreteMaps.cancellationMoment_le_rowD_all
#print axioms Spin.Structured.ConcreteMaps.shellCancellationMoment_le_rowS_all
#print axioms Spin.Structured.ConcreteMaps.fixedStep_liveDominates
#print axioms Spin.Structured.ConcreteMaps.bernoulliStep_eq_mixture
#print axioms Spin.Structured.ConcreteMaps.bernoulliStep_liveDominates
#print axioms Spin.Structured.ConcreteMaps.bernoulli_moment_bound
#print axioms Spin.Structured.ConcreteMaps.sparse_bernoulli_moment
