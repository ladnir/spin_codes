import SpinCodes.Structured.ConcreteTransfer

/-! Supplemental statement pins and axiom audit for the concrete fixed-weight step. -/

namespace Spin.Structured.ConcreteMaps
open Spin.Imt

example (j : Fin 129) (z : ℝ) (q : Finset (Fin 19)) :
    emittedMoment j z q = Occupation.hyperMoment (Aset q).card j z :=
  emittedMoment_eq_hyper j z q

example {z : ℝ} (hz : 0 ≤ z) (j : Fin 129) :
    LiveDominates actualShellSystem
      (Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowZ (stepRow j z ∅) :=
  stepRow_dominates_Z hz j

example {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (j : Fin 129)
    (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2) (c : Coords 5) (μ : Finset (Fin 19) → ℝ)
    (hc : c.Nonneg) (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).apply c) (fixedStep j z μ) :=
  fixedStep_liveDominates_generic hz hz1 j hj1 hj2 c μ hc hμ

end Spin.Structured.ConcreteMaps

#print axioms Spin.Structured.ConcreteMaps.emittedMoment_eq_hyper
#print axioms Spin.Structured.ConcreteMaps.cancellationMoment_le_cap
#print axioms Spin.Structured.ConcreteMaps.cancellationMoment_le_live
#print axioms Spin.Structured.ConcreteMaps.shellCancellationMoment_le_cap
#print axioms Spin.Structured.ConcreteMaps.stepRow_nonzero
#print axioms Spin.Structured.ConcreteMaps.stepRow_from_zero
#print axioms Spin.Structured.ConcreteMaps.refreshRow_zero_le_min
#print axioms Spin.Structured.ConcreteMaps.stepRow_dominates_point
#print axioms Spin.Structured.ConcreteMaps.lowPatterns_none_checked
#print axioms Spin.Structured.ConcreteMaps.stepRow_dominates_D_generic
#print axioms Spin.Structured.ConcreteMaps.stepRow_dominates_Z
#print axioms Spin.Structured.ConcreteMaps.live_zero_syndrome
#print axioms Spin.Structured.ConcreteMaps.shellStepRow_dominates_generic
#print axioms Spin.Imt.kernelApply_liveDominates
#print axioms Spin.Structured.ConcreteMaps.actual_fixed_nonneg
#print axioms Spin.Structured.ConcreteMaps.fixedStep_liveDominates_of_cancellation
#print axioms Spin.Structured.ConcreteMaps.fixedStep_liveDominates_generic
