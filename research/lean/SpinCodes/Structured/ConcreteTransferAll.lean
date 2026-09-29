import SpinCodes.Structured.ConcreteLowCancellation
import SpinCodes.Structured.ConcreteTransfer

/-! Actual fixed-input-weight transfer, with no exceptional weights or
unproved cancellation hypotheses. -/
noncomputable section
namespace Spin.Structured.ConcreteMaps
open Spin.Imt

theorem fixedStep_liveDominates {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (c : Coords 5) (μ : Finset (Fin 19) → ℝ) (hc : c.Nonneg)
    (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).apply c) (fixedStep j z μ) :=
  fixedStep_liveDominates_of_cancellation hz hz1 j
    (fun _ hq => cancellationMoment_le_rowD_all hz hz1 j hq)
    (shellCancellationMoment_le_rowS_all hz hz1 j) c μ hc hμ

end Spin.Structured.ConcreteMaps
