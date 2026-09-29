import SpinCodes.Structured.ConcreteClosedShells
import SpinCodes.Structured.SparseNonneg

/-! One-step domination for the actual encoder at each fixed input weight
other than the two weights that use special cancellation tables. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Spin.Imt

def fixedStep (j : ℕ) (z : ℝ) :
    (Finset (Fin 19) → ℝ) → (Finset (Fin 19) → ℝ) := kernelApply (stepRow j z)

theorem actual_fixed_nonneg {z : ℝ} (hz : 0 ≤ z) (j : Fin 129) :
    (Occupation.fixed Occupation.Sparse.count 524287
      (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).Nonneg := by
  apply Occupation.Sparse.fixed_row_nonneg (Nat.le_of_lt_succ j.isLt) _ _ hz
  rw [polyChoose_eq, ← kernel_card_frozen j]
  exact kernelLayer_card_le j

theorem fixedStep_liveDominates_of_cancellation {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129)
    (hD : ∀ q ≠ ∅, cancellationMoment j z q ≤
      (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelD)
    (hS : ∀ i, shellCancellationMoment j z i ≤
      (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i)
    (c : Coords 5) (μ : Finset (Fin 19) → ℝ) (hc : c.Nonneg)
    (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).apply c) (fixedStep j z μ) := by
  apply kernelApply_liveDominates actualShellSystem _ (stepRow j z)
    (stepRow_nonneg hz j) (actual_fixed_nonneg hz j).2.1 (stepRow_dominates_Z hz j)
    (fun q hq => stepRow_dominates_D hz hz1 j hq (hD q hq)) _ c μ hc hμ
  intro i
  change LiveDominates actualShellSystem _
    (fun r => (∑ q ∈ weightShell i, stepRow j z q r) / ((weightShell i).card : ℝ))
  rw [weightShell_card_actual]
  by_cases hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0
  · exact shellStepRow_dominates_closed hz j i hl (hS i)
  · exact shellStepRow_dominates_live hz hz1 j i hl (hS i)

theorem fixedStep_liveDominates_generic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2)
    (c : Coords 5) (μ : Finset (Fin 19) → ℝ) (hc : c.Nonneg)
    (hμ : LiveDominates actualShellSystem c μ) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).apply c) (fixedStep j z μ) := by
  apply kernelApply_liveDominates actualShellSystem _ (stepRow j z)
    (stepRow_nonneg hz j) (actual_fixed_nonneg hz j).2.1
    (stepRow_dominates_Z hz j)
    (fun q hq => stepRow_dominates_D_generic hz hz1 j hj1 hj2 hq) _ c μ hc hμ
  intro i
  change LiveDominates actualShellSystem _
    (fun r => (∑ q ∈ weightShell i, stepRow j z q r) / ((weightShell i).card : ℝ))
  rw [weightShell_card_actual]
  exact shellStepRow_dominates_generic hz hz1 j hj1 hj2 i

end Spin.Structured.ConcreteMaps
