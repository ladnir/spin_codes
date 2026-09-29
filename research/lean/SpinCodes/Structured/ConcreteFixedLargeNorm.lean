import SpinCodes.Structured.ConcreteFixedLargeLaplaceActual
import SpinCodes.Structured.ConcreteFixedInsertionNorm

/-! The paper's large-occupation weighted norm for the actual continuum encoder. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder

def paperLargeTheta (Q : ℕ) : ℝ := (Q:ℝ)*(133/125)/(262144/524287)

theorem paperLargeTheta_pos {Q : ℕ} (hQ : 0 < Q) : 0 < paperLargeTheta Q := by
  unfold paperLargeTheta
  positivity

theorem paperLargeTheta_epochMean (Q : ℕ) :
    paperLargeTheta Q * epochMean / 128 = (Q:ℝ)*(133/125) := by
  unfold paperLargeTheta epochMean
  ring

theorem fugacityContinuum_paper_rowNorm {Q : ℕ} (hQ : 3 ≤ Q)
    (u : Fin Q → ℝ) (hu : ∀ i, 0 ≤ u i) :
    Spin.RowNormLe (3/1600)
      (Real.exp ((Q:ℝ)*paperLargeExponent) * ∏ i, Spin.mv (127/250) (3/1600) (1/524287) (u i))
      (fugacityContinuum (paperLargeTheta Q) u) := by
  apply fugacityContinuum_rowNorm_of_spacing (by norm_num) (by norm_num) (by norm_num)
    (Real.exp_nonneg _) ?_ u hu
  intro S
  simpa only [paperLargeTheta_epochMean, neg_mul] using paper_liveDuration_laplace_le hQ S

#print axioms paperLargeTheta_pos
#print axioms fugacityContinuum_paper_rowNorm
end Spin.Structured.Placement
