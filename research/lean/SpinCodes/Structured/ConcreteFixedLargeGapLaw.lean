import SpinCodes.Structured.ConcreteFixedLargeGapPermutation

/-! Exact finite permutation symmetry, before any continuum limit. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset

private theorem permuted_gap_sum {R Q : ℕ} (π : Equiv.Perm (Fin (Q+1))) (B : BlockSubset R Q) :
    (∑ i, emptyGaps B (π i))+Q = R := by
  rw [Equiv.sum_comp π (emptyGaps B)]
  exact emptyGaps_sum B

def gapPermutation {R Q : ℕ} (π : Equiv.Perm (Fin (Q+1))) : BlockSubset R Q ≃ BlockSubset R Q where
  toFun B := gapBlocks (fun i => emptyGaps B (π i)) (permuted_gap_sum π B)
  invFun B := gapBlocks (fun i => emptyGaps B (π.symm i)) (permuted_gap_sum π.symm B)
  left_inv B := by
    apply emptyGaps_injective
    funext i
    simp only [emptyGaps_gapBlocks, Equiv.apply_symm_apply]
  right_inv B := by
    apply emptyGaps_injective
    funext i
    simp only [emptyGaps_gapBlocks, Equiv.symm_apply_apply]

@[simp] theorem emptyGaps_gapPermutation {R Q : ℕ} (π : Equiv.Perm (Fin (Q+1)))
    (B : BlockSubset R Q) (i : Fin (Q+1)) :
    emptyGaps (gapPermutation π B) i = emptyGaps B (π i) := by
  exact congrFun (emptyGaps_gapBlocks _ _) i

/-- The actual uniform ordered-block law has exchangeable integer empty gaps. -/
theorem sum_emptyGaps_permutation {R Q : ℕ} (π : Equiv.Perm (Fin (Q+1)))
    (F : (Fin (Q+1) → ℕ) → ℝ) :
    (∑ B : BlockSubset R Q, F (fun i => emptyGaps B (π i))) =
      ∑ B : BlockSubset R Q, F (emptyGaps B) := by
  have h := Equiv.sum_comp (gapPermutation (R := R) π) (fun B => F (emptyGaps B))
  simpa only [funext (emptyGaps_gapPermutation π _)] using h

#print axioms gapPermutation
#print axioms sum_emptyGaps_permutation
end Spin.Structured.Placement
