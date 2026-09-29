import SpinCodes.Structured.ConcreteFixedInsertionSites
import SpinCodes.Structured.ConcreteFixedSubsetMarginal

/-! Exact finite insertion marginal, before taking a continuum limit. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder
attribute [local instance] Classical.propDecidable

theorem selectedBlocks_permuted_eq_picked {R Q : Nat} (A : Finset (Fin Q))
    (B : BlockSubset R Q) (π : Equiv.Perm (Fin Q)) :
    selectedBlocks B (A.map π.toEmbedding) (card_map _) = picked A B π := by
  apply Subtype.ext
  simp only [selectedBlocks, picked, labelledEmbedding, Finset.map_map]

theorem fullSiteProduct_normalized_picked {R Q : Nat} (γ : ℝ) (B : BlockSubset R Q)
    (u : Fin Q → ℝ) (π : Equiv.Perm (Fin Q)) :
    fullSiteProduct γ impulseMatrix (normalizedSites B) (fun i => u (π.symm i)) =
      ∑ A : Finset (Fin Q), markCoefficient u A • siteProduct γ (normalizedSites (picked A B π)) := by
  rw [fullSiteProduct_normalized_permuted]
  simp only [selectedBlocks_permuted_eq_picked]

def finiteInsertionAverage (R : Nat) {Q : Nat} (θ : ℝ) (u : Fin Q → ℝ) (i j : Fin 2) : ℝ :=
  (∑ B : BlockSubset R Q, ∑ π : Equiv.Perm (Fin Q),
    fullSiteProduct (θ * epochMean / 128) impulseMatrix (normalizedSites B) (fun k => u (π.symm k)) i j) /
      ((R.choose Q : ℝ) * (Q.factorial : ℝ))

theorem finiteInsertionAverage_eq {R Q : Nat} (hQR : Q ≤ R) (θ : ℝ) (u : Fin Q → ℝ) (i j : Fin 2) :
    finiteInsertionAverage R θ u i j = ∑ A : Finset (Fin Q), markCoefficient u A *
      ((∑ C : BlockSubset R A.card, siteProduct (θ * epochMean / 128) (normalizedSites C) i j) /
        (R.choose A.card : ℝ)) := by
  unfold finiteInsertionAverage
  simp only [fullSiteProduct_normalized_picked, Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]
  have hs : (∑ B : BlockSubset R Q, ∑ π : Equiv.Perm (Fin Q), ∑ A : Finset (Fin Q),
      markCoefficient u A * siteProduct (θ * epochMean / 128) (normalizedSites (picked A B π)) i j) =
      ∑ A : Finset (Fin Q), ∑ B : BlockSubset R Q, ∑ π : Equiv.Perm (Fin Q),
        markCoefficient u A * siteProduct (θ * epochMean / 128) (normalizedSites (picked A B π)) i j := by
    calc
      _ = ∑ B : BlockSubset R Q, ∑ A : Finset (Fin Q), ∑ π : Equiv.Perm (Fin Q),
          markCoefficient u A * siteProduct (θ * epochMean / 128) (normalizedSites (picked A B π)) i j := by
        apply sum_congr rfl
        intro B hB
        exact sum_comm
      _ = _ := sum_comm
  rw [hs, sum_div]
  apply sum_congr rfl
  intro A hA
  simp only [← mul_sum, mul_div_assoc]
  exact congrArg (fun z => markCoefficient u A * z)
    (picked_blockSubset_average hQR A (fun C => siteProduct (θ * (epochMean / 128)) (normalizedSites C) i j))

end Spin.Structured.Placement

