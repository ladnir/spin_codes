import SpinCodes.Structured.ConcreteFixedInsertionOrdering

/-! Full-site products, permuted fugacities, and exact selected block coordinates. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder
attribute [local instance] Classical.propDecidable

def fullSiteProduct {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (x u : Fin Q → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  fullGapProduct γ P (fun i => positionGaps x i.castSucc) u (positionGaps x (Fin.last Q))

theorem continuous_fullSiteProduct {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (u : Fin Q → ℝ) : Continuous (fun x => fullSiteProduct γ P x u) := by
  unfold fullSiteProduct fullGapProduct
  simp_rw [List.ofFn_eq_map]
  apply Continuous.mul
  · apply continuous_list_prod
    intro i hi
    exact ((continuous_timeEmpty γ).comp (continuous_positionGaps i.castSucc)).mul continuous_const
  · exact (continuous_timeEmpty γ).comp (continuous_positionGaps (Fin.last Q))

theorem markCoefficient_permuted {Q : Nat} (u : Fin Q → ℝ) (π : Equiv.Perm (Fin Q)) (A : Finset (Fin Q)) :
    markCoefficient (fun i => u (π.symm i)) (A.map π.toEmbedding) = markCoefficient u A := by
  simp [markCoefficient]

theorem fullSiteProduct_permuted {Q : Nat} (γ : ℝ) (P : Matrix (Fin 2) (Fin 2) ℝ)
    (x u : Fin Q → ℝ) (π : Equiv.Perm (Fin Q)) :
    fullSiteProduct γ P x (fun i => u (π.symm i)) =
      ∑ A : Finset (Fin Q), markCoefficient u A •
        selectedGapProduct γ P (fun i => positionGaps x i.castSucc) (positionGaps x (Fin.last Q)) (A.map π.toEmbedding) := by
  rw [fullSiteProduct, fullGapProduct_subset_expansion, ← Equiv.sum_comp (Equiv.finsetCongr π)]
  simp only [Equiv.finsetCongr_apply, markCoefficient_permuted]

def selectedBlocks {R Q a : Nat} (B : BlockSubset R Q) (A : Finset (Fin Q)) (hA : A.card = a) : BlockSubset R a :=
  ⟨A.map (orderedBlocks B), by rw [card_map,hA]⟩

theorem orderedBlocks_selectedBlocks {R Q a : Nat} (B : BlockSubset R Q) (A : Finset (Fin Q)) (hA : A.card = a) :
    (orderedBlocks (selectedBlocks B A hA) : Fin a → Fin R) =
      fun i => orderedBlocks B (A.orderEmbOfFin hA i) := by
  symm
  apply Finset.orderEmbOfFin_unique
  · intro i
    exact mem_map.mpr ⟨_,A.orderEmbOfFin_mem hA i,rfl⟩
  · exact (orderedBlocks_strictMono B).comp (A.orderEmbOfFin hA).strictMono

theorem normalizedSites_selectedBlocks {R Q a : Nat} (B : BlockSubset R Q) (A : Finset (Fin Q)) (hA : A.card = a) :
    normalizedSites (selectedBlocks B A hA) = fun i => normalizedSites B (A.orderEmbOfFin hA i) := by
  funext i
  simp only [normalizedSites, orderedBlocks_selectedBlocks]

theorem selectedGapProduct_normalized {R Q a : Nat} (γ : ℝ) (B : BlockSubset R Q)
    (A : Finset (Fin Q)) (hA : A.card = a) :
    selectedGapProduct γ impulseMatrix (fun i => positionGaps (normalizedSites B) i.castSucc)
      (positionGaps (normalizedSites B) (Fin.last Q)) A =
        siteProduct γ (normalizedSites (selectedBlocks B A hA)) := by
  subst a
  rw [selectedGapProduct_eq_subtuple, normalizedSites_selectedBlocks]

theorem fullSiteProduct_normalized_permuted {R Q : Nat} (γ : ℝ) (B : BlockSubset R Q)
    (u : Fin Q → ℝ) (π : Equiv.Perm (Fin Q)) :
    fullSiteProduct γ impulseMatrix (normalizedSites B) (fun i => u (π.symm i)) =
      ∑ A : Finset (Fin Q), markCoefficient u A •
        siteProduct γ (normalizedSites (selectedBlocks B (A.map π.toEmbedding) (card_map _))) := by
  rw [fullSiteProduct_permuted]
  apply sum_congr rfl
  intro A hA
  rw [selectedGapProduct_normalized γ B (A.map π.toEmbedding) (card_map _)]

end Spin.Structured.Placement
