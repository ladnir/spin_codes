import SpinCodes.Structured.ConcretePlacementSimplexLaw
import Mathlib.Data.Fintype.CardEmbedding

noncomputable section
namespace Spin.Structured.Placement
open Finset
attribute [local instance] Classical.propDecidable

def labelledEmbedding {R Q : ℕ} (B : BlockSubset R Q) (π : Equiv.Perm (Fin Q)) : Fin Q ↪ Fin R :=
  π.toEmbedding.trans (orderedBlocks B)

theorem labelledEmbedding_image {R Q : ℕ} (B : BlockSubset R Q) (π : Equiv.Perm (Fin Q)) :
    univ.map (labelledEmbedding B π)=B.val := by
  rw [labelledEmbedding,←Finset.map_map,Finset.univ_map_equiv_to_embedding]
  exact Finset.map_orderEmbOfFin_univ B.val B.property

theorem labelledEmbedding_bijective (R Q : ℕ) :
    Function.Bijective (fun p : BlockSubset R Q × Equiv.Perm (Fin Q) => labelledEmbedding p.1 p.2) := by
  constructor
  · rintro ⟨B,π⟩ ⟨C,ρ⟩ he
    change labelledEmbedding B π = labelledEmbedding C ρ at he
    have hB : B=C := by
      apply Subtype.ext
      rw [←labelledEmbedding_image B π,←labelledEmbedding_image C ρ,he]
    subst C
    have hp : π=ρ := by
      apply Equiv.ext
      intro i
      apply (orderedBlocks B).injective
      exact congrArg (fun e : Fin Q ↪ Fin R => e i) he
    simp only [hp]
  · intro e
    let B : BlockSubset R Q := ⟨univ.map e, by simp⟩
    let f : Fin Q ↪ Fin Q := {
      toFun := fun i => (B.val.orderIsoOfFin B.property).symm
        ⟨e i, Finset.mem_map.mpr ⟨i,mem_univ _,rfl⟩⟩
      inj' := by
        intro i j hij
        apply e.injective
        have hh := congrArg (B.val.orderIsoOfFin B.property) hij
        simpa using congrArg Subtype.val hh }
    let π := f.equivOfFiniteSelfEmbedding
    refine ⟨(B,π),?_⟩
    apply Function.Embedding.ext
    intro i
    change (B.val.orderIsoOfFin B.property ((B.val.orderIsoOfFin B.property).symm
      ⟨e i,?_⟩)).val=e i
    simp

/-- A uniform ordered Q-subset and uniform labelling enumerate every Q-site injection once. -/
def labelledEmbeddingEquiv (R Q : ℕ) : (BlockSubset R Q × Equiv.Perm (Fin Q)) ≃ (Fin Q ↪ Fin R) :=
  Equiv.ofBijective _ (labelledEmbedding_bijective R Q)

#print axioms labelledEmbedding_bijective
end Spin.Structured.Placement
