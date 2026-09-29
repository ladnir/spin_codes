import SpinCodes.Structured.ConcreteFixedSubsetMarginalUniform

noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

theorem uniform_expect_sum {α : Type*} [Fintype α] [Nonempty α] (f : α→ℝ) :
    (FinPMF.uniform α).expect f=(∑ x, f x)/(Fintype.card α:ℝ) := by
  change (∑ x, (1/(Fintype.card α:ℝ))*f x)=_
  rw [←mul_sum]
  ring

theorem shuffleLaw_expect_layer {R a : ℕ} (S : Finset (Fin R)) (hS : S.card=a)
    (F : Finset (Fin R)→ℝ) :
    (shuffleLaw S).expect F=(∑ C:BlockSubset R a, F C.val)/(R.choose a:ℝ) := by
  simp only [FinPMF.expect,shuffleLaw_apply,hS,ite_mul,zero_mul]
  rw [←sum_filter]
  simp only [one_div,mul_comm ((R.choose a:ℝ)⁻¹),←sum_mul]
  congr 1
  exact Finset.sum_subtype _ (by intro T; simp) F

/-- Picking any labelled subcollection from uniform Q distinct sites gives uniform distinct sites. -/
theorem picked_subset_average {R Q : ℕ} (hQR : Q≤R) (A : Finset (Fin Q))
    (F : Finset (Fin R)→ℝ) :
    (∑ B:BlockSubset R Q, ∑ π:Equiv.Perm (Fin Q),
      F (A.map (labelledEmbedding B π)))/((R.choose Q:ℝ)*(Q.factorial:ℝ)) =
      (∑ C:BlockSubset R A.card,F C.val)/(R.choose A.card:ℝ) := by
  let e0 : Fin Q ↪ Fin R := Fin.castLEEmb hQR
  letI : Nonempty (Fin Q ↪ Fin R) := ⟨e0⟩
  letI : Nonempty (BlockSubset R Q) := ⟨⟨univ.map e0,by simp⟩⟩
  have hh := (uniform_expect_equiv (labelledEmbeddingEquiv R Q)
    (fun e : Fin Q ↪ Fin R => F (A.map e))).trans
      (uniform_embedding_subset_expect hQR A F)
  rw [uniform_expect_sum,Fintype.card_prod,card_blockSubset,Fintype.card_perm,Fintype.card_fin,
    Nat.cast_mul,Fintype.sum_prod_type,
    shuffleLaw_expect_layer _ (Finset.card_map _) F] at hh
  simpa only [labelledEmbeddingEquiv,Equiv.ofBijective_apply] using hh

def picked {R Q : ℕ} (A : Finset (Fin Q)) (B : BlockSubset R Q) (π : Equiv.Perm (Fin Q)) :
    BlockSubset R A.card := ⟨A.map (labelledEmbedding B π),Finset.card_map _⟩

theorem picked_blockSubset_average {R Q : ℕ} (hQR : Q≤R) (A : Finset (Fin Q))
    (φ : BlockSubset R A.card→ℝ) :
    (∑ B:BlockSubset R Q, ∑ π:Equiv.Perm (Fin Q),φ (picked A B π))/
      ((R.choose Q:ℝ)*(Q.factorial:ℝ)) = (∑ C:BlockSubset R A.card,φ C)/(R.choose A.card:ℝ) := by
  let F : Finset (Fin R)→ℝ := fun T => if h:T.card=A.card then φ ⟨T,h⟩ else 0
  have h := picked_subset_average hQR A F
  have hleft (B : BlockSubset R Q) (π : Equiv.Perm (Fin Q)) :
      F (A.map (labelledEmbedding B π))=φ (picked A B π) := by simp [F,picked]
  have hright (C : BlockSubset R A.card) : F C.val=φ C := by simp [F,C.property]
  simpa only [hleft,hright] using h

#print axioms picked_subset_average
#print axioms picked_blockSubset_average
end Spin.Structured.Placement
