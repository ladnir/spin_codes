import SpinCodes.Structured.ConcreteFixedSubsetMarginalDefs
import SpinCodes.Structured.ConcretePlacementSimplexSum
import SpinCodes.Structured.ConcreteShufflePermutation

noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

theorem uniform_expect_equiv {α β : Type*} [Fintype α] [Fintype β] [Nonempty α] [Nonempty β]
    (e : α≃β) (f : β→ℝ) :
    (FinPMF.uniform α).expect (fun x => f (e x))=(FinPMF.uniform β).expect f := by
  change (∑ x:α, (1/(Fintype.card α:ℝ))*f (e x)) = ∑ y:β, (1/(Fintype.card β:ℝ))*f y
  rw [Fintype.card_congr e]
  exact Fintype.sum_equiv e _ _ (fun _ => rfl)

theorem expect_comm {α β : Type*} [Fintype α] [Fintype β]
    (P : FinPMF α) (T : FinPMF β) (f : α→β→ℝ) :
    P.expect (fun x => T.expect (f x))=T.expect (fun y => P.expect (fun x => f x y)) := by
  simp only [FinPMF.expect,mul_sum]
  rw [sum_comm]
  apply sum_congr rfl
  intro y _
  apply sum_congr rfl
  intro x _
  ring

/-- A marked subset of a uniform injection is itself uniform on its cardinality layer. -/
theorem uniform_embedding_subset_expect {R Q : ℕ} (hQR : Q≤R)
    (A : Finset (Fin Q)) (F : Finset (Fin R)→ℝ) :
    letI : Nonempty (Fin Q ↪ Fin R) := ⟨Fin.castLEEmb hQR⟩
    (FinPMF.uniform (Fin Q ↪ Fin R)).expect (fun e => F (A.map e)) =
      (shuffleLaw (A.map (Fin.castLEEmb hQR))).expect F := by
  letI : Nonempty (Fin Q ↪ Fin R) := ⟨Fin.castLEEmb hQR⟩
  let U := FinPMF.uniform (Fin Q ↪ Fin R)
  let V := FinPMF.uniform (Equiv.Perm (Fin R))
  have hi (σ : Equiv.Perm (Fin R)) :
      U.expect (fun e => F (shuffleSupport σ (A.map e)))=U.expect (fun e => F (A.map e)) := by
    have he := uniform_expect_equiv (Equiv.embeddingCongr (Equiv.refl (Fin Q)) σ.symm)
      (fun e : Fin Q ↪ Fin R => F (A.map e))
    have heq (e : Fin Q ↪ Fin R) :
        Equiv.embeddingCongr (Equiv.refl (Fin Q)) σ.symm e=e.trans σ.symm.toEmbedding := by
      ext i
      rfl
    simpa only [U,shuffleSupport,Finset.map_map,heq] using he
  calc
    _ = V.expect (fun σ => U.expect (fun e => F (shuffleSupport σ (A.map e)))) := by
      simp only [hi,FinPMF.expect_const]
      rfl
    _ = U.expect (fun e => V.expect (fun σ => F (shuffleSupport σ (A.map e)))) :=
      expect_comm V U _
    _ = U.expect (fun _ => (shuffleLaw (A.map (Fin.castLEEmb hQR))).expect F) := by
      congr 1
      funext e
      rw [←FinPMF.expect_map (FinPMF.uniform (Equiv.Perm (Fin R)))
        (fun σ => shuffleSupport σ (A.map e)) F,uniform_permutation_shuffleLaw]
      congr 1
      apply FinPMF.ext
      funext T
      simp only [shuffleLaw_apply,card_map]
    _ = _ := FinPMF.expect_const _ _

#print axioms uniform_embedding_subset_expect
end Spin.Structured.Placement
