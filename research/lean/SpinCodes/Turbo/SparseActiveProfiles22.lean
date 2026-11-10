import SpinCodes.Turbo.SparseSupportPrune22
import SpinCodes.Turbo.SparseSupportArray22

/-! Reindex a fixed actual active set by one or two labelled groups.
Inactive supports are exactly empty; every regional selected set is the
image of the same active-group embedding, so selected slots are distinct. -/
noncomputable section
namespace Spin.Turbo
open Finset

def sparseActiveEmbedding22 {R q : ℕ} (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S}) : Fin q ↪ Fin (R * 16) :=
  ⟨fun a => (e a).val, fun a b h => e.injective (Subtype.ext h)⟩

def sparseEnumeratedSupport22 {R q : ℕ} (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin q → Finset (Fin 64)) : Fin (R * 16) → Finset (Fin 64) :=
  sparseExtendSupport22 S (fun s => A (e.symm s))

theorem sparseExtendSupport22_mem {R : ℕ} (S : Finset (Fin (R * 16)))
    (a : {i : Fin (R * 16) // i ∈ S} → Finset (Fin 64))
    (i : {i : Fin (R * 16) // i ∈ S}) : sparseExtendSupport22 S a i.val = a i := by
  simp [sparseExtendSupport22, Equiv.piEquivPiSubtypeProd, i.property]

theorem sparseExtendSupport22_not_mem {R : ℕ} (S : Finset (Fin (R * 16)))
    (a : {i : Fin (R * 16) // i ∈ S} → Finset (Fin 64))
    (i : Fin (R * 16)) (hi : i ∉ S) : sparseExtendSupport22 S a i = ∅ := by
  simp [sparseExtendSupport22, Equiv.piEquivPiSubtypeProd, hi]

theorem sparseEnumeratedSupport22_apply {R q : ℕ} (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin q → Finset (Fin 64)) (a : Fin q) :
    sparseEnumeratedSupport22 S e A (e a).val = A a := by
  unfold sparseEnumeratedSupport22
  rw [sparseExtendSupport22_mem, e.symm_apply_apply]
theorem sparseEnumeratedSupport22_region {R q : ℕ} (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin q → Finset (Fin 64)) (region : Fin 64) :
    sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region =
      (Finset.univ.filter (fun a : Fin q => region ∈ A a)).map (sparseActiveEmbedding22 S e) := by
  ext i
  simp only [sparseRegionalSupport22, Finset.mem_filter, Finset.mem_univ, true_and]
  constructor
  · intro hr
    have hi : i ∈ S := by
      by_contra hn
      have he := sparseExtendSupport22_not_mem S (fun s => A (e.symm s)) i hn
      change region ∈ sparseExtendSupport22 S (fun s => A (e.symm s)) i at hr
      rw [he] at hr
      simpa using hr
    let a := e.symm ⟨i, hi⟩
    have he : (e a).val = i := congrArg Subtype.val (e.apply_symm_apply ⟨i, hi⟩)
    refine Finset.mem_map.mpr ⟨a, ?_, he⟩
    apply Finset.mem_filter.mpr
    refine ⟨Finset.mem_univ _, ?_⟩
    have hw := sparseExtendSupport22_mem S (fun s => A (e.symm s)) ⟨i, hi⟩
    change region ∈ sparseExtendSupport22 S (fun s => A (e.symm s)) i at hr
    rw [hw] at hr
    exact hr
  · intro hr
    obtain ⟨a, ha, he⟩ := Finset.mem_map.mp hr
    change (e a).val = i at he
    rw [← he]
    change region ∈ sparseExtendSupport22 S (fun s => A (e.symm s)) (e a).val
    rw [sparseExtendSupport22_mem, e.symm_apply_apply]
    exact (Finset.mem_filter.mp ha).2

theorem sparseEnumeratedSupport22_region_card {R q : ℕ}
    (S : Finset (Fin (R * 16))) (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S})
    (A : Fin q → Finset (Fin 64)) (region : Fin 64) :
    (sparseRegionalSupport22 (sparseEnumeratedSupport22 S e A) region).card =
      ∑ a : Fin q, if region ∈ A a then 1 else 0 := by
  rw [sparseEnumeratedSupport22_region, Finset.card_map, Finset.card_filter]

theorem sparseSupportProfileCount22_reindex {R q : ℕ}
    (points : Fin 8 → Field16) (S : Finset (Fin (R * 16)))
    (e : Fin q ≃ {i : Fin (R * 16) // i ∈ S})
    (F : (Fin (R * 16) → Finset (Fin 64)) → ℝ) :
    (∑ supports : Fin (R * 16) → Finset (Fin 64),
      (∏ i, sparseSupportCount22 points (i ∈ S) (supports i)) * F supports) =
      ∑ A : Fin q → Finset (Fin 64),
        (∏ a : Fin q, rs8ExpectedPacketCount points (A a).card /
          (Nat.choose 64 (A a).card : ℝ)) * F (sparseEnumeratedSupport22 S e A) := by
  letI : Fintype {i : Fin (R * 16) // i ∈ S} :=
    Subtype.fintype (fun i : Fin (R * 16) => i ∈ S)
  rw [sparseSupportProfileCount22_prune]
  let ep := Equiv.arrowCongr e (Equiv.refl (Finset (Fin 64)))
  rw [← Equiv.sum_comp ep]
  apply Finset.sum_congr rfl
  intro A _
  have hw : (∏ i : {i : Fin (R * 16) // i ∈ S},
      rs8ExpectedPacketCount points ((ep A) i).card /
        (Nat.choose 64 ((ep A) i).card : ℝ)) =
      ∏ a : Fin q, rs8ExpectedPacketCount points (A a).card / (Nat.choose 64 (A a).card : ℝ) := by
    rw [← Equiv.prod_comp e]
    simp only [ep, Equiv.arrowCongr_apply, Function.comp_def, Equiv.symm_apply_apply, Equiv.refl_apply]
  rw [hw]
  rfl

#print axioms sparseEnumeratedSupport22_region_card
#print axioms sparseSupportProfileCount22_reindex
end Spin.Turbo


