import SpinCodes.Structured.ConcretePlacementCoordinates

/-! Intrinsic characterization of singleton fibers and coordinate-independent good gaps. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

theorem mem_singletonFiber_iff {R a : Nat} (blocks : Fin a ↪ Fin R)
    (T : Finset (Fin (128 * R))) : T ∈ singletonFiber blocks ↔
      T.card = a ∧ ∀ i, ∃ x ∈ T, x.val / 128 = (blocks i).val := by
  constructor
  · intro h
    obtain ⟨coords, _, rfl⟩ := mem_image.mp h
    refine ⟨singletonSupport_card _ _, fun i => ?_⟩
    exact ⟨blockBit (blocks i) (coords i), mem_image.mpr ⟨i, mem_univ _, rfl⟩, blockBit_div _ _⟩
  · rintro ⟨hc, hw⟩
    choose x hx hd using hw
    let coords : Fin a → Fin 128 := fun i => ⟨(x i).val % 128, Nat.mod_lt _ (by omega)⟩
    have he (i : Fin a) : blockBit (blocks i) (coords i) = x i := by
      apply Fin.ext
      dsimp [blockBit, coords]
      rw [← hd i]
      omega
    have hs : singletonSupport blocks coords ⊆ T := by
      intro y hy
      obtain ⟨i, _, rfl⟩ := mem_image.mp hy
      rw [he]
      exact hx i
    have ht : singletonSupport blocks coords = T := eq_of_subset_of_card_le hs (by simp [hc])
    exact mem_image.mpr ⟨coords, mem_univ _, ht⟩

def badBlockIndices {R a : Nat} (blocks : Fin a ↪ Fin R) (H : Nat) : Prop :=
  (∃ i, (blocks i).val < H ∨ R ≤ (blocks i).val + H) ∨
  (∃ i j, i ≠ j ∧ Nat.dist (blocks i).val (blocks j).val ≤ H)

theorem badBlockPlacement_singletonSupport {R a : Nat} (blocks : Fin a ↪ Fin R)
    (coords : Fin a → Fin 128) (H : Nat) :
    badBlockPlacement H (singletonSupport blocks coords) ↔ badBlockIndices blocks H := by
  unfold badBlockPlacement badBlockIndices
  constructor
  · rintro (⟨x, hx, hd⟩ | ⟨x, hx, y, hy, hxy, hd⟩)
    · obtain ⟨i, _, rfl⟩ := mem_image.mp hx
      exact Or.inl ⟨i, by simpa only [blockBit_div] using hd⟩
    · obtain ⟨i, _, rfl⟩ := mem_image.mp hx
      obtain ⟨j, _, rfl⟩ := mem_image.mp hy
      exact Or.inr ⟨i, j, fun h => hxy (by rw [h]), by simpa only [blockBit_div] using hd⟩
  · rintro (⟨i, hi⟩ | ⟨i, j, hij, hd⟩)
    · exact Or.inl ⟨blockBit (blocks i) (coords i), mem_image.mpr ⟨i, mem_univ _, rfl⟩,
        by simpa only [blockBit_div] using hi⟩
    · exact Or.inr ⟨blockBit (blocks i) (coords i), mem_image.mpr ⟨i, mem_univ _, rfl⟩,
        blockBit (blocks j) (coords j), mem_image.mpr ⟨j, mem_univ _, rfl⟩,
        fun h => hij (blockBit_index_injective blocks coords h),
        by simpa only [blockBit_div] using hd⟩

theorem shuffle_good_singletonFiber_coordinates {R a : Nat}
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (blocks : Fin a ↪ Fin R)
    (H : Nat) (hg : ¬ badBlockIndices blocks H) (f : Finset (Fin (128 * R)) → ℝ) :
    (shuffleLaw S).expect (fun T => if T ∈ singletonFiber blocks ∧ ¬ badBlockPlacement H T then f T else 0) =
      (shuffleLaw S).prob (fun T => T ∈ singletonFiber blocks) *
        (Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128))).expect
          (fun coords => f (singletonSupport blocks coords)) := by
  have he : (fun T => if T ∈ singletonFiber blocks ∧ ¬ badBlockPlacement H T then f T else 0) =
      (fun T => if T ∈ singletonFiber blocks then f T else 0) := by
    funext T
    by_cases ht : T ∈ singletonFiber blocks
    · obtain ⟨coords, _, rfl⟩ := mem_image.mp ht
      have hm : singletonSupport blocks coords ∈ singletonFiber blocks :=
        mem_image.mpr ⟨coords, mem_univ _, rfl⟩
      simp only [hm, badBlockPlacement_singletonSupport, hg, not_false_eq_true, and_self, ite_true]
    · simp only [ht, false_and, ite_false]
  rw [he]
  exact shuffle_singletonFiber_independent_coordinates S hS blocks f

end Spin.Structured.Placement
