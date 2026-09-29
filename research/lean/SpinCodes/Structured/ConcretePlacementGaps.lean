import SpinCodes.Structured.ConcretePlacementBounds

/-! Boundary and block-gap errors for uniformly shuffled fixed-size supports. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

def boundarySites (n D : Nat) : Finset (Fin n) :=
  univ.filter (fun x => x.val < D ∨ n ≤ x.val + D)

theorem boundarySites_card_bound (n D : Nat) : (boundarySites n D).card ≤ 2 * D := by
  let A := univ.filter (fun x : Fin n => x.val < D)
  let B := univ.filter (fun x : Fin n => n ≤ x.val + D)
  have hA : A.card ≤ D := by
    have hs : A.image Fin.val ⊆ range D := by
      intro y hy
      obtain ⟨x, hx, rfl⟩ := mem_image.mp hy
      exact mem_range.mpr (mem_filter.mp hx).2
    calc A.card = (A.image Fin.val).card := (card_image_of_injective _ Fin.val_injective).symm
      _ ≤ (range D).card := card_le_card hs
      _ = D := card_range D
  have hB : B.card ≤ D := by
    have hs : B.image Fin.val ⊆ Ico (n - D) n := by
      intro y hy
      obtain ⟨x, hx, rfl⟩ := mem_image.mp hy
      have h := (mem_filter.mp hx).2
      simp only [mem_Ico]
      exact ⟨by omega, x.isLt⟩
    calc B.card = (B.image Fin.val).card := (card_image_of_injective _ Fin.val_injective).symm
      _ ≤ (Ico (n - D) n).card := card_le_card hs
      _ ≤ D := by rw [Nat.card_Ico]; omega
  have he : boundarySites n D = A ∪ B := by ext x; simp [boundarySites, A, B]
  rw [he]
  exact (card_union_le _ _).trans (by omega)

theorem shuffle_prob_boundary {n : Nat} (S : Finset (Fin n)) (D : Nat) :
    (shuffleLaw S).prob (fun T => ∃ x ∈ T, x.val < D ∨ n ≤ x.val + D) ≤
      2 * (D : ℝ) * S.card / n := by
  have he : (shuffleLaw S).prob (fun T => ∃ x ∈ T, x.val < D ∨ n ≤ x.val + D) =
      (shuffleLaw S).prob (fun T => ∃ x ∈ boundarySites n D, x ∈ T) := by
    apply Spin.FinPMF.prob_congr
    intro T
    simp only [boundarySites, mem_filter, mem_univ, true_and]
    exact ⟨fun ⟨x, hx, hb⟩ => ⟨x, hb, hx⟩, fun ⟨x, hb, hx⟩ => ⟨x, hx, hb⟩⟩
  rw [he]
  have hc : ((boundarySites n D).card : ℝ) ≤ 2 * (D : ℝ) := by
    exact_mod_cast boundarySites_card_bound n D
  exact (shuffle_prob_hits S _).trans (div_le_div_of_nonneg_right
    (mul_le_mul_of_nonneg_right hc (Nat.cast_nonneg _)) (Nat.cast_nonneg _))

def badBlockPlacement {R : Nat} (H : Nat) (T : Finset (Fin (128 * R))) : Prop :=
  (∃ x ∈ T, x.val / 128 < H ∨ R ≤ x.val / 128 + H) ∨
  (∃ x ∈ T, ∃ y ∈ T, x ≠ y ∧ Nat.dist (x.val / 128) (y.val / 128) ≤ H)

theorem bit_dist_of_block_dist (x y H : Nat)
    (h : Nat.dist (x / 128) (y / 128) ≤ H) : Nat.dist x y ≤ 128 * (H + 1) := by
  have hx := Nat.mod_lt x (by omega : 0 < 128)
  have hy := Nat.mod_lt y (by omega : 0 < 128)
  have hxx := Nat.div_add_mod x 128
  have hyy := Nat.div_add_mod y 128
  unfold Nat.dist at *
  omega

theorem shuffle_prob_badBlockPlacement {R : Nat} (hR : 1 ≤ R)
    (S : Finset (Fin (128 * R))) (H : Nat) :
    (shuffleLaw S).prob (badBlockPlacement H) ≤
      256 * (H : ℝ) * S.card / (128 * R) +
      (256 * ((H : ℝ) + 1) + 1) * (S.card : ℝ)^2 / (128 * (R : ℝ) - 1) := by
  have hb := shuffle_prob_boundary S (128 * H)
  have hc := shuffle_prob_close (by omega : 2 ≤ 128 * R) S (128 * (H + 1))
  have hmono : (shuffleLaw S).prob (badBlockPlacement H) ≤ (shuffleLaw S).prob
      (fun T => (∃ x ∈ T, x.val < 128 * H ∨ 128 * R ≤ x.val + 128 * H) ∨
        (∃ x ∈ T, ∃ y ∈ T, x ≠ y ∧ Nat.dist x.val y.val ≤ 128 * (H + 1))) := by
    apply Spin.FinPMF.prob_mono
    intro T h
    rcases h with ⟨x, hx, hl | hr⟩ | ⟨x, hx, y, hy, hxy, hd⟩
    · exact Or.inl ⟨x, hx, Or.inl (by omega)⟩
    · exact Or.inl ⟨x, hx, Or.inr (by omega)⟩
    · exact Or.inr ⟨x, hx, y, hy, hxy, bit_dist_of_block_dist _ _ _ hd⟩
  calc _ ≤ _ := hmono
    _ ≤ _ := Spin.FinPMF.prob_or_le _ _ _
    _ ≤ _ := add_le_add hb hc
    _ = _ := by push_cast; ring

end Spin.Structured.Placement
