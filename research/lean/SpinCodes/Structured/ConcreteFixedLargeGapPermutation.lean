import SpinCodes.Structured.ConcretePlacementSimplexGaps
import Mathlib.Order.Interval.Finset.Fin

/-! The exact finite empty-gap law is invariant under every gap permutation. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset
attribute [local instance] Classical.propDecidable

def gapPrefix {Q : ℕ} (g : Fin (Q+1) → ℕ) (i : Fin (Q+1)) : ℕ := ∑ j ∈ Iic i, g j

theorem gapPrefix_zero {Q : ℕ} (g : Fin (Q+1) → ℕ) : gapPrefix g 0 = g 0 := by
  have h : Iic (0:Fin (Q+1)) = {0} := by
    ext j
    simp only [mem_Iic, mem_singleton]
    constructor
    · intro hj; apply Fin.ext; exact Nat.eq_zero_of_le_zero hj
    · rintro rfl; exact le_rfl
  simp only [gapPrefix, h, sum_singleton]

theorem gapPrefix_step {Q : ℕ} (g : Fin (Q+1) → ℕ) (i : Fin Q) :
    gapPrefix g i.succ = gapPrefix g i.castSucc + g i.succ := by
  have he : Iic i.succ = insert i.succ (Iic i.castSucc) := by
    ext j
    simp only [mem_Iic, mem_insert]
    constructor
    · intro h; by_cases he : j = i.succ
      · exact Or.inl he
      · right; have := j.isLt; change j.val ≤ i.val; change j.val ≤ i.val+1 at h
        have hn : j.val ≠ i.val+1 := fun hh => he (Fin.ext hh)
        omega
    · rintro (rfl|h)
      · exact le_rfl
      · exact h.trans (Fin.castSucc_lt_succ.le)
  unfold gapPrefix
  rw [he, sum_insert (by simp), add_comm]

theorem gapPrefix_last {Q : ℕ} (g : Fin (Q+1) → ℕ) : gapPrefix g (Fin.last Q) = ∑ j, g j := by
  have h : Iic (Fin.last Q) = univ := by ext j; simp only [mem_Iic, mem_univ, iff_true]; exact Fin.le_last j
  simp only [gapPrefix, h]

theorem gapPrefix_emptyGaps {R Q : ℕ} (B : BlockSubset R Q) (i : Fin (Q+1)) :
    gapRight B i = gapPrefix (emptyGaps B) i + i.val := by
  refine Fin.induction ?_ (fun j ih => ?_) i
  · rw [gapPrefix_zero]
    simp [emptyGaps, gapLeft]
  · have hg := Nat.sub_add_cancel (gapLeft_le_right B j.succ)
    change emptyGaps B j.succ + gapLeft B j.succ = gapRight B j.succ at hg
    simp only [gapLeft, Fin.cons_succ] at hg
    rw [gapRight, Fin.snoc_castSucc] at ih
    rw [gapPrefix_step]
    simp only [Fin.val_succ, Fin.val_castSucc] at *
    omega

theorem emptyGaps_injective {R Q : ℕ} : Function.Injective (emptyGaps (R := R) (a := Q)) := by
  intro B C h
  have ho : orderedBlocks B = orderedBlocks C := by
    ext i
    have hB := gapPrefix_emptyGaps B i.castSucc
    have hC := gapPrefix_emptyGaps C i.castSucc
    simp only [gapRight, Fin.snoc_castSucc] at hB hC
    rw [h] at hB
    omega
  apply Subtype.ext
  rw [← orderedBlocks_image B, ← orderedBlocks_image C, ho]

def gapSites {R Q : ℕ} (g : Fin (Q+1) → ℕ) (hg : (∑ i, g i)+Q=R) (i : Fin Q) : Fin R :=
  ⟨gapPrefix g i.castSucc + i.val, by
    have hp : gapPrefix g i.castSucc ≤ ∑ j, g j :=
      sum_le_sum_of_subset_of_nonneg (subset_univ _) (fun _ _ _ => Nat.zero_le _)
    have hi := i.isLt
    omega⟩

theorem gapSites_strictMono {R Q : ℕ} (g : Fin (Q+1) → ℕ) (hg : (∑ i, g i)+Q=R) :
    StrictMono (gapSites g hg) := by
  intro i j hij
  have hp : gapPrefix g i.castSucc ≤ gapPrefix g j.castSucc :=
    sum_le_sum_of_subset_of_nonneg (Iic_subset_Iic.mpr (Fin.castSucc_le_castSucc_iff.mpr hij.le))
      (fun _ _ _ => Nat.zero_le _)
  change gapPrefix g i.castSucc+i.val < gapPrefix g j.castSucc+j.val
  change i.val < j.val at hij
  omega

def gapBlocks {R Q : ℕ} (g : Fin (Q+1) → ℕ) (hg : (∑ i, g i)+Q=R) : BlockSubset R Q :=
  ⟨univ.image (gapSites g hg), by rw [card_image_of_injective _ (gapSites_strictMono g hg).injective]; simp⟩

theorem orderedBlocks_gapBlocks {R Q : ℕ} (g : Fin (Q+1) → ℕ) (hg : (∑ i, g i)+Q=R) :
    orderedBlocks (gapBlocks g hg) = gapSites g hg := by
  symm
  exact Finset.orderEmbOfFin_unique (gapBlocks g hg).property
    (fun i => mem_image.mpr ⟨i, mem_univ _, rfl⟩) (gapSites_strictMono g hg)

theorem emptyGaps_gapBlocks {R Q : ℕ} (g : Fin (Q+1) → ℕ) (hg : (∑ i, g i)+Q=R) :
    emptyGaps (gapBlocks g hg) = g := by
  have hp (i : Fin (Q+1)) : gapPrefix (emptyGaps (gapBlocks g hg)) i = gapPrefix g i := by
    refine Fin.lastCases ?_ (fun j => ?_) i
    · have he := emptyGaps_sum (gapBlocks g hg)
      rw [gapPrefix_last, gapPrefix_last]
      omega
    · have he := gapPrefix_emptyGaps (gapBlocks g hg) j.castSucc
      rw [gapRight, Fin.snoc_castSucc, orderedBlocks_gapBlocks] at he
      change gapPrefix g j.castSucc+j.val = gapPrefix (emptyGaps (gapBlocks g hg)) j.castSucc+j.val at he
      omega
  funext i
  refine Fin.cases ?_ (fun j => ?_) i
  · simpa only [gapPrefix_zero] using hp 0
  · have h1 := hp j.succ
    have h0 := hp j.castSucc
    rw [gapPrefix_step, gapPrefix_step] at h1
    omega

#print axioms emptyGaps_injective
#print axioms emptyGaps_gapBlocks
end Spin.Structured.Placement
