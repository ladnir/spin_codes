import SpinCodes.Structured.ConcretePlacementFibers
import Mathlib.Data.Finset.Sort

/-! The actual distinct-block placement law as a uniform ordered-block finite sum. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

def blockOf {R : Nat} (x : Fin (128 * R)) : Fin R :=
  ⟨x.val / 128, by have := x.isLt; omega⟩

@[simp] theorem blockOf_blockBit {R : Nat} (b : Fin R) (c : Fin 128) :
    blockOf (blockBit b c) = b := Fin.ext (blockBit_div b c)

def occupiedBlocks {R : Nat} (T : Finset (Fin (128 * R))) : Finset (Fin R) :=
  T.image blockOf

abbrev BlockSubset (R a : Nat) := {B : Finset (Fin R) // B.card = a}

def orderedBlocks {R a : Nat} (B : BlockSubset R a) : Fin a ↪ Fin R :=
  (B.val.orderEmbOfFin B.property).toEmbedding

theorem orderedBlocks_strictMono {R a : Nat} (B : BlockSubset R a) :
    StrictMono (orderedBlocks B) := (B.val.orderEmbOfFin B.property).strictMono

@[simp] theorem orderedBlocks_image {R a : Nat} (B : BlockSubset R a) :
    univ.image (orderedBlocks B) = B.val := Finset.image_orderEmbOfFin_univ _ _

@[simp] theorem occupiedBlocks_singletonSupport {R a : Nat} (blocks : Fin a ↪ Fin R)
    (coords : Fin a → Fin 128) :
    occupiedBlocks (singletonSupport blocks coords) = univ.image blocks := by
  simp only [occupiedBlocks, singletonSupport, Finset.image_image, Function.comp_def, blockOf_blockBit]

theorem singletonFiber_iff_occupied {R a : Nat} (B : BlockSubset R a)
    (T : Finset (Fin (128 * R))) :
    T ∈ singletonFiber (orderedBlocks B) ↔ T.card = a ∧ occupiedBlocks T = B.val := by
  constructor
  · intro h
    obtain ⟨coords, _, rfl⟩ := mem_image.mp h
    exact ⟨singletonSupport_card _ _, by simp⟩
  · rintro ⟨hc, hb⟩
    rw [mem_singletonFiber_iff]
    refine ⟨hc, fun i => ?_⟩
    have hi : orderedBlocks B i ∈ B.val := B.val.orderEmbOfFin_mem B.property i
    rw [← hb] at hi
    obtain ⟨x, hx, hxi⟩ := mem_image.mp hi
    exact ⟨x, hx, congrArg Fin.val hxi⟩

theorem good_occupiedBlocks_card {R : Nat} (T : Finset (Fin (128 * R))) (H : Nat)
    (h : ¬ badBlockPlacement H T) : (occupiedBlocks T).card = T.card := by
  apply Finset.card_image_iff.mpr
  intro x hx y hy he
  by_contra hxy
  apply h
  right
  refine ⟨x, hx, y, hy, hxy, ?_⟩
  have hv := congrArg Fin.val he
  change x.val / 128 = y.val / 128 at hv
  rw [hv, Nat.dist_self]
  exact Nat.zero_le _

theorem sum_blockSubset_indicator {R a : Nat} (U : Finset (Fin R)) :
    (∑ B : BlockSubset R a, if U = B.val then (1 : ℝ) else 0) =
      if U.card = a then 1 else 0 := by
  by_cases h : U.card = a
  · rw [if_pos h]
    let B : BlockSubset R a := ⟨U, h⟩
    rw [Finset.sum_eq_single B]
    · simp [B]
    · intro C _ hCB
      rw [if_neg]
      intro he
      apply hCB
      exact Subtype.ext he.symm
    · simp
  · rw [if_neg h]
    apply sum_eq_zero
    intro B _
    rw [if_neg]
    intro he
    exact h (he ▸ B.property)

theorem shuffle_distinct_blocks_sum {R a : Nat} (S : Finset (Fin (128 * R)))
    (hS : S.card = a) (f : Finset (Fin (128 * R)) → ℝ) :
    (shuffleLaw S).expect (fun T => if (occupiedBlocks T).card = a then f T else 0) =
      ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
        ∑ B : BlockSubset R a,
          (Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128))).expect
            (fun coords => f (singletonSupport (orderedBlocks B) coords)) := by
  have he : (shuffleLaw S).expect (fun T => if (occupiedBlocks T).card = a then f T else 0) =
      ∑ B : BlockSubset R a, (shuffleLaw S).expect
        (fun T => if T ∈ singletonFiber (orderedBlocks B) then f T else 0) := by
    unfold Spin.FinPMF.expect
    rw [sum_comm]
    apply sum_congr rfl
    intro T _
    by_cases ht : T.card = a
    · simp only [singletonFiber_iff_occupied, ht, true_and]
      rw [← mul_sum]
      have hh := sum_blockSubset_indicator (a := a) (occupiedBlocks T)
      have hf : (∑ B : BlockSubset R a, if occupiedBlocks T = B.val then f T else 0) =
          (if (occupiedBlocks T).card = a then f T else 0) := by
        calc _ = (∑ B : BlockSubset R a, if occupiedBlocks T = B.val then (1 : ℝ) else 0) * f T := by
              rw [sum_mul]; apply sum_congr rfl; intro B _; split_ifs <;> simp
          _ = _ := by rw [hh]; split_ifs <;> simp
      rw [hf]
    · simp [shuffleLaw_apply, hS, ht]
  rw [he, mul_sum]
  apply sum_congr rfl
  intro B _
  rw [shuffle_singletonFiber_independent_coordinates S hS, shuffle_prob_singletonFiber S hS]

end Spin.Structured.Placement
