import SpinCodes.Structured.ConcretePlacementGaps

/-! Exact uniform local coordinates in each collision-free block-placement fiber. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing
attribute [local instance] Classical.propDecidable

def blockBit {R : Nat} (b : Fin R) (c : Fin 128) : Fin (128 * R) :=
  ⟨128 * b.val + c.val, by have := b.isLt; have := c.isLt; omega⟩

@[simp] theorem blockBit_div {R : Nat} (b : Fin R) (c : Fin 128) :
    (blockBit b c).val / 128 = b.val := by
  simp only [blockBit]
  omega

@[simp] theorem blockBit_mod {R : Nat} (b : Fin R) (c : Fin 128) :
    (blockBit b c).val % 128 = c.val := by
  simp only [blockBit]
  omega

def singletonSupport {R a : Nat} (blocks : Fin a ↪ Fin R) (coords : Fin a → Fin 128) :
    Finset (Fin (128 * R)) := univ.image (fun i => blockBit (blocks i) (coords i))

theorem blockBit_index_injective {R a : Nat} (blocks : Fin a ↪ Fin R)
    (coords : Fin a → Fin 128) : Function.Injective (fun i => blockBit (blocks i) (coords i)) := by
  intro i j h
  have hh := congrArg (fun x : Fin (128 * R) => x.val / 128) h
  simp only [blockBit_div] at hh
  exact blocks.injective (Fin.ext hh)

@[simp] theorem singletonSupport_card {R a : Nat} (blocks : Fin a ↪ Fin R)
    (coords : Fin a → Fin 128) : (singletonSupport blocks coords).card = a := by
  rw [singletonSupport, card_image_of_injective _ (blockBit_index_injective blocks coords)]
  exact Fintype.card_fin a

theorem singletonSupport_injective {R a : Nat} (blocks : Fin a ↪ Fin R) :
    Function.Injective (singletonSupport blocks) := by
  intro c d h
  funext i
  have hm : blockBit (blocks i) (c i) ∈ singletonSupport blocks c := mem_image.mpr ⟨i, mem_univ _, rfl⟩
  rw [h] at hm
  obtain ⟨j, _, hj⟩ := mem_image.mp hm
  have hb := congrArg (fun x : Fin (128 * R) => x.val / 128) hj
  simp only [blockBit_div] at hb
  have hji := blocks.injective (Fin.ext hb)
  subst j
  have hc := congrArg (fun x : Fin (128 * R) => x.val % 128) hj
  simp only [blockBit_mod] at hc
  exact Fin.ext hc.symm

def singletonFiber {R a : Nat} (blocks : Fin a ↪ Fin R) : Finset (Finset (Fin (128 * R))) :=
  univ.image (singletonSupport blocks)

@[simp] theorem singletonFiber_card {R a : Nat} (blocks : Fin a ↪ Fin R) :
    (singletonFiber blocks).card = 128 ^ a := by
  rw [singletonFiber, card_image_of_injective _ (singletonSupport_injective blocks)]
  simp

theorem shuffle_expect_singletonFiber {R a : Nat} (S : Finset (Fin (128 * R)))
    (hS : S.card = a) (blocks : Fin a ↪ Fin R) (f : Finset (Fin (128 * R)) → ℝ) :
    (shuffleLaw S).expect (fun T => if T ∈ singletonFiber blocks then f T else 0) =
      (1 / ((128 * R).choose a : ℝ)) * ∑ coords : Fin a → Fin 128, f (singletonSupport blocks coords) := by
  unfold Spin.FinPMF.expect
  simp only [mul_ite, mul_zero]
  rw [← sum_filter]
  have he : univ.filter (fun T => T ∈ singletonFiber blocks) = singletonFiber blocks := by ext T; simp
  rw [he, singletonFiber, sum_image]
  · rw [mul_sum]
    apply sum_congr rfl
    intro coords _
    simp [shuffleLaw_apply, singletonSupport_card, hS]
  · exact fun _ _ _ _ h => singletonSupport_injective blocks h

theorem shuffle_prob_singletonFiber {R a : Nat} (S : Finset (Fin (128 * R)))
    (hS : S.card = a) (blocks : Fin a ↪ Fin R) :
    (shuffleLaw S).prob (fun T => T ∈ singletonFiber blocks) =
      (128 : ℝ)^a / ((128 * R).choose a : ℝ) := by
  have h := shuffle_expect_singletonFiber S hS blocks (fun _ => 1)
  have he : (shuffleLaw S).prob (fun T => T ∈ singletonFiber blocks) =
      (shuffleLaw S).expect (fun T => if T ∈ singletonFiber blocks then 1 else 0) := by
    simp [Spin.FinPMF.prob, Spin.FinPMF.expect, Finset.sum_filter]
  rw [he]
  simpa only [sum_const, card_univ, Fintype.card_fun, Fintype.card_fin,
    nsmul_eq_mul, mul_one, Nat.cast_pow, Nat.cast_ofNat, one_div_mul_eq_div] using h

theorem coordinate_tuple_law (a : Nat) :
    Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128)) =
      Spin.FinPMF.uniform (Fin a → Fin 128) := by
  apply Spin.FinPMF.ext
  funext coords
  simp [Spin.piPMF_apply, Spin.FinPMF.uniform]

theorem shuffle_singletonFiber_independent_coordinates {R a : Nat}
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (blocks : Fin a ↪ Fin R)
    (f : Finset (Fin (128 * R)) → ℝ) :
    (shuffleLaw S).expect (fun T => if T ∈ singletonFiber blocks then f T else 0) =
      (shuffleLaw S).prob (fun T => T ∈ singletonFiber blocks) *
        (Spin.piPMF (fun _ : Fin a => Spin.FinPMF.uniform (Fin 128))).expect
          (fun coords => f (singletonSupport blocks coords)) := by
  rw [shuffle_expect_singletonFiber S hS, shuffle_prob_singletonFiber S hS, coordinate_tuple_law]
  simp only [Spin.FinPMF.expect, Spin.FinPMF.uniform, Fintype.card_fun, Fintype.card_fin,
    Nat.cast_pow, Nat.cast_ofNat, ← mul_sum]
  have hp : (128 : ℝ)^a ≠ 0 := pow_ne_zero _ (by norm_num)
  field_simp

end Spin.Structured.Placement

