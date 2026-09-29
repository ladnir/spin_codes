import SpinCodes.Structured.ConcretePlacementEncoderMoment

/-! Unconditional actual shuffled encoder moment versus the explicit finite coarse-product sum. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder FiniteKernel
attribute [local instance] Classical.propDecidable

def placementOmissionError (a R H : Nat) : ℝ :=
  256 * (H : ℝ) * a / (128 * R) +
    (256 * ((H : ℝ) + 1) + 1) * (a : ℝ)^2 / (128 * (R : ℝ) - 1)

def finiteGoodMoment (θ : ℝ) (R a H : Nat) (q : State) : ℝ :=
  ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
    ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then coarseRegionMoment θ B q else 0

theorem good_placement_mass_le_one {R a : Nat} (S : Finset (Fin (128 * R)))
    (hS : S.card = a) (H : Nat) :
    ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
      (∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then (1 : ℝ) else 0) ≤ 1 := by
  have he := shuffle_good_blocks_sum S hS H (fun _ => 1)
  simp only [Spin.FinPMF.expect_const] at he
  rw [← he]
  calc
    _ ≤ (shuffleLaw S).expect (fun _ => 1) := by
      apply Spin.FinPMF.expect_mono
      intro T
      split_ifs <;> norm_num
    _ = 1 := Spin.FinPMF.expect_const _ _

theorem weighted_restricted_sum_error {Ω : Type*} [Fintype Ω]
    (E : Ω → Prop) [DecidablePred E] (f g : Ω → ℝ) {w ε : ℝ} (hw : 0 ≤ w) (hε : 0 ≤ ε)
    (hm : w * (∑ x, if E x then (1 : ℝ) else 0) ≤ 1)
    (h : ∀ x, E x → |f x - g x| ≤ ε) :
    |w * (∑ x, if E x then f x else 0) - w * ∑ x, if E x then g x else 0| ≤ ε := by
  rw [← mul_sub, abs_mul, abs_of_nonneg hw, ← sum_sub_distrib]
  calc _ ≤ w * ∑ x, |(if E x then f x else 0) - if E x then g x else 0| :=
      mul_le_mul_of_nonneg_left (Finset.abs_sum_le_sum_abs _ _) hw
    _ ≤ w * ∑ x, (if E x then (1 : ℝ) else 0) * ε := by
      apply mul_le_mul_of_nonneg_left _ hw
      apply sum_le_sum
      intro x _
      by_cases hx : E x
      · simpa only [hx, ite_true, one_mul] using h x hx
      · simp only [hx, ite_false, sub_self, abs_zero, zero_mul, le_refl]
    _ = (w * ∑ x, if E x then (1 : ℝ) else 0) * ε := by rw [← sum_mul]; ring
    _ ≤ ε := by simpa only [one_mul] using mul_le_mul_of_nonneg_right hm hε

theorem shuffled_regionMoment_coarse_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (H : Nat) (hHR : H ≤ R) (q : State) :
    |(shuffleLaw S).expect (fun T => regionMoment (Real.exp (-(θ / (128 * R)))) T q) -
      finiteGoodMoment θ R a H q| ≤ placementOmissionError a R H + regionKernelError a R H θ := by
  let z := Real.exp (-(θ / (128 * (R : ℝ))))
  let F : BlockSubset R a → ℝ := fun B =>
    (Spin.piPMF (fun _ : Fin a => coordinateLaw)).expect
      (fun coords => regionMoment z (singletonSupport (orderedBlocks B) coords) q)
  let w : ℝ := (128 : ℝ)^a / ((128 * R).choose a : ℝ)
  let M : ℝ := w * ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then F B else 0
  have hz0 : 0 ≤ z := (Real.exp_pos _).le
  have hz1 : z ≤ 1 := Real.exp_le_one_iff.mpr (neg_nonpos.mpr (by positivity))
  have hfirst : |(shuffleLaw S).expect (fun T => regionMoment z T q) - M| ≤ placementOmissionError a R H := by
    exact shuffle_good_blocks_sum_error hR S hS H (fun T => regionMoment z T q)
      (fun T => regionMoment_unit hz0 hz1 T q)
  have hsecond : |M - finiteGoodMoment θ R a H q| ≤ regionKernelError a R H θ := by
    have hh := weighted_restricted_sum_error
      (fun B : BlockSubset R a => ¬ badBlockIndices (orderedBlocks B) H)
      F (fun B => coarseRegionMoment θ B q) (w := w) (ε := regionKernelError a R H θ)
      (by positivity) (by unfold regionKernelError; positivity)
      (good_placement_mass_le_one S hS H)
      (fun B hg => regionMoment_good_coordinates_error hR hθ B H hHR hg q)
    exact hh
  exact (abs_sub_le _ M _).trans (add_le_add hfirst hsecond)

end Spin.Structured.Placement



