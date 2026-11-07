import SpinCodes.UniformMixer.ScalableFirstMomentBound

/-!
# Weighted rows of the iid step kernel

These formulas expose the two events relevant to the paper's
three-coordinate envelope: staying at zero and returning from a live
state. They retain the actual 22-bit state and production maps.
-/

noncomputable section

namespace Spin.UniformMixer

def scalableLiveMoment (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (q : State22) : ℝ :=
  (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => z ^ stepWeight x q)

def scalableZeroHold (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) : ℝ :=
  (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => if stepFeedback x = 0 then z ^ stepWeight x 0 else 0)

def scalableFeedbackActive (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) : ℝ :=
  (bernoulliInputLaw β hβ0 hβ1).prob
    (fun x => stepFeedback x ≠ 0)

theorem bernoulliStepKernel_zero_entry (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (r : State22) :
    bernoulliStepKernel β z hβ0 hβ1 0 r =
      (bernoulliInputLaw β hβ0 hβ1).expect
        (fun x => if r = stepFeedback x then z ^ stepWeight x 0 else 0) := by
  change (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => stepKernel z x 0 r) = _
  congr 1
  funext x
  simp [stepKernel, update_zero]

theorem bernoulliStepKernel_zero_zero (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    bernoulliStepKernel β z hβ0 hβ1 0 0 =
      scalableZeroHold β z hβ0 hβ1 := by
  rw [bernoulliStepKernel_zero_entry]
  change (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => if (0 : State22) = stepFeedback x then
      z ^ stepWeight x 0 else 0) =
    (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => if stepFeedback x = 0 then z ^ stepWeight x 0 else 0)
  have hfun :
      (fun x : Input64 => if (0 : State22) = stepFeedback x then
        z ^ stepWeight x 0 else 0) =
      (fun x : Input64 => if stepFeedback x = 0 then
        z ^ stepWeight x 0 else 0) := by
    funext x
    by_cases h : stepFeedback x = 0 <;> simp only [h, eq_comm, ite_true, ite_false]
  exact congrArg (bernoulliInputLaw β hβ0 hβ1).expect hfun

theorem bernoulliStepKernel_live_entry (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (q r : State22) (hq : q ≠ 0) :
    bernoulliStepKernel β z hβ0 hβ1 q r =
      (bernoulliInputLaw β hβ0 hβ1).expect
        (fun x => if r = stepFeedback x then 0 else
          z ^ stepWeight x q / (liveCount State22 : ℝ)) := by
  change (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => stepKernel z x q r) = _
  congr 1
  funext x
  change z ^ stepWeight x q * update (stepFeedback x) q r = _
  rw [update_live _ q r hq]
  split_ifs <;> ring

theorem bernoulliStepKernel_live_atom_le (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z)
    (q r : State22) (hq : q ≠ 0) :
    bernoulliStepKernel β z hβ0 hβ1 q r ≤
      scalableLiveMoment β z hβ0 hβ1 q / (liveCount State22 : ℝ) := by
  rw [bernoulliStepKernel_live_entry β z hβ0 hβ1 q r hq]
  calc
    _ ≤ (bernoulliInputLaw β hβ0 hβ1).expect
        (fun x => z ^ stepWeight x q / (liveCount State22 : ℝ)) := by
      apply FinPMF.expect_mono
      intro x
      split_ifs
      · exact div_nonneg (pow_nonneg hz _) (Nat.cast_nonneg _)
      · exact le_refl _
    _ = _ := by
      unfold scalableLiveMoment FinPMF.expect
      simp only [div_eq_mul_inv, Finset.sum_mul, mul_assoc]

theorem bernoulliStepKernel_live_return_le (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (q : State22) (hq : q ≠ 0) :
    bernoulliStepKernel β z hβ0 hβ1 q 0 ≤
      scalableFeedbackActive β hβ0 hβ1 / (liveCount State22 : ℝ) := by
  rw [bernoulliStepKernel_live_entry β z hβ0 hβ1 q 0 hq]
  calc
    _ ≤ (bernoulliInputLaw β hβ0 hβ1).expect
        (fun x => if stepFeedback x = 0 then 0 else
          1 / (liveCount State22 : ℝ)) := by
      apply FinPMF.expect_mono
      intro x
      by_cases hx : stepFeedback x = 0
      · simp only [hx, ite_true]
        exact le_refl _
      · simp only [hx, ite_false, eq_comm]
        exact div_le_div_of_nonneg_right (pow_le_one₀ hz0 hz1)
          (Nat.cast_nonneg _)
    _ = _ := by
      unfold scalableFeedbackActive
      rw [FinPMF.prob_eq_expect_indicator]
      unfold FinPMF.expect
      simp only [div_eq_mul_inv, Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro x _
      by_cases hx : stepFeedback x = 0
      · simp only [hx, ne_eq, not_true_eq_false, ite_false, ite_true,
          zero_mul, mul_zero]
      · simp only [hx, ne_eq, not_false_eq_true, ite_true, ite_false,
          one_mul, mul_one]

theorem bernoulliStepKernel_row_sum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (q : State22) :
    ∑ r : State22, bernoulliStepKernel β z hβ0 hβ1 q r =
      scalableLiveMoment β z hβ0 hβ1 q := by
  calc
    _ = (bernoulliInputLaw β hβ0 hβ1).expect
        (fun x => ∑ r : State22, stepKernel z x q r) := by
      rw [FinPMF.expect_sum]
      rfl
    _ = _ := by
      unfold scalableLiveMoment
      have hfun :
          (fun x : Input64 => ∑ r : State22, stepKernel z x q r) =
          (fun x : Input64 => z ^ stepWeight x q) := by
        funext x
        simp only [stepKernel, ← Finset.mul_sum, update_sum, mul_one]
      exact congrArg (bernoulliInputLaw β hβ0 hβ1).expect hfun

theorem bernoulliStepKernel_zero_live_mass (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    ∑ r ∈ (Finset.univ : Finset State22).erase 0,
      bernoulliStepKernel β z hβ0 hβ1 0 r =
      scalableLiveMoment β z hβ0 hβ1 0 -
        scalableZeroHold β z hβ0 hβ1 := by
  have h := Finset.add_sum_erase (Finset.univ : Finset State22)
    (fun r => bernoulliStepKernel β z hβ0 hβ1 0 r)
    (Finset.mem_univ (0 : State22))
  rw [bernoulliStepKernel_row_sum, bernoulliStepKernel_zero_zero] at h
  linarith

end Spin.UniformMixer
