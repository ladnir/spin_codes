import SpinCodes.UniformMixer.ScalableFeedbackFourier

/-!
# Exact zero-feedback Bernoulli moment from the production spectrum

The 22-bit character sum evaluates the zero-feedback indicator. The iid
product law then turns each character into a two-factor weight polynomial.
This avoids enumerating the 42-dimensional feedback kernel.
-/

noncomputable section

namespace Spin.UniformMixer

theorem scalableZeroHold_fourier_sum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    (2 ^ 22 : ℝ) * scalableZeroHold β z hβ0 hβ1 =
      ∑ q : State22,
        feedbackWeightMoment β z (coordinateWeight coordMap q) := by
  let P := bernoulliInputLaw β hβ0 hβ1
  change (2 ^ 22 : ℝ) *
      (∑ x : Input64,
        P.p x * (if stepFeedback x = 0 then z ^ stepWeight x 0 else 0)) = _
  rw [Finset.mul_sum]
  calc
    (∑ x : Input64,
      (2 ^ 22 : ℝ) *
        (P.p x * (if stepFeedback x = 0 then z ^ stepWeight x 0 else 0))) =
      ∑ x : Input64,
        P.p x * (z ^ stepWeight x 0 *
          ∑ q : State22,
            (binarySign (statePairing q (stepFeedback x)) : ℝ)) := by
      apply Finset.sum_congr rfl
      intro x _
      rw [← feedback_indicator_character x]
      split_ifs <;> ring
    _ = ∑ q : State22,
          P.expect (fun x => z ^ stepWeight x 0 *
            (binarySign (statePairing q (stepFeedback x)) : ℝ)) := by
      simp only [FinPMF.expect, Finset.mul_sum]
      rw [Finset.sum_comm]
    _ = _ := by
      apply Finset.sum_congr rfl
      intro q _
      exact feedbackCharacterWeighted_moment β z hβ0 hβ1 q

theorem scalableZeroHold_exact_spectrum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableZeroHold β z hβ0 hβ1 =
      (∑ w ∈ spectrumWeights,
        (productionWeightCount w : ℝ) * feedbackWeightMoment β z w) /
        (2 ^ 22 : ℝ) := by
  have h := scalableZeroHold_fourier_sum β z hβ0 hβ1
  rw [productionSpectrum_sum (feedbackWeightMoment β z)] at h
  apply (eq_div_iff (by norm_num : (2 ^ 22 : ℝ) ≠ 0)).mpr
  simpa only [mul_comm] using h

theorem scalableZeroHold_at_one (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableZeroHold β 1 hβ0 hβ1 =
      (bernoulliInputLaw β hβ0 hβ1).prob
        (fun x => stepFeedback x = 0) := by
  unfold scalableZeroHold
  rw [FinPMF.prob_eq_expect_indicator]
  apply congrArg ((bernoulliInputLaw β hβ0 hβ1).expect)
  funext x
  simp only [one_pow]

theorem scalableFeedbackActive_complement (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableFeedbackActive β hβ0 hβ1 =
      1 - scalableZeroHold β 1 hβ0 hβ1 := by
  change (bernoulliInputLaw β hβ0 hβ1).prob
    (fun x => ¬ stepFeedback x = 0) = _
  rw [FinPMF.prob_compl]
  rw [scalableZeroHold_at_one]

theorem feedbackWeightMoment_at_one (β : ℝ) (w : Nat) :
    feedbackWeightMoment β 1 w = (1 - 2 * β) ^ w := by
  unfold feedbackWeightMoment
  ring

theorem scalableFeedbackActive_exact_spectrum (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableFeedbackActive β hβ0 hβ1 =
      1 - (∑ w ∈ spectrumWeights,
        (productionWeightCount w : ℝ) * (1 - 2 * β) ^ w) /
        (2 ^ 22 : ℝ) := by
  rw [scalableFeedbackActive_complement,
    scalableZeroHold_exact_spectrum]
  apply congrArg (fun v : ℝ => 1 - v / (2 ^ 22 : ℝ))
  apply Finset.sum_congr rfl
  intro w _
  rw [feedbackWeightMoment_at_one]

end Spin.UniformMixer
