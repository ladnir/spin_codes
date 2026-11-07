import SpinCodes.UniformMixer.ScalableDenseTransferExactDefs
import SpinCodes.UniformMixer.ScalableDenseTransferEval
import SpinCodes.Structured.DenseOccupationFixedVertex

/-! Soundness of exact rational transfer checks against the production matrix. -/

noncomputable section

namespace Spin.UniformMixer

private theorem ratG0_cast (b z : ℚ) :
    ((transferRatG0 b z : ℚ) : ℝ) =
      1 - (b : ℝ) + (b : ℝ) * (z : ℝ) := by
  simp [transferRatG0]

private theorem ratG1_cast (b z : ℚ) :
    ((transferRatG1 b z : ℚ) : ℝ) =
      (b : ℝ) + (1 - (b : ℝ)) * (z : ℝ) := by
  simp [transferRatG1]

private theorem ratH1_cast (b z : ℚ) :
    ((transferRatH1 b z : ℚ) : ℝ) =
      1 - (b : ℝ) - (b : ℝ) * (z : ℝ) := by
  simp [transferRatH1]

theorem ratWeight_cast (b z : ℚ) (w : Nat) :
    ((transferRatWeight b z w : ℚ) : ℝ) =
      productionWeightMoment (b : ℝ) (z : ℝ) w := by
  simp [transferRatWeight, productionWeightMoment, ratG0_cast, ratG1_cast]

theorem ratFeedback_cast (b z : ℚ) (w : Nat) :
    ((transferRatFeedback b z w : ℚ) : ℝ) =
      feedbackWeightMoment (b : ℝ) (z : ℝ) w := by
  simp [transferRatFeedback, feedbackWeightMoment, ratG0_cast, ratH1_cast]

theorem ratNeutral_cast (b : ℚ) (w : Nat) :
    ((transferRatNeutral b w : ℚ) : ℝ) =
      (1 - 2 * (b : ℝ)) ^ w := by
  simp [transferRatNeutral]

theorem ratSpectrum_cast (f : Nat → ℚ) :
    ((transferRatSpectrum f : ℚ) : ℝ) =
      productionSpectrumPolynomial (fun w => ((f w : ℚ) : ℝ)) := by
  simp [transferRatSpectrum, transferSpectrumTerms,
    productionSpectrumPolynomial]

theorem ratH_cast (b z : ℚ) :
    ((transferRatH b z : ℚ) : ℝ) = spectrumH (b : ℝ) (z : ℝ) := by
  rw [spectrumH_polynomial]
  simp [transferRatH, ratSpectrum_cast, ratFeedback_cast]
  norm_num

theorem ratG_cast (b z : ℚ) :
    ((transferRatG b z : ℚ) : ℝ) = spectrumG (b : ℝ) (z : ℝ) := by
  simp [transferRatG, spectrumG, ratG0_cast, ratH_cast]

theorem ratU_cast (b z : ℚ) :
    ((transferRatU b z : ℚ) : ℝ) = spectrumU (b : ℝ) (z : ℝ) := by
  rw [spectrumU_polynomial, state22_liveCount]
  simp [transferRatU, ratSpectrum_cast, ratWeight_cast]

theorem ratu_cast (b : ℚ) :
    ((transferRatu b : ℚ) : ℝ) = spectrumu (b : ℝ) := by
  rw [spectrumu_polynomial]
  simp [transferRatu, ratSpectrum_cast, ratNeutral_cast]
  norm_num

theorem ratF_cast (b z : ℚ)
    (hb0 : 0 ≤ (b : ℝ)) (hb1 : (b : ℝ) ≤ 1)
    (hz0 : 0 ≤ (z : ℝ)) :
    ((transferRatF b z : ℚ) : ℝ) = spectrumF (b : ℝ) (z : ℝ) := by
  simp [transferRatF, transferLiveWeights, spectrumF,
    scalableSpectrumMaxMoment, liveSpectrumWeights, spectrumWeights,
    ratWeight_cast, max_assoc]
  left
  unfold productionWeightMoment
  have hg0 : 0 ≤ 1 - (b : ℝ) + (b : ℝ) * (z : ℝ) := by
    nlinarith [mul_nonneg hb0 hz0]
  have hg1 : 0 ≤ (b : ℝ) + (1 - (b : ℝ)) * (z : ℝ) := by
    have : 0 ≤ 1 - (b : ℝ) := by linarith
    positivity
  positivity

theorem ratRows_cast (b z : ℚ)
    (hb0 : 0 ≤ (b : ℝ)) (hb1 : (b : ℝ) ≤ 1)
    (hz0 : 0 ≤ (z : ℝ)) (v : Fin 3 → ℚ) (i : Fin 3) :
    ((transferRatRows b z v i : ℚ) : ℝ) =
      (paperSpectrumTransfer (b : ℝ) (z : ℝ)).mulVec
        (fun j => (v j : ℝ)) i := by
  rw [paperSpectrumTransfer_components]
  fin_cases i <;>
    simp [transferRatRows, ratH_cast, ratG_cast, ratU_cast, ratu_cast,
      ratF_cast b z hb0 hb1 hz0, state22_liveCount]

theorem transferRatCheck_sound (b z ρ : ℚ) (v : Fin 3 → ℚ)
    (hb0 : 0 ≤ (b : ℝ)) (hb1 : (b : ℝ) ≤ 1)
    (hz0 : 0 ≤ (z : ℝ))
    (hcheck : transferRatCheck b z ρ v = true) :
    ∀ i, (paperSpectrumTransfer (b : ℝ) (z : ℝ)).mulVec
      (fun j => (v j : ℝ)) i ≤ (ρ : ℝ) * (v i : ℝ) := by
  simp only [transferRatCheck, Bool.and_eq_true, decide_eq_true_eq]
    at hcheck
  obtain ⟨⟨h0, h1⟩, h2⟩ := hcheck
  have hq (i : Fin 3) : transferRatRows b z v i ≤ ρ * v i := by
    fin_cases i
    · simpa using h0
    · simpa using h1
    · simpa using h2
  intro i
  rw [← ratRows_cast b z hb0 hb1 hz0 v i]
  exact_mod_cast hq i

theorem transferQInputRat_cast
    (a : Spin.Structured.DenseOccupationFixed.QInput) :
    ((transferQInputRat a : ℚ) : ℝ) = a.real := by
  simp [transferQInputRat,
    Spin.Structured.DenseOccupationFixed.QInput.real]

#print axioms transferRatCheck_sound

end Spin.UniformMixer
