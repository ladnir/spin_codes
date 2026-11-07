import SpinCodes.UniformMixer.ScalableDualLowCoefficient
import SpinCodes.UniformMixer.ScalableFeedbackFourier
import SpinCodes.Structured.ConcreteFourierProduct

/-!
# Dual low-weight certificate for the actual feedback map

The existing Krawtchouk character-sum theorem counts feedback-kernel words
of each weight from the certified production spectrum. This identifies the
checked low transform coefficients with actual kernel counts.
-/

noncomputable section

namespace Spin.UniformMixer

def feedbackKernelLayer (j : ℕ) : Finset (Finset (Fin 64)) :=
  (Spin.layer 64 j).filter (fun S => stepFeedback (supportInput S) = 0)

theorem feedback_character_chi (q : State22) (S : Finset (Fin 64)) :
    binarySign (statePairing q (stepFeedback (supportInput S))) =
      Spin.chi (expansionSupport q) S := by
  rw [feedback_character_support]
  apply Int.cast_injective (α := ℝ)
  push_cast
  rw [Spin.Structured.ConcreteFourier.chi_prod]
  apply Finset.prod_congr rfl
  intro j _
  by_cases hj : coordMap j q = 1
  · simp [expansionSupport, hj, binarySign_one]
  · have hz := (zmod2_zero_or_one (coordMap j q)).resolve_right hj
    simp [expansionSupport, hz, binarySign_zero]

theorem feedbackKernelLayer_fourier (j : ℕ) :
    (2 ^ 22 : ℤ) * ((feedbackKernelLayer j).card : ℤ) =
      ∑ q : State22,
        Spin.krawtchouk 64 j (coordinateWeight coordMap q) := by
  apply Int.cast_injective (α := ℝ)
  push_cast
  calc
    4194304 * (feedbackKernelLayer j).card =
        ∑ S ∈ Spin.layer 64 j,
          (4194304 : ℝ) *
            (if stepFeedback (supportInput S) = 0 then 1 else 0) := by
      simp [feedbackKernelLayer, Finset.sum_ite]
      ring
    _ = ∑ S ∈ Spin.layer 64 j,
          ∑ q : State22,
            (binarySign
              (statePairing q (stepFeedback (supportInput S))) : ℝ) := by
      apply Finset.sum_congr rfl
      intro S _
      simpa only [show (2 ^ 22 : ℝ) = 4194304 by norm_num] using
        feedback_indicator_character (supportInput S)
    _ = ∑ q : State22,
          ∑ S ∈ Spin.layer 64 j,
            (Spin.chi (expansionSupport q) S : ℝ) := by
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro q _
      apply Finset.sum_congr rfl
      intro S _
      rw [feedback_character_chi]
    _ = ∑ q : State22,
          (Spin.krawtchouk 64 j (coordinateWeight coordMap q) : ℝ) := by
      apply Finset.sum_congr rfl
      intro q _
      rw [← expansionSupport_card q,
        ← Spin.sum_chi_layer (expansionSupport q) j]
      push_cast
      rfl

theorem feedbackKernelLayer_spectrum (j : ℕ) :
    (2 ^ 22 : ℤ) * ((feedbackKernelLayer j).card : ℤ) =
      dualSpectrumNumerator j := by
  rw [feedbackKernelLayer_fourier]
  apply Int.cast_injective (α := ℝ)
  push_cast
  calc
    (∑ q : State22,
      (Spin.krawtchouk 64 j (coordinateWeight coordMap q) : ℝ)) =
        ∑ w ∈ spectrumWeights,
          (productionWeightCount w : ℝ) *
            (Spin.krawtchouk 64 j w : ℝ) :=
      productionSpectrum_sum (fun w => (Spin.krawtchouk 64 j w : ℝ))
    _ = ∑ w ∈ spectrumWeights,
          (productionSpectrumMultiplicity w : ℝ) *
            (Spin.krawtchouk 64 j w : ℝ) := by
      apply Finset.sum_congr rfl
      intro w _
      rw [productionWeightCount_exact]
    _ = (dualSpectrumNumerator j : ℝ) := by
      simp only [dualSpectrumNumerator, Int.cast_sum, Int.cast_mul,
        Int.cast_natCast]

theorem feedbackKernelLayer_low (j : Fin 8) (hj : (j : ℕ) ≠ 0) :
    (feedbackKernelLayer j).card = 0 := by
  have h := feedbackKernelLayer_spectrum (j : ℕ)
  rw [dualSpectrumNumerator_low j] at h
  simp only [hj, ite_false] at h
  omega

theorem feedbackKernel_distance (S : Finset (Fin 64))
    (hzero : stepFeedback (supportInput S) = 0)
    (hne : S.Nonempty) : 8 ≤ S.card := by
  by_contra hsmall
  have hpos : 0 < S.card := Finset.card_pos.mpr hne
  have hj8 : S.card < 8 := by omega
  let j : Fin 8 := ⟨S.card, hj8⟩
  have hmem : S ∈ feedbackKernelLayer j := by
    simp [feedbackKernelLayer, Spin.layer, hzero, j]
  have hempty : feedbackKernelLayer j = ∅ :=
    Finset.card_eq_zero.mp (feedbackKernelLayer_low j (by dsimp [j]; omega))
  rw [hempty] at hmem
  simp at hmem

end Spin.UniformMixer
