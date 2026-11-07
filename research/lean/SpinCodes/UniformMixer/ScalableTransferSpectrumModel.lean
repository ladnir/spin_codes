import SpinCodes.UniformMixer.ScalableLiveMaxSpectrum
import SpinCodes.UniformMixer.ScalableTransferFirstMoment

/-!
# Finite spectrum model of the actual Scalable SPIN transfer

All entries in the three-coordinate transfer now depend on the checked
64-coordinate production spectrum, a maximum over eight weights, and
elementary real arithmetic. This model is the target for the sparse and
dense numerical certificates.
-/

noncomputable section

namespace Spin.UniformMixer

def spectrumF (β z : ℝ) : ℝ := scalableSpectrumMaxMoment β z

def spectrumU (β z : ℝ) : ℝ :=
  ((∑ w ∈ spectrumWeights,
      (productionWeightCount w : ℝ) * productionWeightMoment β z w) -
    productionWeightMoment β z 0) / (liveCount State22 : ℝ)

def spectrumH (β z : ℝ) : ℝ :=
  (∑ w ∈ spectrumWeights,
    (productionWeightCount w : ℝ) * feedbackWeightMoment β z w) /
    (2 ^ 22 : ℝ)

def spectrumu (β : ℝ) : ℝ :=
  1 - (∑ w ∈ spectrumWeights,
    (productionWeightCount w : ℝ) * (1 - 2 * β) ^ w) /
    (2 ^ 22 : ℝ)

def spectrumG (β z : ℝ) : ℝ :=
  (1 - β + β * z) ^ 64 - spectrumH β z

def spectrumTransfer (β z : ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  weightedTransfer (spectrumH β z) (spectrumG β z)
    (spectrumF β z) (spectrumU β z)
    (min (spectrumF β z) (spectrumu β))
    (min (spectrumU β z) (spectrumu β))

theorem scalableLiveMoment_zero_formula (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableLiveMoment β z hβ0 hβ1 0 =
      (1 - β + β * z) ^ 64 := by
  rw [scalableLiveMoment_weight_formula, coordinateWeight_zero]
  simp

theorem scalableThreeTransfer_eq_spectrum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableThreeTransfer β z hβ0 hβ1 = spectrumTransfer β z := by
  unfold scalableThreeTransfer spectrumTransfer spectrumF spectrumU
    spectrumH spectrumu spectrumG
  rw [scalableMaxLiveMoment_eq_spectrum β z hβ0 hβ1,
    scalableMeanLiveMoment_exact_spectrum β z hβ0 hβ1,
    scalableZeroHold_exact_spectrum β z hβ0 hβ1,
    scalableFeedbackActive_exact_spectrum β hβ0 hβ1,
    scalableLiveMoment_zero_formula β z hβ0 hβ1]
  rfl

end Spin.UniformMixer
