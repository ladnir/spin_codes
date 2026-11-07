import SpinCodes.UniformMixer.ScalableTransferSpectrumModel

/-!
# Explicit nine-term polynomial for the production spectrum

This removes the finite weight-index sum from the transfer formulas. The
result is an arithmetic expression in nine evaluations of its argument,
with coefficients already proved to be the actual production counts.
-/

noncomputable section

namespace Spin.UniformMixer

def productionSpectrumPolynomial (f : Nat → ℝ) : ℝ :=
  f 0 + 2604 * f 16 + 291648 * f 24 + 888832 * f 28 +
    1828134 * f 32 + 888832 * f 36 + 291648 * f 40 +
    2604 * f 48 + f 64

theorem productionSpectrum_weighted_sum (f : Nat → ℝ) :
    (∑ w ∈ spectrumWeights,
      (productionWeightCount w : ℝ) * f w) =
        productionSpectrumPolynomial f := by
  calc
    (∑ w ∈ spectrumWeights,
      (productionWeightCount w : ℝ) * f w) =
        ∑ w ∈ spectrumWeights,
          (productionSpectrumMultiplicity w : ℝ) * f w := by
      apply Finset.sum_congr rfl
      intro w _
      rw [productionWeightCount_exact]
    _ = productionSpectrumPolynomial f := by
      simp [spectrumWeights, productionSpectrumMultiplicity,
        productionSpectrumPolynomial]
      ring

theorem spectrumU_polynomial (β z : ℝ) :
    spectrumU β z =
      (productionSpectrumPolynomial (productionWeightMoment β z) -
        productionWeightMoment β z 0) / (liveCount State22 : ℝ) := by
  unfold spectrumU
  rw [productionSpectrum_weighted_sum]

theorem spectrumH_polynomial (β z : ℝ) :
    spectrumH β z =
      productionSpectrumPolynomial (feedbackWeightMoment β z) /
        (2 ^ 22 : ℝ) := by
  unfold spectrumH
  rw [productionSpectrum_weighted_sum]

theorem spectrumu_polynomial (β : ℝ) :
    spectrumu β =
      1 - productionSpectrumPolynomial (fun w => (1 - 2 * β) ^ w) /
        (2 ^ 22 : ℝ) := by
  unfold spectrumu
  rw [productionSpectrum_weighted_sum]

end Spin.UniformMixer
