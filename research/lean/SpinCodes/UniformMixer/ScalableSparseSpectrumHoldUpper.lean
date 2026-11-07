import SpinCodes.UniformMixer.ScalableSparseBounds

/-! The sparse zero-to-zero estimate in the paper's transfer coordinates. -/

noncomputable section

namespace Spin.UniformMixer

theorem spectrumH_sparse_quadratic (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    spectrumH β z ≤
      1 - 64 * β + 2016 * β ^ 2 + (2 ^ 64 : ℝ) * β ^ 8 := by
  have heq : spectrumH β z = scalableZeroHold β z hβ0 hβ1 :=
    (scalableZeroHold_exact_spectrum β z hβ0 hβ1).symm
  rw [heq]
  exact scalableZeroHold_sparse_quadratic β z hβ0 hβ1 hz0 hz1

end Spin.UniformMixer
