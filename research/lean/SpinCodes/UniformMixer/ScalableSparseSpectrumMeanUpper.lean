import SpinCodes.UniformMixer.ScalableSparseMeanUpper

/-! Connect the live-state mean bound to the paper's spectrum transfer. -/

noncomputable section

namespace Spin.UniformMixer

theorem spectrumU_sparse_upper (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    spectrumU β z ≤
      1 - 32 * ((1 - β) * (1 - z)) +
        2016 * ((1 - β) * (1 - z)) ^ 2 := by
  have heq : spectrumU β z = scalableMeanLiveMoment β z hβ0 hβ1 :=
    (scalableMeanLiveMoment_exact_spectrum β z hβ0 hβ1).symm
  rw [heq]
  exact scalableMeanLiveMoment_sparse_upper β z hβ0 hβ1 hz0 hz1

end Spin.UniformMixer
