import SpinCodes.UniformMixer.QuadraticSpectrumSupport
import SpinCodes.Structured.Krawtchouk

/-!
# Low coefficients of the dual spectrum transform

This finite integer calculation checks the first eight coefficients of the
MacWilliams transform of the certified production spectrum. A separate
coefficient-extraction theorem is still needed to identify these integers
with the actual feedback-kernel weight counts.
-/

namespace Spin.UniformMixer

def dualSpectrumNumerator (j : ℕ) : ℤ :=
  ∑ w ∈ spectrumWeights,
    (productionSpectrumMultiplicity w : ℤ) * Spin.krawtchouk 64 j w

theorem dualSpectrumNumerator_low (j : Fin 8) :
    dualSpectrumNumerator j =
      if (j : ℕ) = 0 then 2 ^ 22 else 0 := by
  fin_cases j <;> decide +kernel

end Spin.UniformMixer
