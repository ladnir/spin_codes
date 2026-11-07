import SpinCodes.UniformMixer.QuadraticSpectrumSupport

/-!
# Low coefficients of the dual spectrum transform

This finite integer calculation checks the first eight coefficients of the
MacWilliams transform of the certified production spectrum. A separate
coefficient-extraction theorem is still needed to identify these integers
with the actual feedback-kernel weight counts.
-/

namespace Spin.UniformMixer

def dualKrawtchouk (w j : ℕ) : ℤ :=
  ∑ i ∈ Finset.range (j + 1),
    (if i % 2 = 0 then (1 : ℤ) else -1) *
      (Nat.choose w i : ℤ) * (Nat.choose (64 - w) (j - i) : ℤ)

def dualSpectrumNumerator (j : ℕ) : ℤ :=
  ∑ w ∈ spectrumWeights,
    (productionSpectrumMultiplicity w : ℤ) * dualKrawtchouk w j

theorem dualSpectrumNumerator_low (j : Fin 8) :
    dualSpectrumNumerator j =
      if (j : ℕ) = 0 then 2 ^ 22 else 0 := by
  fin_cases j <;> decide +kernel

end Spin.UniformMixer
