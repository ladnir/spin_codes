import SpinCodes.UniformMixer.QuadraticSpectrumSupport
import SpinCodes.UniformMixer.ScalableMomentSpectrum

/-!
# Exact production-spectrum summation interface

This transfers the checked multiplicities to any real-valued function of
the 64-coordinate output weight. It is the bridge needed for the mean
Bernoulli row moment in the Scalable SPIN transfer matrix.
-/

noncomputable section

namespace Spin.UniformMixer

theorem productionSpectrum_sum (f : Nat → ℝ) :
    (∑ q : State22, f (coordinateWeight coordMap q)) =
      ∑ w ∈ spectrumWeights, (productionWeightCount w : ℝ) * f w := by
  calc
    (∑ q : State22, f (coordinateWeight coordMap q)) =
        ∑ q : State22,
          ∑ w ∈ spectrumWeights,
            if coordinateWeight coordMap q = w then f w else 0 := by
      apply Finset.sum_congr rfl
      intro q _
      have hmem := productionWeight_mem_spectrumWeights q
      simp only [eq_comm, Finset.sum_ite_eq', hmem, ite_true]
    _ = ∑ w ∈ spectrumWeights,
          ∑ q : State22,
            if coordinateWeight coordMap q = w then f w else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ w ∈ spectrumWeights,
          (productionWeightCount w : ℝ) * f w := by
      apply Finset.sum_congr rfl
      intro w _
      rw [Finset.sum_ite]
      simp only [Finset.sum_const_zero, add_zero, Finset.sum_const,
        nsmul_eq_mul]
      rw [productionWeightCount_eq_card_filter]

def productionWeightMoment (β z : ℝ) (w : Nat) : ℝ :=
  (1 - β + β * z) ^ (64 - w) * (β + (1 - β) * z) ^ w

theorem scalableMeanLiveMoment_exact_spectrum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    scalableMeanLiveMoment β z hβ0 hβ1 =
      ((∑ w ∈ spectrumWeights,
        (productionWeightCount w : ℝ) * productionWeightMoment β z w) -
        productionWeightMoment β z 0) /
      (liveCount State22 : ℝ) := by
  have hpoint (q : State22) :
      scalableLiveMoment β z hβ0 hβ1 q =
        productionWeightMoment β z (coordinateWeight coordMap q) := by
    exact scalableLiveMoment_weight_formula β z hβ0 hβ1 q
  have htotal :
      (∑ q : State22, scalableLiveMoment β z hβ0 hβ1 q) =
        ∑ w ∈ spectrumWeights,
          (productionWeightCount w : ℝ) * productionWeightMoment β z w := by
    calc
      _ = ∑ q : State22,
          productionWeightMoment β z (coordinateWeight coordMap q) := by
        apply Finset.sum_congr rfl
        intro q _
        exact hpoint q
      _ = _ := productionSpectrum_sum (productionWeightMoment β z)
  have hsplit := Finset.add_sum_erase (Finset.univ : Finset State22)
    (fun q => scalableLiveMoment β z hβ0 hβ1 q)
    (Finset.mem_univ (0 : State22))
  change scalableLiveMoment β z hβ0 hβ1 0 +
      (∑ q ∈ liveStateSet, scalableLiveMoment β z hβ0 hβ1 q) =
        ∑ q : State22, scalableLiveMoment β z hβ0 hβ1 q at hsplit
  rw [htotal, hpoint 0, coordinateWeight_zero] at hsplit
  unfold scalableMeanLiveMoment
  rw [show (∑ q ∈ liveStateSet, scalableLiveMoment β z hβ0 hβ1 q) =
      (∑ w ∈ spectrumWeights,
        (productionWeightCount w : ℝ) * productionWeightMoment β z w) -
        productionWeightMoment β z 0 from by linarith]

end Spin.UniformMixer
