import SpinCodes.UniformMixer.ScalableSparseBounds

/-! The minimum nonzero production weight bounds the fresh-live row. -/

noncomputable section

namespace Spin.UniformMixer

theorem productionWeightMoment_le_minimum (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (w : ℕ) (hw : 16 ≤ w) :
    productionWeightMoment β z w ≤
      (β + (1 - β) * z) ^ 16 := by
  have hg0nonneg : 0 ≤ 1 - β + β * z := by nlinarith
  have hg0one : 1 - β + β * z ≤ 1 := by nlinarith
  have hg1nonneg : 0 ≤ β + (1 - β) * z := by nlinarith
  have hg1one : β + (1 - β) * z ≤ 1 := by nlinarith
  have hp0 : (1 - β + β * z) ^ (64 - w) ≤ 1 :=
    pow_le_one₀ hg0nonneg hg0one
  have hp1 : (β + (1 - β) * z) ^ w ≤
      (β + (1 - β) * z) ^ 16 :=
    pow_le_pow_of_le_one hg1nonneg hg1one hw
  unfold productionWeightMoment
  calc
    (1 - β + β * z) ^ (64 - w) *
        (β + (1 - β) * z) ^ w ≤
      1 * (β + (1 - β) * z) ^ w := by
      exact mul_le_mul_of_nonneg_right hp0 (pow_nonneg hg1nonneg _)
    _ = (β + (1 - β) * z) ^ w := one_mul _
    _ ≤ (β + (1 - β) * z) ^ 16 := hp1

theorem spectrumF_sparse_upper (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    spectrumF β z ≤
      1 - 16 * ((1 - β) * (1 - z)) +
        120 * ((1 - β) * (1 - z)) ^ 2 := by
  have ht0 : 0 ≤ (1 - β) * (1 - z) :=
    mul_nonneg (by linarith) (by linarith)
  have ht1 : (1 - β) * (1 - z) ≤ 1 := by nlinarith
  have hchoose : (Nat.choose 16 2 : ℝ) = 120 := by
    exact_mod_cast (show Nat.choose 16 2 = 120 by decide)
  have hquad := one_sub_pow_le_quadratic 16
    ((1 - β) * (1 - z)) ht0 ht1
  have hquad' :
      (β + (1 - β) * z) ^ 16 ≤
        1 - 16 * ((1 - β) * (1 - z)) +
          120 * ((1 - β) * (1 - z)) ^ 2 := by
    have heq : β + (1 - β) * z =
        1 - (1 - β) * (1 - z) := by ring
    rw [heq]
    simpa only [hchoose, Nat.cast_ofNat] using hquad
  unfold spectrumF scalableSpectrumMaxMoment
  apply Finset.sup'_le
  intro w hw
  obtain ⟨q, hq, hweight⟩ := productionWeight_exists_of_mem_live w hw
  have hmin : 16 ≤ w := by
    rw [← hweight]
    exact productionWeight_nonzero_lower_bound q hq
  exact (productionWeightMoment_le_minimum β z hβ0 hβ1 hz0 hz1 w hmin).trans hquad'

end Spin.UniformMixer
