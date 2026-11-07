import SpinCodes.UniformMixer.ScalableSparseLiveUpper

/-! Quadratic upper bound for each actual live state. -/

noncomputable section

namespace Spin.UniformMixer

theorem scalableLiveMoment_sparse_quadratic (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1)
    (q : State22) :
    scalableLiveMoment β z hβ0 hβ1 q ≤
      1 - (coordinateWeight coordMap q : ℝ) *
        ((1 - β) * (1 - z)) +
        2016 * ((1 - β) * (1 - z)) ^ 2 := by
  let w := coordinateWeight coordMap q
  let t := (1 - β) * (1 - z)
  have ht0 : 0 ≤ t := mul_nonneg (by linarith) (by linarith)
  have ht1 : t ≤ 1 := by dsimp [t]; nlinarith
  have hg0nonneg : 0 ≤ 1 - β + β * z := by nlinarith
  have hg0one : 1 - β + β * z ≤ 1 := by nlinarith
  have hg1nonneg : 0 ≤ β + (1 - β) * z := by nlinarith
  have hbase :
      productionWeightMoment β z w ≤ (1 - t) ^ w := by
    have hp0 : (1 - β + β * z) ^ (64 - w) ≤ 1 :=
      pow_le_one₀ hg0nonneg hg0one
    have heq : β + (1 - β) * z = 1 - t := by dsimp [t]; ring
    unfold productionWeightMoment
    rw [← heq]
    exact (mul_le_mul_of_nonneg_right hp0 (pow_nonneg hg1nonneg _)).trans_eq
      (one_mul _)
  have hw64 : w ≤ 64 := coordinateWeight_le coordMap q
  have hchoose : (Nat.choose w 2 : ℝ) ≤ 2016 := by
    exact_mod_cast (Nat.choose_le_choose 2 hw64).trans
      (le_of_eq (show Nat.choose 64 2 = 2016 by decide))
  have hquad := one_sub_pow_le_quadratic w t ht0 ht1
  have hsq : 0 ≤ t ^ 2 := sq_nonneg t
  have htail : (Nat.choose w 2 : ℝ) * t ^ 2 ≤
      2016 * t ^ 2 := mul_le_mul_of_nonneg_right hchoose hsq
  have hfull : productionWeightMoment β z w ≤
      1 - (w : ℝ) * t + 2016 * t ^ 2 := by
    exact hbase.trans (hquad.trans (by linarith only [htail]))
  rw [scalableLiveMoment_weight_formula]
  exact hfull

end Spin.UniformMixer
