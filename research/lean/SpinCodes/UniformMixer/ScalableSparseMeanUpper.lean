import SpinCodes.UniformMixer.ScalableSparseMeanPointwise

/-! The average live-row bound uses the production map's exact mean weight. -/

noncomputable section

namespace Spin.UniformMixer

private theorem liveSet_card_generic {S : Type*}
    [Fintype S] [DecidableEq S] [AddCommGroup S] [Nontrivial S] :
    ((Finset.univ : Finset S).erase 0).card = liveCount S := by
  rw [Finset.card_erase_of_mem (Finset.mem_univ (0 : S)),
    Finset.card_univ]
  rfl

set_option maxRecDepth 8192 in
theorem scalableMeanLiveMoment_sparse_upper (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    scalableMeanLiveMoment β z hβ0 hβ1 ≤
      1 - 32 * ((1 - β) * (1 - z)) +
        2016 * ((1 - β) * (1 - z)) ^ 2 := by
  let t := (1 - β) * (1 - z)
  let D : ℝ := (liveCount State22 : ℝ)
  have ht0 : 0 ≤ t := mul_nonneg (by linarith) (by linarith)
  have hDpos : 0 < D := state22_liveCount_real_pos
  have hcard : (liveStateSet.card : ℝ) = D := by
    exact congrArg (fun n : ℕ => (n : ℝ))
      (liveSet_card_generic (S := State22))
  have hmean :
      32 ≤ (∑ q ∈ liveStateSet,
        (coordinateWeight coordMap q : ℝ)) / D := by
    have hnative : 32 ≤ liveWeightAverage coordMap := by
      rw [production_liveWeightAverage]
      norm_num
    rw [liveWeightAverage_erase] at hnative
    exact hnative
  have hsum :
      (∑ q ∈ liveStateSet, scalableLiveMoment β z hβ0 hβ1 q) ≤
        ∑ q ∈ liveStateSet,
          (1 - (coordinateWeight coordMap q : ℝ) * t +
            2016 * t ^ 2) := by
    apply Finset.sum_le_sum
    intro q _
    exact scalableLiveMoment_sparse_quadratic β z hβ0 hβ1 hz0 hz1 q
  have hsumform :
      (∑ q ∈ liveStateSet,
          (1 - (coordinateWeight coordMap q : ℝ) * t +
            2016 * t ^ 2)) =
        D * (1 + 2016 * t ^ 2) -
          (∑ q ∈ liveStateSet,
            (coordinateWeight coordMap q : ℝ)) * t := by
    calc
      _ = ∑ q ∈ liveStateSet,
            ((1 + 2016 * t ^ 2) -
              (coordinateWeight coordMap q : ℝ) * t) := by
        apply Finset.sum_congr rfl
        intro q _
        ring
      _ = _ := by
        rw [Finset.sum_sub_distrib, ← Finset.sum_mul]
        rw [Finset.sum_const, nsmul_eq_mul, hcard]
  have hdivide :
      (∑ q ∈ liveStateSet,
        scalableLiveMoment β z hβ0 hβ1 q) / D ≤
        1 - ((∑ q ∈ liveStateSet,
          (coordinateWeight coordMap q : ℝ)) / D) * t +
          2016 * t ^ 2 := by
    have h := div_le_div_of_nonneg_right hsum hDpos.le
    rw [hsumform] at h
    have heq :
        (D * (1 + 2016 * t ^ 2) -
          (∑ q ∈ liveStateSet,
            (coordinateWeight coordMap q : ℝ)) * t) / D =
          1 - ((∑ q ∈ liveStateSet,
            (coordinateWeight coordMap q : ℝ)) / D) * t +
            2016 * t ^ 2 := by
      field_simp [ne_of_gt hDpos]
      ring
    rwa [heq] at h
  unfold scalableMeanLiveMoment
  have hmt := mul_le_mul_of_nonneg_right hmean ht0
  change (∑ q ∈ liveStateSet,
      scalableLiveMoment β z hβ0 hβ1 q) / D ≤
      1 - 32 * t + 2016 * t ^ 2
  linarith only [hdivide, hmt]

end Spin.UniformMixer
