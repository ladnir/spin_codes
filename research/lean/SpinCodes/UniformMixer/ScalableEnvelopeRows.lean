import SpinCodes.UniformMixer.ScalableEnvelope

/-!
# Row parameters for the actual Bernoulli mixer

The paper's three-coordinate matrix uses the largest weighted live
row, its average over nonzero entering states, and the probability
of nonzero feedback. These definitions keep the actual production
maps and the full state space.
-/

noncomputable section

namespace Spin.UniformMixer

private theorem live_sum_constant_generic {S : Type*}
    [Fintype S] [DecidableEq S] [AddCommGroup S] [Nontrivial S]
    (u : ℝ) :
    ∑ _q ∈ (Finset.univ : Finset S).erase 0,
      u / (liveCount S : ℝ) = u := by
  have hc : ((Finset.univ : Finset S).erase 0).card =
      liveCount S := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ (0 : S)),
      Finset.card_univ]
    rfl
  rw [Finset.sum_const, nsmul_eq_mul, hc]
  exact mul_div_cancel₀ u
    (Nat.cast_ne_zero.mpr (Nat.ne_of_gt (liveCount_pos (S := S))))

theorem liveStateSet_sum_constant (u : ℝ) :
    ∑ _q ∈ liveStateSet, u / (liveCount State22 : ℝ) = u :=
  live_sum_constant_generic (S := State22) u

def scalableMaxLiveMoment (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) : ℝ :=
  liveStateSet.sup' liveStateSet_nonempty
    (scalableLiveMoment β z hβ0 hβ1)

def scalableMeanLiveMoment (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) : ℝ :=
  (∑ q ∈ liveStateSet, scalableLiveMoment β z hβ0 hβ1 q) /
    (liveCount State22 : ℝ)

theorem scalableLiveMoment_nonneg (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z)
    (q : State22) :
    0 ≤ scalableLiveMoment β z hβ0 hβ1 q := by
  exact (bernoulliInputLaw β hβ0 hβ1).expect_nonneg
    (fun x => pow_nonneg hz _)

theorem scalableMaxLiveMoment_nonneg (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) :
    0 ≤ scalableMaxLiveMoment β z hβ0 hβ1 := by
  obtain ⟨q, hq⟩ := liveStateSet_nonempty
  exact (scalableLiveMoment_nonneg β z hβ0 hβ1 hz q).trans
    (Finset.le_sup' _ hq)

theorem scalableLiveMoment_le_max (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (q : State22) (hq : q ≠ 0) :
    scalableLiveMoment β z hβ0 hβ1 q ≤
      scalableMaxLiveMoment β z hβ0 hβ1 :=
  Finset.le_sup' _ (Finset.mem_erase.mpr
    ⟨hq, Finset.mem_univ _⟩)

theorem scalableMeanLiveMoment_nonneg (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) :
    0 ≤ scalableMeanLiveMoment β z hβ0 hβ1 := by
  apply div_nonneg
  · exact Finset.sum_nonneg (fun q _ =>
      scalableLiveMoment_nonneg β z hβ0 hβ1 hz q)
  · exact state22_liveCount_real_pos.le

theorem bernoulliStepKernel_nonneg (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z)
    (q r : State22) :
    0 ≤ bernoulliStepKernel β z hβ0 hβ1 q r := by
  change 0 ≤ (bernoulliInputLaw β hβ0 hβ1).expect
    (fun x => stepKernel z x q r)
  apply FinPMF.expect_nonneg
  intro x
  exact mul_nonneg (pow_nonneg hz _)
    (update_nonneg (stepFeedback x) q r)

theorem scalable_actual_weighted_rows (β z : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1)
    (hz0 : 0 ≤ z) (hz1 : z ≤ 1) :
    WeightedRowBounds (bernoulliStepKernel β z hβ0 hβ1)
      (scalableZeroHold β z hβ0 hβ1)
      (scalableLiveMoment β z hβ0 hβ1 0 -
        scalableZeroHold β z hβ0 hβ1)
      (scalableMaxLiveMoment β z hβ0 hβ1)
      (scalableMeanLiveMoment β z hβ0 hβ1)
      (min (scalableMaxLiveMoment β z hβ0 hβ1)
        (scalableFeedbackActive β hβ0 hβ1))
      (min (scalableMeanLiveMoment β z hβ0 hβ1)
        (scalableFeedbackActive β hβ0 hβ1)) := by
  let K := bernoulliStepKernel β z hβ0 hβ1
  let F := scalableMaxLiveMoment β z hβ0 hβ1
  let U := scalableMeanLiveMoment β z hβ0 hβ1
  let u := scalableFeedbackActive β hβ0 hβ1
  have hd : (0 : ℝ) ≤ (liveCount State22 : ℝ) :=
    state22_liveCount_real_pos.le
  have hK : ∀ q r, 0 ≤ K q r :=
    bernoulliStepKernel_nonneg β z hβ0 hβ1 hz0
  have hF : 0 ≤ F := scalableMaxLiveMoment_nonneg β z hβ0 hβ1 hz0
  have hU : 0 ≤ U := scalableMeanLiveMoment_nonneg β z hβ0 hβ1 hz0
  have hu : 0 ≤ u := (bernoulliInputLaw β hβ0 hβ1).prob_nonneg _
  have hmean (r : State22) :
      ∑ q ∈ liveStateSet, K q r ≤ U := by
    calc
      _ ≤ ∑ q ∈ liveStateSet,
          scalableLiveMoment β z hβ0 hβ1 q /
            (liveCount State22 : ℝ) := by
        apply Finset.sum_le_sum
        intro q hq
        exact bernoulliStepKernel_live_atom_le β z hβ0 hβ1 hz0
          q r (Finset.ne_of_mem_erase hq)
      _ = U := by
        simp only [U, scalableMeanLiveMoment, Finset.sum_div]
  refine ⟨hK, ?_, ?_, hF, hU, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [← bernoulliStepKernel_zero_zero β z hβ0 hβ1]
    exact hK 0 0
  · rw [← bernoulliStepKernel_zero_live_mass β z hβ0 hβ1]
    exact Finset.sum_nonneg (fun r _ => hK 0 r)
  · exact le_min hF hu
  · exact le_min hU hu
  · exact (bernoulliStepKernel_zero_zero β z hβ0 hβ1).le
  · exact (bernoulliStepKernel_zero_live_mass β z hβ0 hβ1).le
  · intro q r hq _hr
    exact (bernoulliStepKernel_live_atom_le β z hβ0 hβ1 hz0 q r hq).trans
      (div_le_div_of_nonneg_right
        (scalableLiveMoment_le_max β z hβ0 hβ1 q hq) hd)
  · intro q hq
    rw [← min_div_div_right hd]
    apply le_min
    · exact (bernoulliStepKernel_live_atom_le β z hβ0 hβ1 hz0 q 0 hq).trans
        (div_le_div_of_nonneg_right
          (scalableLiveMoment_le_max β z hβ0 hβ1 q hq) hd)
    · exact bernoulliStepKernel_live_return_le β z hβ0 hβ1
        hz0 hz1 q hq
  · exact fun r _ => hmean r
  · apply le_min
    · exact hmean 0
    · calc
        _ ≤ ∑ _q ∈ liveStateSet,
            u / (liveCount State22 : ℝ) := by
          apply Finset.sum_le_sum
          intro q hq
          exact bernoulliStepKernel_live_return_le β z hβ0 hβ1
            hz0 hz1 q (Finset.ne_of_mem_erase hq)
        _ = u := liveStateSet_sum_constant u

end Spin.UniformMixer
