import SpinCodes.UniformMixer.ScalableSparseBounds

/-! The feedback-active bound uses the same empty-word lower bound as the
zero-to-live estimate. It is isolated to keep Lean's memory use bounded. -/

noncomputable section

namespace Spin.UniformMixer

theorem spectrumu_sparse_upper (β : ℝ)
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) :
    spectrumu β ≤ 64 * β := by
  have hH := scalableZeroHold_sparse_lower β 1 hβ0 hβ1 (by norm_num)
  have hbern : 1 - 64 * β ≤ (1 - β) ^ 64 := by
    have h := one_add_mul_le_pow (a := -β) (by linarith : -2 ≤ -β) 64
    simpa only [Nat.cast_ofNat, mul_neg, sub_eq_add_neg] using h
  calc
    spectrumu β = scalableFeedbackActive β hβ0 hβ1 :=
      (scalableFeedbackActive_exact_spectrum β hβ0 hβ1).symm
    _ = 1 - scalableZeroHold β 1 hβ0 hβ1 :=
      scalableFeedbackActive_complement β hβ0 hβ1
    _ ≤ 1 - (1 - β) ^ 64 := sub_le_sub_left hH _
    _ ≤ 64 * β := by linarith only [hbern]

end Spin.UniformMixer
