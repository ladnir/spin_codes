import SpinCodes.UniformMixer.ScalableFixedLargeTilt
import SpinCodes.Structured.ConcreteNativeFixedLargeCount

/-! Reuse the selected outer's shell count with a changed inner tilt. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily
  Spin.Structured.ConcreteFixedNumeric
  Spin.Structured.Placement Finset

theorem scalable_tilted_counted_profile_bound
    {Q k L d : ℕ} (hk : 0 < k)
    (seed : Seed k) (W : Finset ℕ)
    (hg : Spin.Good (seedLaw k) spectrum (k * 24) W seed)
    (rows : Fin Q → Finset (Fin (k * 24)))
    (u : FixedFugacityChoice Q)
    (hw : ∀ i, (rows i).card ∈ W)
    (hr : ∀ i, ((rows i).card : ℝ) /
      (k * 24 : ℕ) ∈ Set.Icc (13 / 125) (112 / 125))
    (hu : ∀ i, rate (((rows i).card : ℝ) /
      (k * 24 : ℕ)) (fixedFugacity u i) ≤
        -(8679 / 10000000))
    (θNew : ℝ) {P : ℝ}
    (hP : P ≤
      (((fixedNormCost u + fixedNormSlack u) ^
            (k * 24) / min 1 (3 / 1600 : ℝ)) /
          (profileChoices rows *
            profileCoefficient (fixedFugacity u)
              (fun i => (rows i).card))) /
        Real.exp (-(θNew / (L : ℝ))) ^ d) :
    (∏ i, (spectrum seed (rows i).card : ℝ)) * P ≤
      fixedProfileCountBound Q (k * 24) L d *
        Real.exp (((θNew - paperLargeTheta Q) /
          (L : ℝ)) * (d : ℝ)) := by
  let A : ℝ :=
    ((fixedNormCost u + fixedNormSlack u) ^
      (k * 24) / min 1 (3 / 1600 : ℝ)) /
        (profileChoices rows *
          profileCoefficient (fixedFugacity u)
            (fun i => (rows i).card))
  let Pold : ℝ := A /
    Real.exp (-(paperLargeTheta Q / (L : ℝ))) ^ d
  let correction : ℝ :=
    Real.exp (((θNew - paperLargeTheta Q) /
      (L : ℝ)) * (d : ℝ))
  have hfactor : A / Real.exp (-(θNew / (L : ℝ))) ^ d =
      Pold * correction := by
    simpa only [Pold, correction, neg_div] using
      tilt_div_identity A θNew (paperLargeTheta Q) (L : ℝ) d
  have hOld :
      (∏ i, (spectrum seed (rows i).card : ℝ)) * Pold ≤
        fixedProfileCountBound Q (k * 24) L d := by
    apply fixed_counted_profile_bound hk seed W hg rows u
      hw hr hu
    change Pold ≤ Pold
    exact le_rfl
  calc
    _ ≤ (∏ i, (spectrum seed (rows i).card : ℝ)) *
        (Pold * correction) := by
      apply mul_le_mul_of_nonneg_left _
        (prod_nonneg fun _ _ => Nat.cast_nonneg _)
      exact hP.trans (by simpa only [A] using hfactor.le)
    _ = ((∏ i, (spectrum seed (rows i).card : ℝ)) *
        Pold) * correction := by ring
    _ ≤ _ := mul_le_mul_of_nonneg_right hOld
      (Real.exp_pos _).le

end Spin.UniformMixer
