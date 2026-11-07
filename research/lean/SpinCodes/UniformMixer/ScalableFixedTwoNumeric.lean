import SpinCodes.UniformMixer.ScalableFixedTwoComparison
import SpinCodes.Structured.ConcreteFixedTwoBound

/-! Rational two-state margin for two fixed active outer rows. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement

def scalableThetaTwo : ℝ := 2 * scalableThetaOne

theorem regionGamma_thetaTwo : regionGamma scalableThetaTwo = 2 := by
  norm_num [regionGamma, production_liveWeightAverage,
    scalableThetaTwo, scalableThetaOne]

theorem scalable_fixed_two_row_gate :
    scalableTwoRowNormLe (1 / 1000) (101 / 100)
      (scalableFugacityTwo scalableThetaTwo
        (fun _ : Fin 2 => 1)) := by
  let γ : ℝ := regionGamma scalableThetaTwo
  have hγ : 2 ≤ γ := by
    dsimp [γ]
    rw [regionGamma_thetaTwo]
  have hγpos : 0 < γ := by linarith
  have he0 : 0 ≤ Real.exp (-γ) := (Real.exp_pos _).le
  have he : Real.exp (-γ) ≤ 1 / 4 :=
    exp_neg_le_quarter_of_two_le hγ
  have hf : (1 - Real.exp (-γ)) / γ ≤ 1 / 2 := by
    apply (div_le_iff₀ hγpos).mpr
    linarith
  have hg : 2 * (γ - 1 + Real.exp (-γ)) / γ ^ 2 ≤ 1 := by
    apply (div_le_iff₀ (sq_pos_of_pos hγpos)).mpr
    nlinarith [sq_nonneg (γ - 1)]
  have hd : 2 * (1 - (1 + γ) * Real.exp (-γ)) / γ ^ 2 ≤ 1 / 2 := by
    apply (div_le_iff₀ (sq_pos_of_pos hγpos)).mpr
    have hh : 0 ≤ (1 + γ) * Real.exp (-γ) :=
      mul_nonneg (by linarith) he0
    nlinarith [sq_nonneg (γ - 2)]
  have h00 := scalableContinuumTwo_two_le hγpos.ne' 0 0
  have h01 := scalableContinuumTwo_two_le hγpos.ne' 0 1
  have h10 := scalableContinuumTwo_two_le hγpos.ne' 1 0
  have h11 := scalableContinuumTwo_two_le hγpos.ne' 1 1
  change scalableContinuumTwo scalableThetaTwo 2 false false ≤
    (1 / 524287) * (2 * (γ - 1 + Real.exp (-γ)) / γ ^ 2) at h00
  change scalableContinuumTwo scalableThetaTwo 2 false true ≤
    2 * (1 - (1 + γ) * Real.exp (-γ)) / γ ^ 2 at h01
  change scalableContinuumTwo scalableThetaTwo 2 true false ≤
    (1 / 524287) *
      (2 * (1 - (1 + γ) * Real.exp (-γ)) / γ ^ 2) at h10
  change scalableContinuumTwo scalableThetaTwo 2 true true ≤
    (1 / 524287) *
      (2 * (1 - (1 + γ) * Real.exp (-γ)) / γ ^ 2) +
      Real.exp (-γ) at h11
  rw [scalableFugacityTwo_two, scalableContinuumTwo_zero,
    scalableContinuumTwo_one hγpos.ne']
  unfold scalableTwoRowNormLe
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul]
  dsimp only [timeRegionGap, coarseD, scalableOneFormula]
  simp only [Bool.false_eq_true, Bool.true_eq_false,
    ite_true, ite_false, mul_one, zero_add]
  constructor
  · nlinarith only [h00, h01, hf, hg, hd]
  · change 2 * (1 / (liveCount State22 : ℝ) *
        ((1 - Real.exp (-γ)) / γ)) +
        scalableContinuumTwo scalableThetaTwo 2 true false +
        (Real.exp (-γ) +
          2 * ((1 - 1 / (liveCount State22 : ℝ)) *
            Real.exp (-γ)) +
          scalableContinuumTwo scalableThetaTwo 2 true true) *
          (1 / 1000) ≤ (101 / 100) * (1 / 1000)
    rw [state22_liveCount]
    norm_num only [Nat.reducePow, Nat.reduceSub,
      Nat.cast_ofNat]
    nlinarith only [h10, h11, hf, he, hd]

theorem scalable_fixed_two_profile_margin {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ m : ℕ in Filter.atTop,
      ∀ marks : Fin (Spin.Structured.bsched m) →
          (Fin 2 ↪ Fin (Spin.Structured.Lsched m)),
        ∀ rowWeights : Fin 2 → ℕ,
          (∑ A : Fin (Spin.Structured.bsched m) → Finset (Fin 2),
            if Spin.Structured.Placement.rowCount A = rowWeights then
              (Spin.Structured.ConcreteRoute.regionKernel
                (fun j => (A j).map (marks j))).expect
                (fun regions =>
                  ∑ r : State22,
                    inputListKernel
                      (Real.exp (-(scalableThetaTwo /
                        (Spin.Structured.Lsched m : ℝ))))
                      (List.ofFn (scalableFlatInput m regions)) 0 r)
            else 0) ≤
          1000 * (101 / 100 + δ) ^ Spin.Structured.bsched m := by
  have h := scalable_fixed_profile_two_state_margin
    (Q := 2) (θ := scalableThetaTwo) (v := 1 / 1000)
    (c := 101 / 100) (δ := δ)
    (by norm_num [scalableThetaTwo, scalableThetaOne])
    (by norm_num) (by norm_num) hδ (fun _ : Fin 2 => 1)
    (by intro i; norm_num) scalable_fixed_two_row_gate
  filter_upwards [h] with m hm
  intro marks rowWeights
  have hh := hm marks rowWeights
  simp only [Spin.Structured.Placement.profileCoefficient,
    one_pow, Finset.prod_const_one] at hh
  norm_num at hh ⊢
  convert hh using 1
  ring

end Spin.UniformMixer
