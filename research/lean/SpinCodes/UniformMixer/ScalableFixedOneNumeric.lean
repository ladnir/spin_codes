import SpinCodes.UniformMixer.ScalableFixedTwoStateWeight
import SpinCodes.Structured.ConcreteFixedOneBound

/-! A rational two-state margin for one fixed active outer row. -/

noncomputable section

namespace Spin.UniformMixer

def scalableThetaOne : ℝ := 4194303 / 2097152

theorem regionGamma_thetaOne : regionGamma scalableThetaOne = 1 := by
  norm_num [regionGamma, production_liveWeightAverage, scalableThetaOne]

theorem scalable_fixed_one_row_gate :
    scalableTwoRowNormLe (1 / 1000) (101 / 100)
      (scalableFugacityTwo scalableThetaOne
        (fun _ : Fin 1 => 1)) := by
  have hγ : regionGamma scalableThetaOne ≠ 0 := by
    rw [regionGamma_thetaOne]
    norm_num
  rw [scalableFugacityTwo_one, scalableContinuumTwo_zero,
    scalableContinuumTwo_one hγ, regionGamma_thetaOne]
  have he0 : 0 ≤ Real.exp (-1 : ℝ) := (Real.exp_pos _).le
  have he : Real.exp (-1 : ℝ) ≤ 1 / 2 :=
    Spin.Structured.Placement.exp_neg_le_half_of_one_le
      (by norm_num : (1 : ℝ) ≤ 1)
  unfold scalableTwoRowNormLe
  simp only [Matrix.add_apply, one_smul]
  dsimp only [timeRegionGap, coarseD, scalableOneFormula]
  simp only [Bool.false_eq_true, Bool.true_eq_false,
    ite_true, ite_false, zero_add, add_zero, mul_one]
  constructor
  · nlinarith only [he0]
  · change (1 / (liveCount State22 : ℝ)) *
        ((1 - Real.exp (-1 : ℝ)) / 1) +
        (Real.exp (-1 : ℝ) +
          (1 - 1 / (liveCount State22 : ℝ)) * Real.exp (-1 : ℝ)) *
          (1 / 1000) ≤ (101 / 100) * (1 / 1000)
    rw [state22_liveCount]
    norm_num only [Nat.reducePow, Nat.reduceSub, Nat.cast_ofNat,
      div_one]
    nlinarith only [he0, he]

theorem scalable_fixed_one_profile_margin {δ : ℝ} (hδ : 0 < δ) :
    ∀ᶠ m : ℕ in Filter.atTop,
      ∀ marks : Fin (Spin.Structured.bsched m) →
          (Fin 1 ↪ Fin (Spin.Structured.Lsched m)),
        ∀ rowWeights : Fin 1 → ℕ,
          (∑ A : Fin (Spin.Structured.bsched m) → Finset (Fin 1),
            if Spin.Structured.Placement.rowCount A = rowWeights then
              (Spin.Structured.ConcreteRoute.regionKernel
                (fun j => (A j).map (marks j))).expect
                (fun regions =>
                  ∑ r : State22,
                    inputListKernel
                      (Real.exp (-(scalableThetaOne /
                        (Spin.Structured.Lsched m : ℝ))))
                      (List.ofFn (scalableFlatInput m regions)) 0 r)
            else 0) ≤
          1000 * (101 / 100 + δ) ^ Spin.Structured.bsched m := by
  have h := scalable_fixed_profile_two_state_margin
    (Q := 1) (θ := scalableThetaOne) (v := 1 / 1000)
    (c := 101 / 100) (δ := δ)
    (by norm_num [scalableThetaOne]) (by norm_num)
    (by norm_num) hδ (fun _ : Fin 1 => 1)
    (by intro i; norm_num) scalable_fixed_one_row_gate
  filter_upwards [h] with m hm
  intro marks rowWeights
  have hh := hm marks rowWeights
  simp only [Spin.Structured.Placement.profileCoefficient,
    one_pow, Finset.prod_const_one] at hh
  norm_num at hh ⊢
  convert hh using 1
  ring

end Spin.UniformMixer
