import SpinCodes.UniformMixer.ScalableFixedFair

/-! Unconditional fair-reference moments for one and two active outer rows. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteMarked
  Spin.Structured.Placement Filter

theorem scalable_fair_one_marks :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin 1 ↪ Fin (Lsched m),
        (regionsLaw (Finset.univ.map marks) (bsched m)).expect
          (scalableRoutedMoment m
            (Real.exp (-(scalableThetaOne / (Lsched m : ℝ))))) ≤
          1000 * (1011 / 2000 : ℝ) ^ bsched m := by
  have h := scalable_fair_region_margin
    (Q := 1) (θ := scalableThetaOne)
    (v := 1 / 1000) (c := 101 / 100) (δ := 1 / 1000)
    (by norm_num [scalableThetaOne]) (by norm_num)
    (by norm_num) (by norm_num) scalable_fixed_one_row_gate
  filter_upwards [h] with m hm
  intro marks
  convert hm marks using 1
  norm_num
  ring

theorem scalable_fair_two_marks :
    ∀ᶠ m : ℕ in atTop,
      ∀ marks : Fin 2 ↪ Fin (Lsched m),
        (regionsLaw (Finset.univ.map marks) (bsched m)).expect
          (scalableRoutedMoment m
            (Real.exp (-(scalableThetaTwo / (Lsched m : ℝ))))) ≤
          1000 * (1011 / 4000 : ℝ) ^ bsched m := by
  have h := scalable_fair_region_margin
    (Q := 2) (θ := scalableThetaTwo)
    (v := 1 / 1000) (c := 101 / 100) (δ := 1 / 1000)
    (by norm_num [scalableThetaTwo, scalableThetaOne])
    (by norm_num) (by norm_num) (by norm_num)
    scalable_fixed_two_row_gate
  filter_upwards [h] with m hm
  intro marks
  convert hm marks using 1
  norm_num
  ring

theorem scalable_fair_one_support :
    ∀ᶠ m : ℕ in atTop,
      ∀ S : Finset (Fin (Lsched m)), S.card = 1 →
        (regionsLaw S (bsched m)).expect
          (scalableRoutedMoment m
            (Real.exp (-(scalableThetaOne / (Lsched m : ℝ))))) ≤
          1000 * (1011 / 2000 : ℝ) ^ bsched m := by
  filter_upwards [scalable_fair_one_marks] with m hm
  intro S hS
  let e : Fin 1 ↪ Fin (Lsched m) := (S.orderEmbOfFin hS).toEmbedding
  have he : Finset.univ.map e = S := by
    rw [Finset.map_eq_image]
    exact Finset.image_orderEmbOfFin_univ S hS
  rw [← he]
  exact hm e

theorem scalable_fair_two_support :
    ∀ᶠ m : ℕ in atTop,
      ∀ S : Finset (Fin (Lsched m)), S.card = 2 →
        (regionsLaw S (bsched m)).expect
          (scalableRoutedMoment m
            (Real.exp (-(scalableThetaTwo / (Lsched m : ℝ))))) ≤
          1000 * (1011 / 4000 : ℝ) ^ bsched m := by
  filter_upwards [scalable_fair_two_marks] with m hm
  intro S hS
  let e : Fin 2 ↪ Fin (Lsched m) := (S.orderEmbOfFin hS).toEmbedding
  have he : Finset.univ.map e = S := by
    rw [Finset.map_eq_image]
    exact Finset.image_orderEmbOfFin_univ S hS
  rw [← he]
  exact hm e

theorem scalable_fair_probability_of_moment (m : ℕ)
    {θ M : ℝ} (hθ : 0 ≤ θ)
    (S : Finset (Fin (Lsched m)))
    (hM : (regionsLaw S (bsched m)).expect
      (scalableRoutedMoment m
        (Real.exp (-(θ / (Lsched m : ℝ))))) ≤ M)
    (d : ℕ) :
    ((regionsLaw S (bsched m)).prod
      (scalableStateMapLaw m)).prob
        (fun ω =>
          sampledOutputWeight (scalableFlatInput m ω.1)
            (scalableFlatMaps m ω.2) 0 ≤ d) ≤
      M / Real.exp (-(θ / (Lsched m : ℝ))) ^ d := by
  apply FinPMF.prob_weight_le_of_moment _ _ d
    (Real.exp_pos _)
    (Real.exp_le_one_iff.mpr (by
      apply neg_nonpos.mpr
      have hL : (0 : ℝ) < Lsched m := by
        exact_mod_cast Lsched_pos m
      exact div_nonneg hθ hL.le))
  rw [FinPMF.expect_prod]
  simp_rw [scalable_fixed_route_moment]
  exact hM

theorem scalable_fair_one_probability :
    ∀ᶠ m : ℕ in atTop,
      ∀ S : Finset (Fin (Lsched m)), S.card = 1 →
        ∀ d : ℕ,
          ((regionsLaw S (bsched m)).prod
            (scalableStateMapLaw m)).prob
            (fun ω =>
              sampledOutputWeight (scalableFlatInput m ω.1)
                (scalableFlatMaps m ω.2) 0 ≤ d) ≤
          (1000 * (1011 / 2000 : ℝ) ^ bsched m) /
            Real.exp (-(scalableThetaOne / (Lsched m : ℝ))) ^ d := by
  filter_upwards [scalable_fair_one_support] with m hm
  intro S hS d
  exact scalable_fair_probability_of_moment m
    (by norm_num [scalableThetaOne]) S (hm S hS) d

theorem scalable_fair_two_probability :
    ∀ᶠ m : ℕ in atTop,
      ∀ S : Finset (Fin (Lsched m)), S.card = 2 →
        ∀ d : ℕ,
          ((regionsLaw S (bsched m)).prod
            (scalableStateMapLaw m)).prob
            (fun ω =>
              sampledOutputWeight (scalableFlatInput m ω.1)
                (scalableFlatMaps m ω.2) 0 ≤ d) ≤
          (1000 * (1011 / 4000 : ℝ) ^ bsched m) /
            Real.exp (-(scalableThetaTwo / (Lsched m : ℝ))) ^ d := by
  filter_upwards [scalable_fair_two_support] with m hm
  intro S hS d
  exact scalable_fair_probability_of_moment m
    (by norm_num [scalableThetaTwo, scalableThetaOne])
    S (hm S hS) d

end Spin.UniformMixer
