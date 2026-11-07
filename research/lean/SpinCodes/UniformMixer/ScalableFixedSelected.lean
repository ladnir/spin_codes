import SpinCodes.UniformMixer.ScalableOuterCounting

/-!
One- and two-row first moments for the actual Scalable route and sampled
state maps, conditional on selecting a good shared outer constituent.
-/

noncomputable section
attribute [local instance] Classical.propDecidable

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute Spin.Structured.ConcreteNativeFamily
  Finset Filter

private theorem nativeRows_symm_eq (m : ℕ) (seed : NativeSeed m)
    (x : NativeMessage m) :
    nativeRows m seed ((nativeMessageEquiv m).symm x) =
      fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding := by
  funext i
  simp only [nativeRows, Equiv.apply_symm_apply]

def scalableSelectedLayer (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m)) : ℝ :=
  ((nativeSeedLaw m).condition (nativeGood m) hpos).expect (fun seed =>
    ∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
      ((seedLaw (Lsched m) (bsched m)).prod
        (scalableStateMapLaw m)).prob (fun ω =>
          sampledOutputWeight
            (scalableFlatInput m
              (routeEval (nativeRows m seed
                ((nativeMessageEquiv m).symm x)) ω.1))
            (scalableFlatMaps m ω.2) 0 ≤ threshold m))

theorem scalableSelectedLayer_eq_route (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m)) :
    scalableSelectedLayer m Q hpos =
      ((nativeSeedLaw m).condition (nativeGood m) hpos).expect (fun seed =>
        ∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
          (law (nativeRows m seed ((nativeMessageEquiv m).symm x))).expect
            (scalableRegionFailure m (threshold m))) := by
  unfold scalableSelectedLayer
  congr 1
  funext seed
  apply sum_congr rfl
  intro x hx
  rw [FinPMF.prob_prod_eq_expect _ _ (fun route maps =>
    sampledOutputWeight
      (scalableFlatInput m
        (routeEval (nativeRows m seed ((nativeMessageEquiv m).symm x)) route))
      (scalableFlatMaps m maps) 0 ≤ threshold m)]
  change (seedLaw (Lsched m) (bsched m)).expect (fun route =>
      scalableRegionFailure m (threshold m)
        (routeEval (nativeRows m seed ((nativeMessageEquiv m).symm x)) route)) = _
  rw [← FinPMF.expect_map, permutation_law_eq]

theorem scalableSelectedLayer_nonneg (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m)) :
    0 ≤ scalableSelectedLayer m Q hpos := by
  unfold scalableSelectedLayer
  apply FinPMF.expect_nonneg
  intro seed
  exact sum_nonneg fun x hx => FinPMF.prob_nonneg _ _

theorem scalableSelectedLayer_one_bound :
    ∀ᶠ m : ℕ in atTop, ∀ hpos,
      scalableSelectedLayer m 1 hpos ≤ scalableFixedOneBound m := by
  filter_upwards [scalable_native_one_tuple_bound] with m hm
  intro hpos
  rw [scalableSelectedLayer_eq_route]
  calc
    _ ≤ ((nativeSeedLaw m).condition (nativeGood m) hpos).expect
        (fun _ => scalableFixedOneBound m) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      have hh := hm seed hg
      simp only [nativeRows_symm_eq]
      exact hh
    _ = _ := FinPMF.expect_const _ _

theorem scalableSelectedLayer_two_bound :
    ∀ᶠ m : ℕ in atTop, ∀ hpos,
      scalableSelectedLayer m 2 hpos ≤ scalableFixedTwoBound m := by
  filter_upwards [scalable_native_two_tuple_bound] with m hm
  intro hpos
  rw [scalableSelectedLayer_eq_route]
  calc
    _ ≤ ((nativeSeedLaw m).condition (nativeGood m) hpos).expect
        (fun _ => scalableFixedTwoBound m) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      have hh := hm seed hg
      simp only [nativeRows_symm_eq]
      exact hh
    _ = _ := FinPMF.expect_const _ _

def scalableGoodLayer (m Q : ℕ) : ℝ :=
  if hpos : 0 < (nativeSeedLaw m).prob (nativeGood m) then
    scalableSelectedLayer m Q hpos else 0

theorem scalableGoodLayer_nonneg (m Q : ℕ) :
    0 ≤ scalableGoodLayer m Q := by
  unfold scalableGoodLayer
  split_ifs with hpos
  · exact scalableSelectedLayer_nonneg m Q hpos
  · exact le_refl 0

theorem scalableGoodLayer_one_eventually :
    ∀ᶠ m : ℕ in atTop,
      scalableGoodLayer m 1 ≤
        1000 * Real.exp (-(1 / 50 : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [native_good_positive_eventually,
    scalableSelectedLayer_one_bound,
    scalableFixedOneBound_eventually] with m hp hb he
  simp only [scalableGoodLayer, dite_eq_left hp]
  exact (hb hp).trans he

theorem scalableGoodLayer_two_eventually :
    ∀ᶠ m : ℕ in atTop,
      scalableGoodLayer m 2 ≤
        1000 * Real.exp (-(1 / 20 : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [native_good_positive_eventually,
    scalableSelectedLayer_two_bound,
    scalableFixedTwoBound_eventually] with m hp hb he
  simp only [scalableGoodLayer, dite_eq_left hp]
  exact (hb hp).trans he

theorem scalableGoodLayer_one_tendsto :
    Tendsto (fun m : ℕ => scalableGoodLayer m 1) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have he := (Real.tendsto_exp_atBot.comp
    (hb.const_mul_atTop_of_neg (by norm_num : -(1 / 50 : ℝ) < 0))).const_mul 1000
  simp only [mul_zero] at he
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds he
  · exact Eventually.of_forall (fun m => scalableGoodLayer_nonneg m 1)
  · exact scalableGoodLayer_one_eventually

theorem scalableGoodLayer_two_tendsto :
    Tendsto (fun m : ℕ => scalableGoodLayer m 2) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have he := (Real.tendsto_exp_atBot.comp
    (hb.const_mul_atTop_of_neg (by norm_num : -(1 / 20 : ℝ) < 0))).const_mul 1000
  simp only [mul_zero] at he
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds he
  · exact Eventually.of_forall (fun m => scalableGoodLayer_nonneg m 2)
  · exact scalableGoodLayer_two_eventually

end Spin.UniformMixer
