import SpinCodes.UniformMixer.ScalableFixedLargeLayer
import SpinCodes.UniformMixer.ScalableFixedSelected

/-! Fixed occupations of at least three vanish for the actual Scalable setup. -/

noncomputable section
set_option maxHeartbeats 1000000
attribute [local instance] Classical.propDecidable
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute
  Finset Filter

private theorem nativeRows_symm_eq_large (m : ℕ) (seed : NativeSeed m)
    (x : NativeMessage m) :
    nativeRows m seed ((nativeMessageEquiv m).symm x) =
      fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding := by
  funext i
  simp only [nativeRows, Equiv.apply_symm_apply]

theorem scalable_native_large_tuple_eventually {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m in atTop, ∀ seed : NativeSeed m, nativeGood m seed →
      (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
        (law (fun i => (rowSupports seed x i).map
          (finCongr (native_width m)).toEmbedding)).expect
            (scalableRegionFailure m (threshold m))) ≤
        (1600 / 3 : ℝ) * Real.exp
          (-(1 / 2500 : ℝ) * (Q : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [scalable_native_large_tuple_bound hQ,
    fixedLargeBound_eventually Q] with m ht hb
  intro seed hg
  have hc := scalableTiltCorrection_le Q m
  have hnonneg : 0 ≤ fixedLargeBound Q m := by
    unfold fixedLargeBound
    positivity
  calc
    _ ≤ fixedLargeBound Q m * scalableTiltCorrection Q m := ht seed hg
    _ ≤ ((1600 / 3 : ℝ) * Real.exp
          (-(1 / 2000 : ℝ) * (Q : ℝ) * (bsched m : ℝ))) *
        Real.exp ((Q : ℝ) * (bsched m : ℝ) / 2000000) := by
      exact (mul_le_mul_of_nonneg_right hb (Real.exp_pos _).le).trans
        (mul_le_mul_of_nonneg_left hc (by positivity))
    _ = (1600 / 3 : ℝ) * Real.exp
          (-(1 / 2000 : ℝ) * (Q : ℝ) * (bsched m : ℝ) +
            (Q : ℝ) * (bsched m : ℝ) / 2000000) := by
      rw [mul_assoc, ← Real.exp_add]
    _ ≤ _ := by
      apply mul_le_mul_of_nonneg_left _ (by norm_num)
      apply Real.exp_le_exp.mpr
      have hn : 0 ≤ (Q : ℝ) * (bsched m : ℝ) := by positivity
      nlinarith

theorem scalableSelectedLayer_large_bound {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m in atTop, ∀ hpos,
      scalableSelectedLayer m Q hpos ≤
        (1600 / 3 : ℝ) * Real.exp
          (-(1 / 2500 : ℝ) * (Q : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [scalable_native_large_tuple_eventually hQ] with m hm
  intro hpos
  rw [scalableSelectedLayer_eq_route]
  calc
    _ ≤ ((nativeSeedLaw m).condition (nativeGood m) hpos).expect
          (fun _ => (1600 / 3 : ℝ) * Real.exp
            (-(1 / 2500 : ℝ) * (Q : ℝ) * (bsched m : ℝ))) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      simpa only [nativeRows_symm_eq_large] using hm seed hg
    _ = _ := FinPMF.expect_const _ _

theorem scalableGoodLayer_large_eventually {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m in atTop,
      scalableGoodLayer m Q ≤
        (1600 / 3 : ℝ) * Real.exp
          (-(1 / 2500 : ℝ) * (Q : ℝ) * (bsched m : ℝ)) := by
  filter_upwards [native_good_positive_eventually,
    scalableSelectedLayer_large_bound hQ] with m hp hb
  simp only [scalableGoodLayer, dite_eq_left hp]
  exact hb hp

theorem scalableGoodLayer_large_tendsto {Q : ℕ} (hQ : 3 ≤ Q) :
    Tendsto (fun m : ℕ => scalableGoodLayer m Q) atTop (nhds 0) := by
  have hb : Tendsto (fun m => (bsched m : ℝ)) atTop atTop :=
    tendsto_natCast_atTop_atTop.comp bsched_tendsto
  have hneg : -(1 / 2500 : ℝ) * (Q : ℝ) < 0 := by
    have hQr : (0 : ℝ) < Q := by exact_mod_cast (by omega : 0 < Q)
    nlinarith
  have he := (Real.tendsto_exp_atBot.comp
    (hb.const_mul_atTop_of_neg hneg)).const_mul (1600 / 3)
  simp only [mul_zero] at he
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le'
    tendsto_const_nhds he
  · exact Eventually.of_forall (fun m => scalableGoodLayer_nonneg m Q)
  · exact scalableGoodLayer_large_eventually hQ

theorem scalableGoodLayer_fixed_tendsto {Q : ℕ} (hQ : 0 < Q) :
    Tendsto (fun m : ℕ => scalableGoodLayer m Q) atTop (nhds 0) := by
  by_cases h1 : Q = 1
  · subst Q
    exact scalableGoodLayer_one_tendsto
  by_cases h2 : Q = 2
  · subst Q
    exact scalableGoodLayer_two_tendsto
  exact scalableGoodLayer_large_tendsto (by omega)

#print axioms scalableGoodLayer_large_tendsto
#print axioms scalableGoodLayer_fixed_tendsto
end Spin.UniformMixer
