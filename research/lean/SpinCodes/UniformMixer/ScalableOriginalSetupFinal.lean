import SpinCodes.UniformMixer.ScalableDenseConditionalFinal

/-! Original sampled Scalable setup closure. The conditioning and first-moment
argument reuses `Setup.prob_bad_le_cond_sum` from the old Structured proof;
the encoder weight here is the revised one-bit route and sampled GL(22,2) chain. -/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.UniformMixer
open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute Finset Filter

def scalableMomentSetup (m : ℕ) :
    Setup (NativeMessage m) (NativeSeed m × NativeMessage m)
      (NativeSeed m) (ScalableRouteSeed m × ScalableStateMaps m) where
  Pout := nativeSeedLaw m
  Pin := (seedLaw (Lsched m) (bsched m)).prod (scalableStateMapLaw m)
  outer := fun out x => (out, x)
  innerWt := fun ω u => sampledOutputWeight
    (scalableFlatInput m (routeEval
      (nativeRows m u.1 ((nativeMessageEquiv m).symm u.2)) ω.1))
    (scalableFlatMaps m ω.2) 0

def scalableBad (m : ℕ) (ω : ScalableSetupSeed m) : Prop :=
  ∃ x : NativeMessage m, x ≠ 0 ∧
    sampledOutputWeight (scalableFlatInput m
      (routeEval (nativeRows m ω.1 ((nativeMessageEquiv m).symm x)) ω.2.1))
      (scalableFlatMaps m ω.2.2) 0 ≤ threshold m

theorem scalableBad_iff_count (m : ℕ) (ω : ScalableSetupSeed m) :
    scalableBad m ω ↔ 1 ≤ (scalableMomentSetup m).Z (threshold m) ω := by
  simp only [scalableBad, Setup.Z, Nat.one_le_iff_ne_zero, card_ne_zero,
    Finset.nonempty_def, nonzeroMsgs, mem_filter, mem_univ, true_and]
  rfl

theorem scalable_selected_moment_eq_layer (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m)) (hQ : 0 < Q) :
    (((nativeSeedLaw m).condition (nativeGood m) hpos).prod
      (scalableMomentSetup m).Pin).expect
      (fun ω => ((scalableMomentSetup m).ZQ (threshold m)
        (fun _ x => occupation x) Q ω : ℝ)) =
      scalableSelectedLayer m Q hpos := by
  rw [Setup.expect_ZQ_eq_outerLaw]
  unfold scalableSelectedLayer
  congr 1
  funext out
  apply sum_congr
  · ext x
    simp only [nonzeroMsgs, mem_filter, mem_univ, true_and]
    constructor
    · exact fun h => h.2
    · intro hx
      refine ⟨?_, hx⟩
      intro hz
      have : occupation x = 0 := by simp [hz, occupation]
      omega
  · intro x hx
    rfl

theorem scalable_failure_le_bad_outer_add_total (m : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m)) :
    (scalableSetupLaw m).prob (scalableBad m) ≤
      (nativeSeedLaw m).prob (fun out => ¬ nativeGood m out) +
        ∑ Q ∈ Ico 1 (Lsched m + 1), scalableGoodLayer m Q := by
  have h := (scalableMomentSetup m).prob_bad_le_cond_sum (threshold m)
    (fun _ x => occupation x) (nativeGood m) hpos (Lsched m) (by
      intro out x hx
      exact mem_Ico.mpr ⟨occupation_pos x (mem_filter.mp hx).2,
        Nat.lt_succ_of_le (occupation_le x)⟩)
  have he : (scalableSetupLaw m).prob (scalableBad m) =
      (scalableMomentSetup m).joint.prob
        (fun ω => 1 ≤ (scalableMomentSetup m).Z (threshold m) ω) := by
    exact FinPMF.prob_congr _ (scalableBad_iff_count m)
  rw [he]
  convert h using 1
  congr 1
  apply sum_congr rfl
  intro Q hQ
  change scalableGoodLayer m Q = _
  simpa only [scalableGoodLayer, dite_eq_left hpos, scalableMomentSetup] using
    (scalable_selected_moment_eq_layer m Q hpos (mem_Ico.mp hQ).1).symm

theorem scalable_failure_tendsto_of_dense_rates {η C : ℝ}
    (h : ScalableDenseRates η C) (hη : 0 < η) (hC : 0 < C) :
    Tendsto (fun m => (scalableSetupLaw m).prob (scalableBad m))
      atTop (nhds 0) := by
  have ht := (raw_selection_failure_tendsto nativeTail_tendsto).add
    (scalableGoodLayer_total_tendsto_of_dense_rates h hη hC)
  simp only [add_zero] at ht
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds ht
  · exact Eventually.of_forall (fun m => FinPMF.prob_nonneg _ _)
  · filter_upwards [native_good_positive_eventually] with m hm
    exact scalable_failure_le_bad_outer_add_total m hm

#print axioms scalable_failure_le_bad_outer_add_total
#print axioms scalable_failure_tendsto_of_dense_rates
end Spin.UniformMixer
