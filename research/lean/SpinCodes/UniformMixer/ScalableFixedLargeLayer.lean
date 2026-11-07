import SpinCodes.UniformMixer.ScalableFixedLargeProfile
import SpinCodes.Structured.ConcreteNativeFixedLargeLayer

/-! Sum fixed-row profiles over each occupation layer. -/

noncomputable section
set_option maxHeartbeats 1000000
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute
  Spin.Structured.DenseOccupationFixed Finset Filter

theorem scalable_native_large_tuple_bound {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m in atTop, ∀ seed : NativeSeed m,
      nativeGood m seed →
        (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
          (law (fun i => (rowSupports seed x i).map
            (finCongr (native_width m)).toEmbedding)).expect
              (scalableRegionFailure m (threshold m))) ≤
          fixedLargeBound Q m * scalableTiltCorrection Q m := by
  filter_upwards [scalable_fixed_large_profile_bound hQ] with m hm
  intro seed hg
  have hh := sum_occupation_le_profile_bound seed Q
    (fun x => (law (fun i => (rowSupports seed x i).map
      (finCongr (native_width m)).toEmbedding)).expect
        (scalableRegionFailure m (threshold m)))
    (fun S hS w => hm seed hg S hS w)
  simp only [native_width] at hh
  calc
    _ ≤ ((Lsched m).choose Q : ℝ) *
      (((bsched m : ℝ) + 1) ^ Q *
        (fixedProfileCountBound Q (bsched m) (Lsched m)
          (threshold m) * scalableTiltCorrection Q m)) := hh
    _ = (((Lsched m).choose Q : ℝ) *
      (((bsched m : ℝ) + 1) ^ Q *
        fixedProfileCountBound Q (bsched m) (Lsched m)
          (threshold m))) * scalableTiltCorrection Q m := by ring
    _ ≤ _ := mul_le_mul_of_nonneg_right
      (fixedLarge_counting_cost_le Q m) (Real.exp_pos _).le

#print axioms scalable_native_large_tuple_bound
end Spin.UniformMixer
