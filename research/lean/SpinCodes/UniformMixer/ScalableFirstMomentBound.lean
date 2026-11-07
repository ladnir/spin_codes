import SpinCodes.UniformMixer.ScalableIidWiring

/-!
# First-moment comparison for a selected Scalable outer

For each nonzero message and each selected shared outer, the actual
route-and-inner moment is bounded by a combinatorial routing factor
times the iid 64-bit transfer-matrix moment. The remaining numerical
work must bound this matrix power and sum over outer words.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute

theorem scalable_selected_outer_moment_bound (m : ℕ)
    {z : ℝ} (hz : 0 ≤ z)
    {out : NativeSeed m} (hg : nativeGood m out)
    {message : Fin (Nsched m / 2) → Bool} (hm : message ≠ 0) :
    (law (nativeRows m out message)).expect
      (scalableRoutedMoment m z) ≤
    (((bsched m : ℝ) + 1) ^
        activeRows (nativeRows m out message) *
      ((Lsched m : ℝ) + 1) ^ bsched m) *
    ∑ r : State22,
      ((bernoulliStepKernel
        (density (nativeRows m out message)) z
        (density_nonneg (nativeRows m out message))
        (density_le_one (Lsched_pos m) (bsched_pos m)
          (nativeRows m out message))) ^
        (bsched m * scalableStepsPerRegion m)) 0 r := by
  calc
    _ ≤ (((bsched m : ℝ) + 1) ^
          activeRows (nativeRows m out message) *
        ((Lsched m : ℝ) + 1) ^ bsched m) *
        (iidLaw (Lsched m) (bsched m)
          (density (nativeRows m out message))
          (density_nonneg (nativeRows m out message))
          (density_le_one (Lsched_pos m) (bsched_pos m)
            (nativeRows m out message))).expect
          (scalableRoutedMoment m z) :=
      scalable_good_outer_route_iid_bound m hz hg hm
    _ = _ := by
      rw [scalable_iid_moment_power]

end Spin.UniformMixer
