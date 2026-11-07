import SpinCodes.UniformMixer.ScalableSparseFair
import SpinCodes.UniformMixer.ScalableOuterCounting

/-! Exact sparse occupation-layer bound before the numerical exponent check. -/

noncomputable section
set_option maxHeartbeats 1000000
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute
  Spin.Structured.ConcreteMarked Finset

theorem scalable_sparse_tuple_raw (m Q d : ℕ)
    (seed : NativeSeed m) {B α : ℝ}
    (hB : ∀ w ≤ bsched m,
      (spectrum seed w : ℝ) ≤
        B * ((bsched m).choose w : ℝ))
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (∑ x ∈ univ.filter (fun x : NativeMessage m => occupation x = Q),
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d)) ≤
      ((Lsched m).choose Q : ℝ) *
        (((2 : ℝ) ^ bsched m * B) ^ Q *
          (((1 / markMass (Lsched m) Q ((8 / 5) * α)) ^
              (bsched m) *
              (1024 * sparseContraction α ^
                (bsched m * scalableStepsPerRegion m))) /
            (sparseTilt α) ^ d)) := by
  have hprob : ∀ S : Finset (Fin (Lsched m)), S.card = Q →
      (regionsLaw S (bsched m)).expect
        (scalableRegionFailure m d) ≤
      (((1 / markMass (Lsched m) Q ((8 / 5) * α)) ^
          (bsched m) *
          (1024 * sparseContraction α ^
            (bsched m * scalableStepsPerRegion m))) /
        (sparseTilt α) ^ d) := by
    intro S hS
    rw [scalableRegionFailure_fair]
    simpa only [hS] using
      scalable_fair_sparse_probability m S hα0 hα1 d
  exact scalable_occupation_failure_le_fair m seed Q d hB hprob

#print axioms scalable_sparse_tuple_raw
end Spin.UniformMixer
