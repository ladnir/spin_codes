import SpinCodes.UniformMixer.ScalableDenseLayer
import SpinCodes.UniformMixer.ScalableFixedSelected
import SpinCodes.Structured.ConcreteNativeDenseRates

/-! Dense selected-setup first moment, conditional on the numerical rate cover. -/

noncomputable section
set_option maxHeartbeats 1000000
set_option maxRecDepth 8192
attribute [local instance] Classical.propDecidable

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily
  Spin.Structured.DenseOccupationFixed Finset Filter

theorem ScalableDenseRates.selected_layer_bound {η C : ℝ}
    (h : ScalableDenseRates η C) (hC : 0 ≤ C)
    (m Q : ℕ)
    (hpos : 0 < (nativeSeedLaw m).prob (nativeGood m))
    (hQ : Q ∈ Ico (nativeCut m + 1) (Lsched m + 1)) :
    scalableSelectedLayer m Q hpos ≤
      densePrefactor (Lsched m) (bsched m) Q C *
        Real.exp (-η * (Nsched m : ℝ)) := by
  have hQ0 : 0 < Q := by
    have := (mem_Ico.mp hQ).1
    omega
  have hd : (threshold m : ℝ) ≤ (11 / 100) *
      ((Lsched m : ℝ) * bsched m) := by
    simpa only [mul_assoc] using threshold_le m
  rw [scalableSelectedLayer_eq_route]
  calc
    _ ≤ ((nativeSeedLaw m).condition (nativeGood m) hpos).expect
        (fun _ => densePrefactor (Lsched m) (bsched m) Q C *
          Real.exp (-η * (Nsched m : ℝ))) := by
      apply FinPMF.expect_condition_mono
      intro seed hg
      have hh := h.layer_bound hC m Q (threshold m) hQ0 seed hg
        (native_dense_density m Q hQ) hd
      have hrow (x : NativeMessage m) :
          nativeRows m seed ((nativeMessageEquiv m).symm x) =
            fun i => (rowSupports seed x i).map
              (finCongr (native_width m)).toEmbedding := by
        funext i
        simp only [nativeRows, Equiv.apply_symm_apply]
      convert hh using 1
      · apply Finset.sum_congr rfl
        intro x hx
        rw [hrow x]
      · rw [Nsched_eq, Nat.cast_mul]
    _ = _ := FinPMF.expect_const _ _

theorem ScalableDenseRates.good_layer_bound {η C : ℝ}
    (h : ScalableDenseRates η C) (hC : 0 ≤ C)
    (m Q : ℕ)
    (hQ : Q ∈ Ico (nativeCut m + 1) (Lsched m + 1)) :
    scalableGoodLayer m Q ≤
      densePrefactor (Lsched m) (bsched m) Q C *
        Real.exp (-η * (Nsched m : ℝ)) := by
  unfold scalableGoodLayer
  split_ifs with hpos
  · exact h.selected_layer_bound hC m Q hpos hQ
  · have hpre : 0 ≤ densePrefactor (Lsched m) (bsched m) Q C := by
      unfold densePrefactor outerCost
      positivity
    positivity

theorem ScalableDenseRates.good_dense_sum_tendsto {η C : ℝ}
    (h : ScalableDenseRates η C)
    (hη : 0 < η) (hC : 0 < C) :
    Tendsto (fun m : ℕ =>
      ∑ Q ∈ Ico (nativeCut m + 1) (Lsched m + 1),
        scalableGoodLayer m Q) atTop (nhds 0) := by
  apply dense_regime_eventually _ Lsched Nsched η
    (denseRemainder C) hη
  · intro m
    exact sum_nonneg fun Q _ => scalableGoodLayer_nonneg m Q
  · exact Lsched_le_Nsched
  · exact Eventually.of_forall (fun m =>
      native_dense_sum m (nativeCut m) (scalableGoodLayer m)
        hC η (fun Q hQ => h.good_layer_bound hC.le m Q hQ))
  · exact denseRemainder_tendsto C
  · exact Nsched_tendsto

#print axioms ScalableDenseRates.selected_layer_bound
#print axioms ScalableDenseRates.good_dense_sum_tendsto

end Spin.UniformMixer
