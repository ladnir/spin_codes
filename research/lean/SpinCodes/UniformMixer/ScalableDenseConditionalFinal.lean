import SpinCodes.UniformMixer.ScalableDenseSelected
import SpinCodes.UniformMixer.ScalableSparseFinal
import SpinCodes.UniformMixer.ScalableFixedLargeFinal
import SpinCodes.Structured.Regimes

/-! The complete selected-setup first moment follows from the dense rate cover. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteNativeFamily Finset Filter

theorem scalableGoodLayer_total_tendsto_of_dense_rates
    {η C : ℝ} (h : ScalableDenseRates η C)
    (hη : 0 < η) (hC : 0 < C) :
    Tendsto (fun m : ℕ =>
      ∑ Q ∈ Ico 1 (Lsched m + 1), scalableGoodLayer m Q)
      atTop (nhds 0) := by
  apply total_regime_sum scalableGoodLayer Lsched nativeCut
    nativeCut_spec
  · apply fixed_regime
    intro Q hQ
    exact scalableGoodLayer_fixed_tendsto (by
      have := (mem_Ico.mp hQ).1
      omega)
  · exact scalableGoodLayer_sparse_sum_tendsto
  · exact h.good_dense_sum_tendsto hη hC

#print axioms scalableGoodLayer_total_tendsto_of_dense_rates

end Spin.UniformMixer
