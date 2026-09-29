import SpinCodes.Structured.ConcreteNativeFixedAll
import SpinCodes.Structured.DenseOccupationMixedCertified
import SpinCodes.Structured.ConcreteNativeDenseRates

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter ConcreteOuter

/-- The actual native Structured SPIN construction has minimum distance greater
than 0.11 of its block length with probability tending to one. -/
theorem native_minimum_distance_failure_tendsto :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
        (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11:ℝ)*Nsched m⌋₊))
      atTop (nhds 0) := by
  exact minimum_distance_of_fixed_and_dense_rates native_fixed_range_tendsto
    (by norm_num : (0:ℝ)<4/10000000)
    (by norm_num : (0:ℝ)<24000000000000) DenseGeometry.certified_denseRates

/-- The same conclusion through the established concrete first-moment family. -/
theorem native_probBad_tendsto : Tendsto concreteFamily.probBad atTop (nhds 0) := by
  apply native_minimum_distance_failure_tendsto.congr
  intro m
  rw [concrete_probBad_minimumDistance,threshold_floor]

#print axioms native_minimum_distance_failure_tendsto
#print axioms native_probBad_tendsto
end Spin.Structured.ConcreteNativeFamily
