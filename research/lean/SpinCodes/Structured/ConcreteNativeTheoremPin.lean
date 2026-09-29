import SpinCodes.Structured.ConcreteNativeTheorem

noncomputable section
namespace Spin.Structured.ConcreteNativeTheoremPin
open Filter ConcreteOuter ConcreteNativeFamily

/-- A full-statement pin for the unconditional actual minimum-distance limit. -/
theorem actual_minimum_distance :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
        (fun ω => minimumDistance m ω.1 ω.2 ≤ ⌊(0.11:ℝ)*Nsched m⌋₊))
      atTop (nhds 0) := native_minimum_distance_failure_tendsto

/-- Every realization has exact rate one half. -/
theorem actual_rate (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (Module.finrank (ZMod 2) (realizedCode m out inner):ℝ)/(Nsched m:ℝ)=1/2 :=
  realizedCode_rate m out inner

#print actual_minimum_distance
#print axioms actual_minimum_distance
#print actual_rate
#print axioms actual_rate
end Spin.Structured.ConcreteNativeTheoremPin
