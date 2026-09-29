import SpinCodes.Structured.ConcreteNativeDistanceConsequences
import SpinCodes.Structured.ConcreteNativeTheorem

/-! Unconditional success and existence statements for the actual native code. -/
noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter ConcreteOuter
attribute [local instance] Classical.propDecidable

/-- A random actual native realization has relative minimum distance above 11 percent
with probability tending to one. -/
theorem relative_distance_success_tendsto :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
      (fun ω => (11/100:ℝ) < (minimumDistance m ω.1 ω.2 : ℝ)/(Nsched m : ℝ)))
      atTop (nhds 1) :=
  success_probability_tendsto_of_failure native_minimum_distance_failure_tendsto

/-- Every sufficiently large native length admits an actual setup realization
with rate one half and relative minimum distance strictly above 11 percent. -/
theorem eventually_exists_rate_half_distance_gt_eleven_percent :
    ∀ᶠ m in atTop, ∃ out : NativeSeed m, ∃ inner : InnerSeed m,
      (11/100:ℝ) < (minimumDistance m out inner : ℝ)/(Nsched m : ℝ) ∧
      (Module.finrank (ZMod 2) (realizedCode m out inner) : ℝ)/(Nsched m : ℝ) = 1/2 ∧
      0 < ((nativeSeedLaw m).prod
        (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).p (out, inner) :=
  eventually_exists_good_code_of_failure native_minimum_distance_failure_tendsto

end Spin.Structured.ConcreteNativeFamily
