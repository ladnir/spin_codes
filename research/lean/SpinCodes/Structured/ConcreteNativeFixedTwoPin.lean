import SpinCodes.Structured.ConcreteNativeFixedTwoLimit

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter

example : Tendsto (fun m => concreteFamily.EZ m 2) atTop (nhds 0) := native_two_EZ_tendsto
example : ∀ᶠ m in atTop, concreteFamily.EZ m 2 ≤ 1000*Real.exp (-(1/20:ℝ)*(bsched m:ℝ)) := by
  filter_upwards [native_two_EZ_bound,fixedTwoBound_eventually] with m hm hb
  exact hm.trans hb

#print axioms Placement.actual_two_fair_moment
#print axioms Placement.actual_two_fair_probability_regionMajor
#print axioms native_two_tuple_bound
#print axioms native_two_EZ_bound
#print axioms fixed_two_log_gap
#print axioms fixedTwoBound_eventually
#print axioms native_two_EZ_tendsto
end Spin.Structured.ConcreteNativeFamily
