import SpinCodes.Structured.ConcreteNativeFixedLargeFinal

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter

example {Q : ℕ} (hQ : 3≤Q) : Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0) :=
  native_fixed_large_tendsto hQ

example {Q : ℕ} (hQ : 3≤Q) : ∀ᶠ m in atTop,
    concreteFamily.EZ m Q≤(1600/3)*Real.exp (-(1/2000:ℝ)*(Q:ℝ)*(bsched m:ℝ)) :=
  native_fixed_large_eventually hQ

#print axioms choose_fixed_fugacities
#print axioms fixed_profile_uniform
#print axioms weighted_shell_bound
#print axioms normalizedFixedProfile_bound
#print axioms fixed_profile_cancellation
#print axioms fixed_counted_profile_bound
#print axioms native_fixed_large_profile_bound
#print axioms fixedLarge_counting_cost_le
#print axioms fixedLargeBound_tendsto
#print axioms native_large_EZ_tendsto
#print axioms fixedLargeNorm_proved
#print axioms native_fixed_large_tendsto
#print axioms native_fixed_large_eventually
end Spin.Structured.ConcreteNativeFamily
