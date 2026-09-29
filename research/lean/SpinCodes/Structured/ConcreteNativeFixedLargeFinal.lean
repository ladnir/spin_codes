import SpinCodes.Structured.ConcreteNativeFixedLargeClosure
import SpinCodes.Structured.ConcreteFixedLargeNorm

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Filter ConcreteFixedNumeric Placement

theorem fixedLargeNorm_proved {Q : ℕ} (hQ : 3 ≤ Q) : FixedLargeNorm Q := by
  intro u
  have h := fugacityContinuum_paper_rowNorm hQ (fixedFugacity u) (fun i => (fixedFugacity_pos u i).le)
  simpa only [fixedNormCost,gThree_eq,paperLargeTheta] using h

/-- Unconditional decay of every fixed occupation at least three in the actual native family. -/
theorem native_fixed_large_tendsto {Q : ℕ} (hQ : 3 ≤ Q) :
    Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0) :=
  native_large_EZ_tendsto (by omega) (fixedLargeNorm_proved hQ)

theorem native_fixed_large_eventually {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m in atTop, concreteFamily.EZ m Q≤(1600/3)*Real.exp (-(1/2000:ℝ)*(Q:ℝ)*(bsched m:ℝ)) := by
  filter_upwards [native_large_EZ_bound (by omega) (fixedLargeNorm_proved hQ),fixedLargeBound_eventually Q] with m hm hb
  exact hm.trans hb

#print axioms fixedLargeNorm_proved
#print axioms native_fixed_large_tendsto
#print axioms native_fixed_large_eventually
end Spin.Structured.ConcreteNativeFamily
