import SpinCodes.Structured.ConcreteFixedNumericCover

noncomputable section
namespace Spin.Structured.ConcreteFixedNumeric
open Spin.Numeric Set

theorem fugacities_range : ∀ u∈fugacities, (9/10:ℚ)≤u ∧ u≤7 := by decide +kernel

theorem gThree_bound {Q : ℕ} (hQ : 3≤Q) {y : ℝ} (hy : y∈Icc 0 1) :
    -(133/125)*y+(1+1/(Q:ℝ))*Real.log (1-y+y/(127/250))≤gThree := by
  rw [gThree_eq]
  exact Placement.paperLargeExponent_bound hQ y hy

example {x : ℝ} (hx : x∈Icc (13/125) (112/125)) :
    ∃ u : ℚ, u∈fugacities ∧ 0<(u:ℝ) ∧
      Spin.Majorant.refined.toFun x-hEnt x-x*Real.log u+
        Real.log (Spin.mv (127/250) (3/1600) (1/524287) u)+Placement.paperLargeExponent+
        (11/100)*(133/125)/(262144/524287)+(4/39)*Real.log 2 ≤ -(8679/10000000) := by
  simpa only [rate,gThree_eq] using rate_uniform hx

#print axioms rowCheck_sound
#print axioms rowRate_interval
#print axioms rate_uniform
#print axioms fugacities_range
#print axioms gThree_bound
end Spin.Structured.ConcreteFixedNumeric
