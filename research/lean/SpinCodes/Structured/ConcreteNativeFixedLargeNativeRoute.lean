import SpinCodes.Structured.ConcreteNativeFixedLargeRoute

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteRoute ConcreteOuter Placement

theorem native_fixed_profile_routed {Q : ℕ} (hQ : 0 < Q) (hK : FixedLargeNorm Q) :
    ∀ᶠ m in atTop, ∀ u : FixedFugacityChoice Q, ∀ e : Fin Q ↪ Fin (Lsched m),
      ∀ rows : Fin Q → Finset (Fin (nativeBlocks m*24)),
      failureProbability (tupleWiring m) (embedRows e rows) (threshold m) ≤
        (((fixedNormCost u+fixedNormSlack u)^(nativeBlocks m*24)/min 1 (3/1600:ℝ))/
          (profileChoices rows*profileCoefficient (fixedFugacity u) (fun i => (rows i).card)))/
          Real.exp (-(((Q:ℝ)*(133/125)/(262144/524287))/(Lsched m:ℝ)))^(threshold m) := by
  have hlim : Tendsto (fun m : ℕ => m+1) atTop atTop := tendsto_add_atTop_nat 1
  filter_upwards [hlim.eventually (fixed_profile_uniform hQ hK)] with m hm
  intro u e rows
  rw [failureProbability,tupleWiring_eq_regionMajor]
  have hh := hm u (nativeBlocks m*24) (rounds m)
    (by change (nativeBlocks m*24)*Lsched m=rounds m*128
        rw [native_width]; simpa only [Nat.mul_comm] using (rounds_length m).symm)
    e rows (threshold m)
  simp only [Lsched,Nat.cast_mul,Nat.cast_ofNat,Nat.cast_add,Nat.cast_one] at hh ⊢
  convert hh using 1
  apply FinPMF.prob_congr
  intro ω
  rfl

#print axioms native_fixed_profile_routed
end Spin.Structured.ConcreteNativeFamily
