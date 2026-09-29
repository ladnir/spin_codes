import SpinCodes.Structured.ConcreteNativeFixedLargeProfile
import SpinCodes.Structured.ConcreteNativeFixedLargeLayer

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter DenseOccupationFixed

theorem native_large_tuple_bound {Q : ℕ} (hQ : 0 < Q) (hK : FixedLargeNorm Q) :
    ∀ᶠ m in atTop, ∀ seed, nativeGood m seed →
      (∑ x∈univ.filter (fun x : NativeMessage m => occupation x=Q),
        failureProbability (tupleWiring m) (rowSupports seed x) (threshold m))≤fixedLargeBound Q m := by
  filter_upwards [native_fixed_large_profile_bound hQ hK] with m hm
  intro seed hg
  have h := sum_occupation_le_profile_bound seed Q
    (fun x => failureProbability (tupleWiring m) (rowSupports seed x) (threshold m))
    (fun S hS w => hm seed hg S hS w)
  simp only [native_width] at h
  exact h.trans (fixedLarge_counting_cost_le Q m)

theorem native_large_EZ_bound {Q : ℕ} (hQ : 0 < Q) (hK : FixedLargeNorm Q) :
    ∀ᶠ m in atTop, concreteFamily.EZ m Q≤fixedLargeBound Q m := by
  filter_upwards [native_good_positive_eventually,native_large_tuple_bound hQ hK] with m hm hb
  exact EZ_le_tuple_layer_bound m Q hm hb

/-- Every fixed positive occupation layer decays once its actual continuum norm is certified. -/
theorem native_large_EZ_tendsto {Q : ℕ} (hQ : 0 < Q) (hK : FixedLargeNorm Q) :
    Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0) := by
  apply tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds (fixedLargeBound_tendsto hQ)
  · exact Eventually.of_forall (fun m => concreteFamily.EZ_nonneg m Q)
  · exact native_large_EZ_bound hQ hK

#print axioms native_large_tuple_bound
#print axioms native_large_EZ_bound
#print axioms native_large_EZ_tendsto
end Spin.Structured.ConcreteNativeFamily
