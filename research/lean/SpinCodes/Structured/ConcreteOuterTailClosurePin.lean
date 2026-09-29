import SpinCodes.Structured.ConcreteOuterTailClosure

open Spin Spin.Structured Spin.Structured.ConcreteOuter Spin.Structured.ConcreteNativeFamily
open Filter Finset
attribute [local instance] Classical.propDecidable

example : Tendsto nativeTail atTop (nhds 0) := nativeTail_tendsto
example : ∀ᶠ m in atTop, 0 < (nativeSeedLaw m).prob (nativeGood m) := by
  classical
  exact native_good_positive_eventually
example : Tendsto concreteFamily.probNotGood atTop (nhds 0) := native_selection_failure_tendsto
example : ∀ᶠ m in atTop, ∀ Q∈Ico 4096 (nativeCut m+1),
    concreteFamily.EZ m Q ≤ Real.exp (-(0.006:ℝ)*(Q:ℝ)*bsched m) := native_sparse_eventually

#print axioms denseTail_closed_real
#print axioms dense_near_gap
#print axioms dense_pathEntropy_le
#print axioms denseShell_le
#print axioms outerTail_le
#print axioms outerTail_tendsto
#print axioms nativeTail_tendsto
#print axioms native_good_positive_eventually
#print axioms native_selection_failure_tendsto
#print axioms native_sparse_eventually

