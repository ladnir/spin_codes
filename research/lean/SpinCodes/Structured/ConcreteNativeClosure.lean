import SpinCodes.Structured.EventualRegimes
import SpinCodes.Structured.ConcreteOuterDefectEnvelope

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter

/-- The remaining closure frontier, with the outer spectrum and sparse rate discharged.
This is still conditional on the outer tail and the fixed/dense first moments. -/
theorem distance_of_remaining_regimes
    (htail : Tendsto nativeTail atTop (nhds 0))
    (hfixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0))
    (η : ℝ) (hη : 0 < η) (ε : ℕ → ℝ) (hε : Tendsto ε atTop (nhds 0))
    (hdense : ∀ᶠ m in atTop,
      ∑ Q ∈ Ico (nativeCut m+1) (Lsched m+1), concreteFamily.EZ m Q ≤
        (Lsched m : ℝ)*Real.exp (-η*Nsched m+ε m*Nsched m)) :
    Tendsto concreteFamily.probBad atTop (nhds 0) :=
  concreteFamily.distance_whp_native_eventually (fun _ => rfl) (fun _ => rfl)
    (concrete_selection_failure_tendsto htail) hfixed nativeCut nativeCut_spec
    (concrete_sparse_eventually_of_tail htail) η hη ε hε hdense

end Spin.Structured.ConcreteNativeFamily
