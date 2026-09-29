import SpinCodes.Structured.ConcreteNativeClosure
import SpinCodes.Structured.ConcreteOuterTailClosure
import SpinCodes.Structured.ConcreteNativeCodeword

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteEncoder

/-- The actual native distance event now needs only fixed and positive occupation.
The selection event, outer spectrum, and growing sparse regime are discharged. -/
theorem distance_of_fixed_dense
    (hfixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0))
    (η : ℝ) (hη : 0 < η) (ε : ℕ → ℝ) (hε : Tendsto ε atTop (nhds 0))
    (hdense : ∀ᶠ m in atTop,
      ∑ Q ∈ Ico (nativeCut m+1) (Lsched m+1), concreteFamily.EZ m Q ≤
        (Lsched m : ℝ)*Real.exp (-η*Nsched m+ε m*Nsched m)) :
    Tendsto concreteFamily.probBad atTop (nhds 0) :=
  distance_of_remaining_regimes nativeTail_tendsto hfixed η hη ε hε hdense

/-- The same conditional conclusion states the event directly on realized binary codewords. -/
theorem codeword_distance_of_fixed_dense
    (hfixed : ∀ Q ∈ Ico 1 4096, Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0))
    (η : ℝ) (hη : 0 < η) (ε : ℕ → ℝ) (hε : Tendsto ε atTop (nhds 0))
    (hdense : ∀ᶠ m in atTop,
      ∑ Q ∈ Ico (nativeCut m+1) (Lsched m+1), concreteFamily.EZ m Q ≤
        (Lsched m : ℝ)*Real.exp (-η*Nsched m+ε m*Nsched m)) :
    Tendsto (fun m => ((nativeSeedLaw m).prod
      (ConcreteRoutedEncoder.experimentLaw (Lsched m) (bsched m) (rounds m))).prob
        (fun ω => ∃ c ∈ codewords m ω.1 ω.2, c ≠ 0 ∧ binaryWeight c ≤ threshold m))
      atTop (nhds 0) := by
  simpa only [← concrete_probBad_codewords] using distance_of_fixed_dense hfixed η hη ε hε hdense

#print axioms distance_of_fixed_dense
#print axioms codeword_distance_of_fixed_dense

end Spin.Structured.ConcreteNativeFamily
