import SpinCodes.Structured.ConcreteNativeFixedLargeFinal
import SpinCodes.Structured.ConcreteNativeFixedOneLimit

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter

/-- Every fixed positive occupation has a vanishing actual native first moment. -/
theorem native_fixed_EZ_tendsto (Q : ℕ) (hQ : 0<Q) :
    Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0) := by
  by_cases h1 : Q=1
  · subst Q
    exact native_one_EZ_tendsto
  by_cases h2 : Q=2
  · subst Q
    exact native_two_EZ_tendsto
  exact native_fixed_large_tendsto (by omega)

theorem native_fixed_range_tendsto :
    ∀ Q∈Ico 1 4096, Tendsto (fun m => concreteFamily.EZ m Q) atTop (nhds 0) := by
  intro Q hQ
  exact native_fixed_EZ_tendsto Q (by have := (mem_Ico.mp hQ).1; omega)

#print axioms native_fixed_EZ_tendsto
#print axioms native_fixed_range_tendsto
end Spin.Structured.ConcreteNativeFamily
