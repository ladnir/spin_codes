import SpinCodes.Structured.FiberNumericsBridge
import SpinCodes.Structured.FiberNumericsData.Weight2
import SpinCodes.Structured.FiberNumericsData.Weight64

/-! Representative applications of the integer certificates. The spectrum
premise is explicit and remains to be discharged by the full map replay. -/

noncomputable section
namespace Spin.Structured.FiberNumerics

theorem sample_kernel_two
    (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = Data.spectrum.getD w 0) :
    (ConcreteMaps.kernelLayer 2).card = 0 :=
  kernel_count_of_certificate Data.spectrum hs Data.values2_checked Data.signed2_checked

theorem sample_cap_two
    (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = Data.spectrum.getD w 0)
    (q : Finset (Fin 19)) :
    (ConcreteMaps.syndromeFiber 2 q).card ≤ 61 := by
  apply fiber_cap_of_fourier_certificate Data.spectrum hs Data.values2_checked
  rw [Data.absolute2_checked]
  exact Data.cap2_checked

theorem sample_cap_middle
    (hs : ∀ w : Fin 129, weightCounts ConcreteMaps.CtransposeSet w = Data.spectrum.getD w 0)
    (q : Finset (Fin 19)) :
    (ConcreteMaps.syndromeFiber 64 q).card ≤ Data.caps.getD 64 0 := by
  apply fiber_cap_of_fourier_certificate Data.spectrum hs Data.values64_checked
  rw [Data.absolute64_checked]
  exact Data.cap64_checked

end Spin.Structured.FiberNumerics
