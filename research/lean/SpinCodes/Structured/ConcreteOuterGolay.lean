import SpinCodes.Structured.ConcreteOuterAlgebra

set_option maxRecDepth 1000000
set_option maxHeartbeats 0

namespace Spin.Structured.ConcreteOuter

/-- The checked generator is systematic on every 12-bit input. -/
theorem golay_systematic_checked :
    (allWords 12).all (fun m => decide ((Golay.encList m).take 12 = m)) = true := by
  decide +kernel

theorem golay_systematic_list (m : Fin 12 → Bool) :
    (Golay.encList (List.ofFn m)).take 12 = List.ofFn m := by
  have h := List.all_eq_true.mp golay_systematic_checked _ (ofFn_mem_allWords m)
  exact of_decide_eq_true h

theorem golay_injective : Function.Injective Golay.enc := by
  intro x y h
  have hh := congrArg (fun z => (List.ofFn z).take 12) h
  apply List.ofFn_injective
  simpa only [Golay.ofFn_enc, golay_systematic_list] using hh

end Spin.Structured.ConcreteOuter

