import SpinCodes.Structured.ConcreteOuterNative

namespace Spin.Structured.ConcreteOuter

example {k : ℕ} (seed : Seed k) (x : LocalMessage k) :
    encode seed x = accF ((accF (outerWord Golay.enc x ∘ seed.1)) ∘ seed.2) := rfl

example {L k : ℕ} (seed : Seed k) (messages : Fin L → LocalMessage k) (i : Fin L) :
    rowSupports seed messages i = support (encode seed (messages i)) := rfl

example (m : ℕ) (seed : NativeSeed m) : Function.Injective
    (nativeRows m seed : (Fin (Nsched m / 2) → Bool) →
      (Fin (Lsched m) → Finset (Fin (bsched m)))) := native_rows_injective m seed

example (m : ℕ) (seed : NativeSeed m) (message : Fin (Nsched m / 2) → Bool)
    (hne : message ≠ 0) :
    ConcreteRoute.activeRows (nativeRows m seed message) ∈ Finset.Ico 1 (Lsched m + 1) :=
  native_occupation_mem m seed message hne

#print axioms Spin.Structured.ConcreteOuter.accF_injective
#print axioms Spin.Structured.ConcreteOuter.golay_systematic_checked
#print axioms Spin.Structured.ConcreteOuter.golay_injective
#print axioms Spin.Structured.ConcreteOuter.encode_injective
#print axioms Spin.Structured.ConcreteOuter.encode_zero
#print axioms Spin.Structured.ConcreteOuter.activeRows_eq
#print axioms Spin.Structured.ConcreteOuter.expected_spectrum
#print axioms Spin.Structured.ConcreteOuter.native_width
#print axioms Spin.Structured.ConcreteOuter.native_dimension
#print axioms Spin.Structured.ConcreteOuter.native_rows_injective
#print axioms Spin.Structured.ConcreteOuter.native_occupation_mem
#print axioms Spin.Structured.ConcreteOuter.native_totalWeight
#print axioms Spin.Structured.ConcreteOuter.native_round_divisibility

end Spin.Structured.ConcreteOuter
