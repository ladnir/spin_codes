import SpinCodes.UniformMixer.Map22Quadratic
import SpinCodes.UniformMixer.Map22SetParity

/-!
# Six-bit indices as subsets of six coordinates

This equivalence connects the integer submask predicate used by the
production row ANF checker with the Boolean-lattice interval theorem.
-/

namespace Spin.UniformMixer

def bitPositions (j : Fin 64) : Finset (Fin 6) :=
  Finset.univ.filter (fun b => j.val.testBit b.val)

set_option maxRecDepth 8192 in
theorem bitPositions_injective : Function.Injective bitPositions := by
  decide

theorem bitPositions_bijective : Function.Bijective bitPositions := by
  apply (Fintype.bijective_iff_injective_and_card bitPositions).mpr
  constructor
  · exact bitPositions_injective
  · norm_num [Fintype.card_finset]

noncomputable def bitPositionsEquiv : Fin 64 ≃ Finset (Fin 6) :=
  Equiv.ofBijective bitPositions bitPositions_bijective

abbrev isSubmask (a b : Fin 64) : Prop :=
  Nat.land a.val b.val = a.val

set_option maxRecDepth 8192 in
theorem isSubmask_iff_subset :
    ∀ a b : Fin 64,
      isSubmask a b ↔ bitPositions a ⊆ bitPositions b := by
  decide

end Spin.UniformMixer
