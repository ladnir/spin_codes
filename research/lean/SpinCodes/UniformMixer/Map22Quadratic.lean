import SpinCodes.UniformMixer.Map22

/-!
# A finite algebraic check on the production expansion

The 22 production row masks are truth tables of Boolean functions on
six input bits. Their algebraic normal forms have degree at most two.
This finite calculation is a first step toward a structural certificate
for the expansion and feedback spectra; it does not by itself establish
either spectrum or the minimum distances.
-/

namespace Spin.UniformMixer

def monomialDegree (j : Fin 64) : ℕ :=
  ((Finset.univ : Finset (Fin 6)).filter
    (fun b => j.val.testBit b.val)).card

/-- Möbius coefficient of a 64-bit Boolean truth table at monomial `j`. -/
def anfCoefficient (mask : ℕ) (j : Fin 64) : ZMod 2 :=
  ∑ k : Fin 64,
    if Nat.land k.val j.val = k.val ∧ mask.testBit k.val then 1 else 0

set_option maxRecDepth 8192 in
theorem rowMask_degree_le_two :
    ∀ i : Fin 22, ∀ j : Fin 64,
      2 < monomialDegree j → anfCoefficient (rowMask i) j = 0 := by
  decide

end Spin.UniformMixer
