import SpinCodes.Numeric.FixedDefs

namespace Spin.Turbo
open Spin.Numeric

/-- Clearing positive denominators gives one integer comparison. -/
def integerRpowUpperCheck22 (p q : Nat) (an ad un ud : Int) : Bool :=
  decide (0 < q ∧ 0 ≤ an ∧ 0 < ad ∧ 0 ≤ un ∧ 0 < ud ∧
    an^p * ud^q ≤ un^q * ad^p)

def profileFixRpowUpperCheck22 (box : Fix) (p q : Nat) (num den : Int) : Bool :=
  integerRpowUpperCheck22 p q box.hi Spin.Numeric.scale num den

end Spin.Turbo
