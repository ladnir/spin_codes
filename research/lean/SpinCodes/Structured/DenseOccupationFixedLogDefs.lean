import SpinCodes.Numeric.BAEvalDefs

namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric

/-- One x log(x/y) term, including the exact x=0 endpoint. -/
def klTermEval (un ud vn vd : Int) (n : Nat) : Option Fix :=
  if un=0 then some Fix.zero else do
    let lg ← Fix.flogQ (un*vd) (ud*vn) n
    pure (Fix.mul (Fix.ofFrac un ud) lg)

def klEval (un ud vn vd : Int) (n : Nat) : Option Fix := do
  let a ← klTermEval un ud vn vd n
  let b ← klTermEval (ud-un) ud (vd-vn) vd n
  pure (Fix.add a b)

end Spin.Structured.DenseOccupationFixed
