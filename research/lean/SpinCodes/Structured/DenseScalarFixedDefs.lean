import SpinCodes.Numeric.FixedDefs

namespace Spin.Structured.DenseScalarFixed
open Spin.Numeric

def fEval (β z : Fix) (d : Nat) : Fix :=
  let g0 := Fix.add (Fix.sub (Fix.ofInt 1) β) (Fix.mul β z)
  let g1 := Fix.add β (Fix.mul (Fix.sub (Fix.ofInt 1) β) z)
  Fix.mul (Fix.pow g0 (128-d)) (Fix.pow g1 d)

end Spin.Structured.DenseScalarFixed
