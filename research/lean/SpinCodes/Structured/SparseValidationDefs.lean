import SpinCodes.Structured.SparsePolynomialDefs

namespace Spin.Structured.SparsePolynomial

/-- Range checks needed to turn a selected entry into an upper bound for a
minimum. The kernel bound also identifies natural and integer subtraction. -/
def checkData (j : Nat) (d : WeightData) : Bool :=
  d.kernel ≤ polyChoose 128 j &&
  d.choices.momentMax < 5 && d.choices.freshD < 2 &&
  d.choices.cancelD < 3 + (if d.lowPatterns.isSome then 1 else 0) &&
  d.choices.cancelS.length == 5 && d.choices.freshS.length == 5 &&
  (List.range 5).all (fun i =>
    d.choices.cancelS.getD i 0 < 2 + (if d.lowPatterns.isSome then 1 else 0) &&
    d.choices.freshS.getD i 0 < 2) &&
  match d.lowPatterns with
  | none => true
  | some ps => ps.isEmpty || d.choices.lowMax < ps.length

end Spin.Structured.SparsePolynomial
