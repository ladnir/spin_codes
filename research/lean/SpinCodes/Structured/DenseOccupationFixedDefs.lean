import SpinCodes.Numeric.FixedDefs
import SpinCodes.Structured.SparseWeights

namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric

def iminFix (a b : Fix) : Fix := ⟨imin a.lo b.lo, imin a.hi b.hi⟩
def imaxFix (a b : Fix) : Fix := ⟨imax a.lo b.lo, imax a.hi b.hi⟩
def isum (xs : List Fix) : Fix := xs.foldr Fix.add Fix.zero

def checkPowers (z : Fix) (xs : List Fix) (n : Nat) : Bool :=
  (xs.getD 0 Fix.zero == Fix.ofInt 1) &&
    (List.range n).all (fun k => xs.getD (k+1) Fix.zero == Fix.mul (xs.getD k Fix.zero) z)

/-- The fallback makes the table total without requiring bounds on emitted weights. -/
def powers (z : Fix) (xs : List Fix) (n k : Nat) : Fix :=
  if k ≤ n then xs.getD k Fix.zero else Fix.pow z k

/-- Normalize each hypergeometric coefficient before summation, avoiding huge intermediate intervals. -/
def moment (w j : Nat) (zp : Nat → Fix) : Fix :=
  isum ((List.range 129).map fun h =>
    let c := polyChoose w h * polyChoose (128-w) (j-h)
    if h ≤ j ∧ c ≠ 0 then
      Fix.mul (Fix.ofFrac c (polyChoose 128 j)) (zp (w+j-2*h))
    else Fix.zero)

def pattern (j : Nat) (pat : List (Nat × Nat)) (zp : Nat → Fix) : Fix :=
  isum (pat.map fun (out,c) => Fix.mul (Fix.ofFrac c (polyChoose 128 j)) (zp out))

def live (j : Nat) (d : SparsePolynomial.WeightData) : Fix :=
  Fix.ofFrac ((polyChoose 128 j:Int)-d.kernel) (polyChoose 128 j)

def maximumMoment (j : Nat) (zp : Nat → Fix) : Fix :=
  (SparsePolynomial.levels.map fun w => moment w j zp).foldr imaxFix Fix.zero

end Spin.Structured.DenseOccupationFixed
