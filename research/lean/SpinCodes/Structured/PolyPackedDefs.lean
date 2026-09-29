import SpinCodes.Structured.PolyIdentityDefs

/-! Exact checks for sums of polynomial products using integer encoding.
The radix exceeds a bound on the sum of absolute coefficients of the
difference. Soundness is proved in `PolyPacked.lean`. -/

namespace Spin.Structured

def coefficientMass : List Int → Nat
  | [] => 0
  | a :: p => a.natAbs + coefficientMass p

def encodePoly : List Int → Int → Int
  | [], _ => 0
  | a :: p, b => a + b * encodePoly p b

structure PackedTerm where
  left : RatPoly
  right : RatPoly
  factor : Nat

namespace PackedTerm

def coefficients (t : PackedTerm) : List Int :=
  polyScale t.factor (polyMul t.left.num t.right.num)

def massBound (t : PackedTerm) : Nat :=
  t.factor * coefficientMass t.left.num * coefficientMass t.right.num

def valid (D : Nat) (t : PackedTerm) : Bool :=
  0 < t.left.den && 0 < t.right.den && t.left.den * t.right.den * t.factor == D

def encoded (t : PackedTerm) (b : Int) : Int :=
  t.factor * encodePoly t.left.num b * encodePoly t.right.num b

end PackedTerm

def packedCoefficients : List PackedTerm → List Int
  | [] => []
  | t :: ts => polyAdd t.coefficients (packedCoefficients ts)

def packedMass : List PackedTerm → Nat
  | [] => 0
  | t :: ts => t.massBound + packedMass ts

def packedEncoded : List PackedTerm → Int → Int
  | [], _ => 0
  | t :: ts, b => t.encoded b + packedEncoded ts b

def packedBase (ts : List PackedTerm) (r : RatPoly) (rf : Nat) : Nat :=
  packedMass ts + rf * coefficientMass r.num + 1

def checkSumProducts (D : Nat) (ts : List PackedTerm) (r : RatPoly) (rf : Nat) : Bool :=
  0 < D && 0 < r.den && r.den * rf == D && ts.all (PackedTerm.valid D) &&
  packedEncoded ts (packedBase ts r rf) == (rf : Int) * encodePoly r.num (packedBase ts r rf)

end Spin.Structured
