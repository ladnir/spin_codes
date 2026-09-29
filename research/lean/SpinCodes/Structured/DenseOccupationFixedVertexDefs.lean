import SpinCodes.Structured.DenseOccupationFixedLogDefs

namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric

structure QInput where
  num : Int
  den : Int
  deriving Repr

def QInput.interval (a : QInput) : Fix := Fix.ofFrac a.num a.den

def vertexEval (m c p y radius z α x : QInput) (n : Nat) : Option Fix := do
  let a ← klEval α.num α.den p.num p.den n
  let b ← klEval x.num x.den y.num y.den n
  let lr ← Fix.flogQ radius.num radius.den n
  let lz ← Fix.flogQ z.num z.den n
  pure (Fix.sub
    (Fix.add (Fix.add (Fix.add
      (Fix.mul α.interval (Fix.add (Fix.mul m.interval x.interval) c.interval)) a)
      (Fix.mul α.interval b)) (Fix.divInt lr 128))
    (Fix.mul (Fix.ofFrac 11 100) lz))

/-- Successful evaluation with an explicitly checked integer upper endpoint. -/
def vertexCheck (m c p y radius z α x : QInput) (n : Nat) (upper : Int) : Bool :=
  match vertexEval m c p y radius z α x n with
  | none => false
  | some v => decide (v.hi ≤ upper)

end Spin.Structured.DenseOccupationFixed
