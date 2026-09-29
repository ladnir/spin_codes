import SpinCodes.Structured.DenseOccupationFixedVertexDefs

namespace Spin.Structured.ConcreteFixedNumeric
open Spin.Numeric DenseOccupationFixed

def rowEval (x s i u v : QInput) : Option Fix := do
  let h ← Fix.fhEnt x.interval 30
  let lu ← Fix.flogQ u.num u.den 30
  let lv ← Fix.flogQ v.num v.den 30
  let lg ← Fix.flogQ 20500 16891 30
  pure (Fix.add (Fix.add (Fix.add (Fix.add
    (Fix.add (Fix.sub (Fix.sub (Fix.add (Fix.mul s.interval x.interval) i.interval) h)
      (Fix.mul x.interval lu)) lv)
    (Fix.add (Fix.ofFrac (-1203) 5125) (Fix.mul (Fix.ofFrac 4 3) lg)))
    (Fix.ofFrac 767031881 3276800000))
    (Fix.mul (Fix.ofFrac 4 39) Fix.log2)) Fix.zero)

def rowCheck (x s i u v : QInput) : Bool :=
  match rowEval x s i u v with
  | none => false
  | some a => decide (a.hi ≤ (Fix.ofFrac (-8679) 10000000).lo)
end Spin.Structured.ConcreteFixedNumeric
