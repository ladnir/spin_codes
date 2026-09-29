import SpinCodes.Structured.PolyIdentityDefs
import SpinCodes.Structured.PolyCertDefs

namespace Spin.Structured

def dropInitialZeros : List Int → List Int
  | [] => []
  | a :: p => if a = 0 then dropInitialZeros p else a :: p

namespace RatPoly

def checkNonpos (p : RatPoly) : Bool :=
  0 < p.den && tailBound (dropInitialZeros p.num) ≤ 0

def checkLE (p q : RatPoly) : Bool :=
  0 < p.den && 0 < q.den && checkNonpos (add p (mul (constant (-1)) q))

end RatPoly
end Spin.Structured
