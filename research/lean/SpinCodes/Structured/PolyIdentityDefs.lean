/-! Executable integer arithmetic for exact polynomial identities.
No analytic library is needed by the generated certificate modules. -/
namespace Spin.Structured

def polyFactorial : Nat → Nat
  | 0 => 1
  | n + 1 => (n + 1) * polyFactorial n

def polyChoose (n k : Nat) : Nat :=
  if k ≤ n then polyFactorial n / (polyFactorial k * polyFactorial (n - k)) else 0

def polyAdd : List Int → List Int → List Int
  | [], q => q
  | p, [] => p
  | a :: p, b :: q => (a + b) :: polyAdd p q

def polyScale (a : Int) : List Int → List Int
  | [] => []
  | b :: p => (a * b) :: polyScale a p

def polyMul : List Int → List Int → List Int
  | [], _ => []
  | a :: p, q => polyAdd (polyScale a q) (0 :: polyMul p q)

/-- Equality ignores trailing zero coefficients. -/
def polyZero : List Int → Bool
  | [] => true
  | a :: p => a == 0 && polyZero p

def polyEq : List Int → List Int → Bool
  | [], q => polyZero q
  | a :: p, [] => a == 0 && polyZero p
  | a :: p, b :: q => a == b && polyEq p q

structure RatPoly where
  num : List Int
  den : Nat

namespace RatPoly

def add (p q : RatPoly) : RatPoly :=
  let d := Nat.lcm p.den q.den
  ⟨polyAdd (polyScale ((d / p.den : Nat) : Int) p.num) (polyScale ((d / q.den : Nat) : Int) q.num), d⟩

def mul (p q : RatPoly) : RatPoly := ⟨polyMul p.num q.num, p.den * q.den⟩

def constant (a : Int) (d : Nat := 1) : RatPoly := ⟨[a], d⟩

def X : RatPoly := ⟨[0, 1], 1⟩

def pow (p : RatPoly) : Nat → RatPoly
  | 0 => constant 1
  | n + 1 => mul (pow p n) p

def sum : List RatPoly → RatPoly
  | [] => constant 0
  | p :: ps => add p (sum ps)

def checkEq (p q : RatPoly) : Bool :=
  0 < p.den && 0 < q.den && polyEq (polyScale q.den p.num) (polyScale p.den q.num)

def checkAdd (p q r : RatPoly) : Bool :=
  0 < p.den && 0 < q.den && 0 < r.den &&
  polyEq
    (polyAdd (polyScale (q.den * r.den) p.num) (polyScale (p.den * r.den) q.num))
    (polyScale (p.den * q.den) r.num)

def checkMul (p q r : RatPoly) : Bool :=
  0 < p.den && 0 < q.den && 0 < r.den &&
  polyEq (polyScale r.den (polyMul p.num q.num)) (polyScale (p.den * q.den) r.num)

end RatPoly
end Spin.Structured


