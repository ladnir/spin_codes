/-
Coefficient-list arithmetic for the polynomial sign certificates.

Deliberately dependency-free: the seven degree-257 data modules import only
this file.  Importing Mathlib into each of them made seven parallel module
builds load the whole library at once, which exhausted memory.  Keeping the
`decide` checks over plain `Int` recursion also makes kernel reduction cheap —
no intermediate `List.map` is built.
-/

namespace Spin.Structured

/-- Positive part. -/
def posPart (a : Int) : Int := if 0 ≤ a then a else 0

/-- Sum of the positive parts of a coefficient list. -/
def posSum : List Int → Int
  | [] => 0
  | a :: t => posPart a + posSum t

/-- The verifier's `upper`: constant coefficient plus the positive tail. -/
def tailBound : List Int → Int
  | [] => 0
  | a :: t => a + posSum t

end Spin.Structured
