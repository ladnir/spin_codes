import SpinCodes.Turbo.SparseListDefs22

/-! Core-only definitions of entrywise `x ^ (a/b)` integer certificates
(generalizing `sparseRpow45CheckL22`, which is the case `a = 4`, `b = 5`).
Soundness: `SparseRpowABCert22`. -/
namespace Spin.Turbo

/-- One entry: `0 ≤ y` and `x^a · sy^b ≤ y^b · sx^a`. -/
def sparseRpowABEntryL22 (a b : Nat) (sx sy : Int) (x y : Int) : Bool :=
  decide (0 ≤ y) && decide (x ^ a * sy ^ b ≤ y ^ b * sx ^ a)

def sparseRpowABRowL22 (a b : Nat) (sx sy : Int) (xs ys : List Int) : Bool :=
  (List.zipWith (sparseRpowABEntryL22 a b sx sy) xs ys).all id

def sparseRpowABCheckL22 (a b : Nat) (sx sy : Int) (X Y : List (List Int)) : Bool :=
  decide (0 < sx) && decide (0 < sy) && sparseMatShapeL22 X && sparseMatShapeL22 Y &&
    (List.zipWith (sparseRpowABRowL22 a b sx sy) X Y).all id

end Spin.Turbo
