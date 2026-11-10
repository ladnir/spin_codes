import SpinCodes.Turbo.SparseFixedDefs22

/-! Core-only integer description of physical rows j ≤ 2 and the check that
a regional base stage dominates choose(16,d)·B_d at scale den·mult. -/
namespace Spin.Turbo

/-- Entry (i,k) of `profileOutwardRows22` with integer numerators. -/
def profileRowInt22 (birth first returns : List Int) (i k : Nat) : Int :=
  if i = 0 then birth.getD k 0 else if k = 0 then returns.getD i 0
  else if k = 2 then first.getD i 0 else 0

def sparseBinom16Small22 (d : Nat) : Int :=
  match d with
  | 0 => 1
  | 1 => 16
  | 2 => 120
  | _ => 0

def sparseRegionalBaseCheck22 (row : Nat → Nat → Nat → Int) (mult : Int)
    (base : Array Int) : Bool :=
  sparseAll3Bool22 fun d => sparseAll11Bool22 fun i => sparseAll11Bool22 fun k =>
    decide (sparseBinom16Small22 d * row d i k * mult ≤ sparseIntPolyMatrixGet22 base d i k)

end Spin.Turbo
