import SpinCodes.Numeric.FixedDefs

/-! Materialized 11-by-11 fixed-interval arithmetic for finite transfer powers.
Every output entry is computed once. Squaring retains all boundary coordinates;
a positive rational normalization can raise or lower the matrix scale. -/
namespace Spin.Turbo
open Spin.Numeric

def matrixFixSum22 : List Fix → Fix
  | [] => Fix.zero
  | x :: xs => Fix.add x (matrixFixSum22 xs)

def matrixFixEntry22 (a : List (List Fix)) (i j : Nat) : Fix :=
  (a.getD i []).getD j Fix.zero

def matrixFixMul22 (a b : List (List Fix)) : List (List Fix) :=
  List.ofFn (fun i : Fin 11 => List.ofFn (fun j : Fin 11 =>
    matrixFixSum22 (List.ofFn (fun k : Fin 11 =>
      Fix.mul (matrixFixEntry22 a i.val k.val) (matrixFixEntry22 b k.val j.val)))))

def matrixFixScale22 (num den : Int) (a : List (List Fix)) : List (List Fix) :=
  let s := Fix.ofFrac num den
  List.ofFn (fun i : Fin 11 => List.ofFn (fun j : Fin 11 =>
    Fix.mul s (matrixFixEntry22 a i.val j.val)))

def matrixFixSquareStep22 (num den : Int) (a : List (List Fix)) : List (List Fix) :=
  matrixFixScale22 num den (matrixFixMul22 a a)

def matrixFixSquareChain22 (a : List (List Fix)) : List (Int × Int) → List (List Fix)
  | [] => a
  | (num,den) :: ss => matrixFixSquareChain22 (matrixFixSquareStep22 num den a) ss

def matrixFixTerminal22 (a : List (List Fix)) : Fix :=
  matrixFixSum22 (List.ofFn (fun j : Fin 11 => matrixFixEntry22 a 0 j.val))
end Spin.Turbo


namespace Spin.Turbo
open Spin.Numeric

/-- Compare one computed stage against an independently materialized enclosure.
Certificate modules check one stage at a time; never reduce an entire chain. -/
def matrixFixEnclosesCheck22 (outer inner : List (List Fix)) : Bool :=
  (List.range 11).all (fun i => (List.range 11).all (fun j =>
    let o := matrixFixEntry22 outer i j
    let v := matrixFixEntry22 inner i j
    decide (o.lo ≤ v.lo ∧ v.hi ≤ o.hi)))

def matrixFixSquareStageCheck22 (previous proposed : List (List Fix)) (num den : Int) : Bool :=
  decide (0 < num ∧ 0 < den) &&
    matrixFixEnclosesCheck22 proposed (matrixFixSquareStep22 num den previous)
end Spin.Turbo
