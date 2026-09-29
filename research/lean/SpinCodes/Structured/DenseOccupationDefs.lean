import SpinCodes.Structured.SparseWeights
import Mathlib.Data.Rat.Cast.Order
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Data.Fintype.Fin

namespace Spin.Structured.DenseOccupation
open SparsePolynomial
open scoped BigOperators

structure QCoords where
  Z : ℚ
  D : ℚ
  S : Fin 5 → ℚ

def live (j : ℕ) (d : WeightData) : ℚ :=
  ((Nat.choose 128 j : ℚ) - d.kernel) / Nat.choose 128 j

def pattern (j : ℕ) (p : List (ℕ × ℕ)) (z : ℚ) : ℚ :=
  (p.map fun (out, c) => (c : ℚ) * z ^ out).sum / Nat.choose 128 j

def hyperMoment (w j : ℕ) (z : ℚ) : ℚ :=
  ((List.range 129).map fun h => if h ≤ j then
    (w.choose h : ℚ) * ((128-w).choose (j-h):ℚ) * z^(w+j-2*h)
    else 0).sum / Nat.choose 128 j

def maximumMoment (j : ℕ) (z : ℚ) : ℚ :=
  (levels.map fun w => hyperMoment w j z).foldr max 0

def cancellationBoundsD (j : ℕ) (d : WeightData) (z : ℚ) : List ℚ :=
  [maximumMoment j z, live j d,
   (d.cap:ℚ)/Nat.choose 128 j * z^minDistance j] ++
  match d.lowPatterns with
  | none => []
  | some ps => [(ps.map fun p => pattern j p z).foldr max 0]

def cancellationBoundsS (j : ℕ) (d : WeightData) (z : ℚ) (i : Fin 5) : List ℚ :=
  let w := levels.getD i 0
  let c := counts.getD i 1
  [hyperMoment w j z,
   (min (Nat.choose 128 j-d.kernel) (c*d.cap):ℕ) /
     ((c:ℚ)*Nat.choose 128 j) * z^distance w j] ++
  match d.lowPatterns with
  | none => []
  | some _ => [pattern j (d.lowShells.getD i []) z/c]

def liveAverage (v : QCoords) : ℚ :=
  (∑ i : Fin 5, (counts.getD i 1:ℚ)*v.S i)/524287

def fixedColumn (j : ℕ) (d : WeightData) (z : ℚ) (v : QCoords) : QCoords where
  Z := (1-live j d)*z^j*v.Z + live j d*z^j*v.D
  D := (((cancellationBoundsD j d z).foldr min (maximumMoment j z))/2 +
      min (maximumMoment j z) (live j d)/(2*524287))*v.Z +
      maximumMoment j z*(v.D+liveAverage v)/2
  S i := let m := hyperMoment (levels.getD i 0) j z
    (((cancellationBoundsS j d z i).foldr min m)/2 + min m (live j d)/(2*524287))*v.Z +
      m*(liveAverage v+if live j d=0 then v.S i else v.D)/2

def probability (j : ℕ) (q : ℚ) : ℚ :=
  (Nat.choose 128 j:ℚ)*q^j*(1-q)^(128-j)

def column (q z : ℚ) (v : QCoords) : QCoords where
  Z := ∑ j : Fin 129, probability j q*(fixedColumn j (Data.weight j) z v).Z
  D := ∑ j : Fin 129, probability j q*(fixedColumn j (Data.weight j) z v).D
  S i := ∑ j : Fin 129, probability j q*(fixedColumn j (Data.weight j) z v).S i

/-- Exact rational arithmetic; successful checks will be replayed by the kernel. -/
def check (q z radius floor : ℚ) (v : QCoords) : Bool :=
  decide (0≤q ∧ q≤1 ∧ 0<z ∧ z≤1 ∧ 0<floor ∧ floor≤v.Z ∧ floor≤v.D ∧
    (∀ i, floor≤v.S i) ∧
    (column q z v).Z≤radius*v.Z ∧ (column q z v).D≤radius*v.D ∧
    (∀ i, (column q z v).S i≤radius*v.S i))

end Spin.Structured.DenseOccupation

