import SpinCodes.Structured.DenseOccupationFixedDefs

namespace Spin.Structured.DenseOccupationFixed
open Spin.Numeric SparsePolynomial

structure FCoords where
  Z : Fix
  D : Fix
  S : Fin 5 → Fix

def boundsD (j : Nat) (d : WeightData) (zp : Nat → Fix) : List Fix :=
  [maximumMoment j zp, live j d,
   Fix.mul (Fix.ofFrac d.cap (polyChoose 128 j)) (zp (minDistance j))] ++
  match d.lowPatterns with
  | none => []
  | some ps => [(ps.map fun pat => pattern j pat zp).foldr imaxFix Fix.zero]

def boundsS (j : Nat) (d : WeightData) (zp : Nat → Fix) (i : Fin 5) : List Fix :=
  let w := levels.getD i 0
  let c := counts.getD i 1
  [moment w j zp,
   Fix.mul (Fix.ofFrac ((min (polyChoose 128 j-d.kernel) (c*d.cap):Nat):Int) ((c*polyChoose 128 j:Nat):Int))
      (zp (distance w j))] ++
  match d.lowPatterns with
  | none => []
  | some _ => [Fix.divInt (pattern j (d.lowShells.getD i []) zp) c]

def liveAverage (v : FCoords) : Fix :=
  Fix.divInt (isum (List.ofFn fun i : Fin 5 => Fix.mul (Fix.ofInt (counts.getD i 1)) (v.S i))) 524287

def fixedColumn (j : Nat) (d : WeightData) (zp : Nat → Fix) (v : FCoords) : FCoords :=
  let ell := live j d
  let a := maximumMoment j zp
  let avg := liveAverage v
  { Z := Fix.add (Fix.mul (Fix.mul (Fix.sub (Fix.ofInt 1) ell) (zp j)) v.Z)
      (Fix.mul (Fix.mul ell (zp j)) v.D)
    D := Fix.add
      (Fix.mul (Fix.add (Fix.divInt ((boundsD j d zp).foldr iminFix a) 2)
        (Fix.divInt (iminFix a ell) 1048574)) v.Z)
      (Fix.divInt (Fix.mul a (Fix.add v.D avg)) 2)
    S i := let m := moment (levels.getD i 0) j zp
      Fix.add
        (Fix.mul (Fix.add (Fix.divInt ((boundsS j d zp i).foldr iminFix m) 2)
          (Fix.divInt (iminFix m ell) 1048574)) v.Z)
        (Fix.divInt (Fix.mul m (Fix.add avg
          (if polyChoose 128 j=d.kernel then v.S i else v.D))) 2) }

/-- Round the complete binomial mass once; powers cannot underflow before multiplication by choose. -/
def probability (j : Nat) (qn qd : Int) : Fix :=
  Fix.ofFrac ((polyChoose 128 j:Int)*qn^j*(qd-qn)^(128-j)) (qd^128)

def column (qn qd : Int) (zp : Nat → Fix) (v : FCoords) : FCoords where
  Z := isum ((List.range 129).map fun j => Fix.mul (probability j qn qd)
    (fixedColumn j (Data.weights.getD j Data.weight0) zp v).Z)
  D := isum ((List.range 129).map fun j => Fix.mul (probability j qn qd)
    (fixedColumn j (Data.weights.getD j Data.weight0) zp v).D)
  S i := isum ((List.range 129).map fun j => Fix.mul (probability j qn qd)
    ((fixedColumn j (Data.weights.getD j Data.weight0) zp v).S i))

end Spin.Structured.DenseOccupationFixed
