import SpinCodes.Structured.PolyIdentityDefs
import SpinCodes.Structured.SparsePowers

/-! Rational polynomial program for the sparse verifier. All branch choices
are explicit finite data. Their domination checks and the semantic link to
the real occupation matrix are separate from this executable program. -/

namespace Spin.Structured.SparsePolynomial

open RatPoly

structure Choices where
  momentMax : Nat
  lowMax : Nat
  cancelD : Nat
  freshD : Nat
  cancelS : List Nat
  freshS : List Nat

structure WeightData where
  kernel : Nat
  cap : Nat
  lowPatterns : Option (List (List (Nat × Nat)))
  lowShells : List (List (Nat × Nat))
  choices : Choices

def levels : List Nat := [48, 56, 64, 72, 80]
def counts : List Nat := [5166, 110288, 293455, 110128, 5250]
def zero : RatPoly := constant 0
def one : RatPoly := constant 1
def z : RatPoly := ⟨[6250, -1], 6250⟩
def beta : RatPoly := ⟨[0, 1], 12500⟩
def complement : RatPoly := ⟨[12500, -1], 12500⟩

def witness : Nat → RatPoly
  | 0 => one
  | 1 => ⟨[100, 9], 102400⟩
  | 2 => ⟨[22570555350000, 116184452864], 23112248678400000⟩
  | 3 => ⟨[36139102910000, 92653723520], 37006441379840000⟩
  | 4 => ⟨[7692732079250000, -26853143552], 7877357649152000000⟩
  | 5 => ⟨[180433371050000, -461216096896], 184763771955200000⟩
  | _ => ⟨[68812668750000, -351546952448], 70464172800000000⟩

def scale (a : Int) (d : Nat) (p : RatPoly) : RatPoly := mul (constant a d) p

def moment (w j : Nat) : RatPoly :=
  scale 1 (polyChoose 128 j) (sum ((List.range 129).map fun h =>
    let c := polyChoose w h * polyChoose (128 - w) (j - h)
    if h ≤ j ∧ c ≠ 0 then scale c 1 (SparsePowers.zPow (w + j - 2 * h)) else zero))

def pattern (j : Nat) (p : List (Nat × Nat)) : RatPoly :=
  scale 1 (polyChoose 128 j) (sum (p.map fun (out, c) => scale c 1 (SparsePowers.zPow out)))

def live (j : Nat) (d : WeightData) : RatPoly :=
  constant (polyChoose 128 j - d.kernel) (polyChoose 128 j)

def moments (j : Nat) : List RatPoly := levels.map (fun w => moment w j)

def arbitrary (j : Nat) (d : WeightData) : RatPoly :=
  (moments j).getD d.choices.momentMax zero

def distance (a b : Nat) : Nat := (a - b) + (b - a)

def minDistance (j : Nat) : Nat :=
  (levels.map (fun w => distance w j)).foldr min 128

def boundsD (j : Nat) (d : WeightData) : List RatPoly :=
  [arbitrary j d, live j d, scale d.cap (polyChoose 128 j) (SparsePowers.zPow (minDistance j))] ++
  match d.lowPatterns with
  | none => []
  | some ps => [(ps.map (pattern j)).getD d.choices.lowMax zero]

def boundsS (j : Nat) (d : WeightData) (i : Nat) : List RatPoly :=
  let w := levels.getD i 0
  let c := counts.getD i 1
  [moment w j,
   scale (min (polyChoose 128 j - d.kernel) (c * d.cap)) (c * polyChoose 128 j)
     (SparsePowers.zPow (distance w j))] ++
  match d.lowPatterns with
  | none => []
  | some _ => [scale 1 c (pattern j (d.lowShells.getD i []))]

def action (j : Nat) (d : WeightData) : Nat → RatPoly
  | 0 => add (scale d.kernel (polyChoose 128 j) (SparsePowers.zPow j))
      (mul (mul (live j d) (SparsePowers.zPow j)) (witness 1))
  | 1 =>
      let a := arbitrary j d
      let c := (boundsD j d).getD d.choices.cancelD zero
      let f := [a, live j d].getD d.choices.freshD zero
      add (add (scale 1 2 c) (scale 1 1048574 f))
        (scale 1 2 (mul a (add (witness 1) (constant 1 1024))))
  | i + 2 =>
      let m := moment (levels.getD i 0) j
      let c := (boundsS j d i).getD (d.choices.cancelS.getD i 0) zero
      let f := [m, live j d].getD (d.choices.freshS.getD i 0) zero
      let lazy := if polyChoose 128 j = d.kernel then witness (i + 2) else witness 1
      add (add (scale 1 2 c) (scale 1 1048574 f))
        (scale 1 2 (mul m (add (constant 1 1024) lazy)))

def probability (j : Nat) : RatPoly :=
  scale (polyChoose 128 j) 1 (mul (SparsePowers.betaPow j) (SparsePowers.complementPow (128 - j)))

def contribution (j : Nat) (d : WeightData) (i : Nat) : RatPoly :=
  mul (probability j) (action j d i)

end Spin.Structured.SparsePolynomial
