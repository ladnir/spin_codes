import SpinCodes.Turbo.Hi.PhysicalProfileNumericDefs22

namespace Spin.Turbo.Hi
open Spin.NumericHi

def profileListWeight22 (profile : List Nat) : Nat :=
  (List.range 5).foldl (fun w r => w + r * profile.getD r 0) 0

def profileListClassAllows22 (i : Nat) (profile : List Nat) : Bool :=
  let w := profileListWeight22 profile
  if i = 0 then decide (w = 0) else
  if i = 1 ∨ i = 2 then decide (w ≠ 0) else
    decide (w = ([0,0,0,16,24,28,32,36,40,48,64] : List Nat).getD i 0)

def profileFixClassSum22 (entries : List (Nat × List Nat)) (i : Nat)
    (value : List Nat → Fix) : Fix :=
  entries.foldr (fun e acc => if profileListClassAllows22 i e.2 then
    Fix.add (Fix.mul (Fix.ofInt (e.1 : Int)) (value e.2)) acc else acc) Fix.zero

def profileListClassCount22 (entries : List (Nat × List Nat)) (i : Nat) : Nat :=
  entries.foldr (fun e acc => if profileListClassAllows22 i e.2 then e.1+acc else acc) 0

def profileFixClassMean22 (entries : List (Nat × List Nat)) (i : Nat)
    (value : List Nat → Fix) : Fix :=
  Fix.divInt (profileFixClassSum22 entries i value) (profileListClassCount22 entries i : Int)

def profileFixAllSum22 (entries : List (Nat × List Nat)) (value : List Nat → Fix) : Fix :=
  entries.foldr (fun e acc => Fix.add (Fix.mul (Fix.ofInt (e.1 : Int)) (value e.2)) acc) Fix.zero

def profileFixAllMean22 (entries : List (Nat × List Nat)) (value : List Nat → Fix) : Fix :=
  Fix.divInt (profileFixAllSum22 entries value) 4194304

def profileFixClassUpperCheck22 (entries : List (Nat × List Nat)) (i : Nat)
    (value : List Nat → Fix) (num den : Int) : Bool :=
  entries.all (fun e => if profileListClassAllows22 i e.2 then profileFixUpperCheck22 (value e.2) num den else true)

def profileFixRawMoment22 (profile : List Nat) (box : Fix) (j : Nat) : Fix :=
  (profileFixCoefficients22 (profileFixFactors22 profile (profileFixMomentBank22 box))).getD j Fix.zero

def profileFixRawFourier22 (profile : List Nat) (box : Fix) (j : Nat) : Fix :=
  (profileFixCoefficients22 (profileFixFactors22 profile (profileFixFourierBank22 box))).getD j Fix.zero

/-- Common denominators are applied after summing, avoiding interval widening
from separately normalized signed terms. -/
def profileFixAllFourier22 (entries : List (Nat × List Nat)) (box : Fix) (j : Nat) : Fix :=
  Fix.divInt (profileFixAllSum22 entries (fun p => profileFixRawFourier22 p box j))
    (4194304 * profileLayerDenominator22 j)

def profileFixClassMomentMean22 (entries : List (Nat × List Nat)) (i : Nat) (box : Fix) (j : Nat) : Fix :=
  Fix.divInt (profileFixClassSum22 entries i (fun p => profileFixRawMoment22 p box j))
    ((profileListClassCount22 entries i : Int) * profileLayerDenominator22 j)

/-- Compute a whole row once; callers checking all seventeen occupations
should retain this list instead of rebuilding the recurrence for each j. -/
def profileFixRowsFromBank22 (profile : List Nat) (bank : List (Fix × Fix)) : List Fix :=
  (profileFixCoefficients22 (profileFixFactors22 profile bank)).mapIdx
    (fun j c => Fix.divInt c (profileLayerDenominator22 j))

end Spin.Turbo.Hi
