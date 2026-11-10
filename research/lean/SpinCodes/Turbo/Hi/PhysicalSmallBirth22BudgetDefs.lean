import SpinCodes.Turbo.PhysicalSmallBirth22NumericData
import SpinCodes.Turbo.Hi.PhysicalProfileNumericDefs22

/-!
Integer-only birth-budget evaluation. The historical FixedDefs Horner
pattern avoids repeated powers. Each nine-coefficient class row is read
once after dropping its flat-table prefix; no state-space enumeration occurs.
-/
namespace Spin.Turbo.Hi
open Spin.NumericHi

def smallBirthFixClassCoefficients22 (hist : List Nat) (i : Nat) : List Nat :=
  (hist.drop (9*i)).take 9

def smallBirthFixHorner22 : List Nat → Fix → Fix
  | [], _ => Fix.zero
  | c :: cs, box => Fix.add (Fix.ofInt (c : Int)) (Fix.mul box (smallBirthFixHorner22 cs box))

def smallBirthFixBudget22 (hist : List Nat) (norm : Int) (box : Fix) (i : Nat) : Fix :=
  Fix.divInt (smallBirthFixHorner22 (smallBirthFixClassCoefficients22 hist i) box) norm

def smallBirthFixOneBudget22 (box : Fix) (i : Nat) : Fix :=
  smallBirthFixBudget22 smallBirthOneHistogram22 240 box i

def smallBirthFixTwoBudget22 (box : Fix) (i : Nat) : Fix :=
  smallBirthFixBudget22 smallBirthTwoHistogram22 27000 box i

def smallBirthFixZeroBudget22 (i : Nat) : Fix :=
  if i = 0 then Fix.ofInt 1 else Fix.zero

def smallBirthFixVectorUpperCheck22 (values : Nat → Fix) (nums : List Int) (den : Int) : Bool :=
  decide (0 < den) && (List.range 11).all (fun i =>
    profileFixUpperCheck22 (values i) (nums.getD i 0) den)

def smallBirthFixOneBudgetCheck22 (box : Fix) (nums : List Int) (den : Int) : Bool :=
  decide (0 ≤ box.lo) && smallBirthFixVectorUpperCheck22 (smallBirthFixOneBudget22 box) nums den

def smallBirthFixTwoBudgetCheck22 (box : Fix) (nums : List Int) (den : Int) : Bool :=
  decide (0 ≤ box.lo) && smallBirthFixVectorUpperCheck22 (smallBirthFixTwoBudget22 box) nums den

def smallBirthFixZeroBudgetCheck22 (nums : List Int) (den : Int) : Bool :=
  smallBirthFixVectorUpperCheck22 smallBirthFixZeroBudget22 nums den

end Spin.Turbo.Hi
