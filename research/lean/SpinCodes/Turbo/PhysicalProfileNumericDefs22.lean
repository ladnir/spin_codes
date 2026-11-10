import SpinCodes.Numeric.FixedDefs

/-!
Kernel-reducible integer enclosures for the small packet-profile coefficient
calculation. Keep this file free of Mathlib imports: historical certificate
data checks depend on that separation. Every linear factor is processed in
one pass over the stored coefficient list, so no coefficient is recursively
recomputed through the expansion tree.
-/
namespace Spin.Turbo
open Spin.Numeric

def profileFixLinearAux22 (a b prev : Fix) : List Fix → List Fix
  | [] => [Fix.mul b prev]
  | c :: cs => Fix.add (Fix.mul a c) (Fix.mul b prev) :: profileFixLinearAux22 a b c cs

def profileFixLinear22 (a b : Fix) (cs : List Fix) : List Fix :=
  profileFixLinearAux22 a b Fix.zero cs

def profileFixCoefficients22 : List (Fix × Fix) → List Fix
  | [] => [Fix.ofInt 1]
  | (a,b) :: factors => profileFixLinear22 a b (profileFixCoefficients22 factors)

def profileFixFactors22 (profile : List Nat) (bank : List (Fix × Fix)) : List (Fix × Fix) :=
  ((List.range 5).map (fun r => List.replicate (profile.getD r 0)
    (bank.getD r (Fix.zero,Fix.zero)))).flatten

def profileFixMomentBank22 (z : Fix) : List (Fix × Fix) :=
  let one := Fix.ofInt 1
  let z2 := Fix.mul z z
  let z3 := Fix.mul z2 z
  let z4 := Fix.mul z2 z2
  let total := Fix.pow (Fix.add one z) 4
  [(one,Fix.sub total one), (z,Fix.sub total z), (z2,Fix.sub total z2),
    (z3,Fix.sub total z3), (z4,Fix.sub total z4)]

def profileFixFourierBank22 (z : Fix) : List (Fix × Fix) :=
  let one := Fix.ofInt 1
  let plus := Fix.add one z
  let minus := Fix.sub one z
  (List.range 5).map (fun r => (one, Fix.sub
    (Fix.mul (Fix.pow plus (4-r)) (Fix.pow minus r)) one))

def profileChoose16Values22 : List Nat :=
  [1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1]

def profileLayerDenominator22 (j : Nat) : Int :=
  (profileChoose16Values22.getD j 0 : Int) * (15 : Int)^j

def profileFixMoment22 (profile : List Nat) (z : Fix) (j : Nat) : Fix :=
  Fix.divInt ((profileFixCoefficients22
    (profileFixFactors22 profile (profileFixMomentBank22 z))).getD j Fix.zero)
    (profileLayerDenominator22 j)

def profileFixFourier22 (profile : List Nat) (z : Fix) (j : Nat) : Fix :=
  Fix.divInt ((profileFixCoefficients22
    (profileFixFactors22 profile (profileFixFourierBank22 z))).getD j Fix.zero)
    (profileLayerDenominator22 j)

def profileFixSecondMoment22 (profile : List Nat) (z : Fix) (j : Nat) : Fix :=
  profileFixMoment22 profile (Fix.mul z z) j

/-- Exact cross multiplication, with positivity checked by the soundness
theorem. No floating point comparison and no reduction of rational fields. -/
def profileFixUpperCheck22 (x : Fix) (num den : Int) : Bool :=
  decide (x.hi * den ≤ num * Spin.Numeric.scale)

def profileFixLowerCheck22 (x : Fix) (num den : Int) : Bool :=
  decide (num * Spin.Numeric.scale ≤ x.lo * den)

end Spin.Turbo
