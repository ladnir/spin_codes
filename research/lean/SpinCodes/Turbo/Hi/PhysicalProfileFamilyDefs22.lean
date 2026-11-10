import SpinCodes.Turbo.Hi.PhysicalProfileNumericCacheDefs22

namespace Spin.Turbo.Hi
open Spin.NumericHi

def profileCachedCapsCheck22 (entries : List (Nat × List Nat)) (cache : ProfileFixCache22)
    (bank : List (Fix × Fix)) (j : Nat) (caps : List Int) (den : Int) : Bool :=
  decide (caps.length = 11) && caps.all (fun n => decide (0 ≤ n)) &&
    (List.range 11).all (fun i =>
      if i = 2 then profileFixUpperCheck22
        (profileFixClassMean22 entries i (fun p => profileFixCachedMoment22 cache bank p j))
        (caps.getD i 0) den
      else profileFixClassUpperCheck22 entries i
        (fun p => profileFixCachedMoment22 cache bank p j) (caps.getD i 0) den)

def profileReturnCapsCheck22 (first second returns : List Int) (active : Int) : Bool :=
  decide (returns.length = 11) && returns.all (fun n => decide (0 ≤ n)) &&
    (List.range 10).all (fun k =>
      let i := k+1
      let r := returns.getD i 0 * 4194303
      decide (first.getD i 0 ≤ r ∨ second.getD i 0 * active ≤ r*r))

def profileFamilyRowCheck22 (entries : List (Nat × List Nat))
    (firstCache secondCache zeroCache feedbackCache : ProfileFixCache22)
    (firstBank secondBank zeroBank feedbackBank : List (Fix × Fix))
    (j : Nat) (first second returns : List Int) (active zeroHold den : Int) : Bool :=
  decide (0 < den ∧ 0 ≤ active ∧ active ≤ den ∧ 0 ≤ zeroHold) &&
    profileCachedCapsCheck22 entries firstCache firstBank j first den &&
    profileCachedCapsCheck22 entries secondCache secondBank j second den &&
    profileFixLowerCheck22 (profileFixCachedAllFourier22 entries feedbackCache feedbackBank j)
      (den-active) den &&
    profileFixUpperCheck22 (profileFixCachedAllFourier22 entries zeroCache zeroBank j) zeroHold den &&
    profileReturnCapsCheck22 first second returns active

def profileFamilyCheck22 (entries : List (Nat × List Nat)) (box : Fix)
    (first second returns : List (List Int)) (active zeroHold : List Int) (den : Int) : Bool :=
  let firstBank := profileFixMomentBank22 box
  let secondBank := profileFixMomentBank22 (Fix.mul box box)
  let zeroBank := profileFixFourierBank22 box
  let feedbackBank := profileFixFourierBank22 (Fix.ofInt 1)
  let firstCache := profileFixCache22 entries firstBank
  let secondCache := profileFixCache22 entries secondBank
  let zeroCache := profileFixCache22 entries zeroBank
  let feedbackCache := profileFixCache22 entries feedbackBank
  decide (first.length = 17 ∧ second.length = 17 ∧ returns.length = 17 ∧
    active.length = 17 ∧ zeroHold.length = 17) &&
    (List.range 17).all (fun j => profileFamilyRowCheck22 entries
      firstCache secondCache zeroCache feedbackCache firstBank secondBank zeroBank feedbackBank
      j (first.getD j []) (second.getD j []) (returns.getD j [])
      (active.getD j 0) (zeroHold.getD j 0) den)

end Spin.Turbo.Hi
