import SpinCodes.Turbo.Hi.PhysicalProfileNumericAggregateDefs22

namespace Spin.Turbo.Hi
open Spin.NumericHi

abbrev ProfileFixCache22 := List (List Nat × List Fix)

def profileFixCache22 (entries : List (Nat × List Nat)) (bank : List (Fix × Fix)) : ProfileFixCache22 :=
  entries.map (fun e => (e.2,profileFixCoefficients22 (profileFixFactors22 e.2 bank)))

def profileFixCachedRaw22 (cache : ProfileFixCache22) (bank : List (Fix × Fix))
    (profile : List Nat) (j : Nat) : Fix :=
  match cache.lookup profile with
  | some cs => cs.getD j Fix.zero
  | none => (profileFixCoefficients22 (profileFixFactors22 profile bank)).getD j Fix.zero

def profileFixCachedMoment22 (cache : ProfileFixCache22) (bank : List (Fix × Fix))
    (profile : List Nat) (j : Nat) : Fix :=
  Fix.divInt (profileFixCachedRaw22 cache bank profile j) (profileLayerDenominator22 j)

def profileFixCachedAllFourier22 (entries : List (Nat × List Nat))
    (cache : ProfileFixCache22) (bank : List (Fix × Fix)) (j : Nat) : Fix :=
  Fix.divInt (profileFixAllSum22 entries (fun p => profileFixCachedRaw22 cache bank p j))
    (4194304 * profileLayerDenominator22 j)

end Spin.Turbo.Hi
