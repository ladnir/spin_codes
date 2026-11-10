import SpinCodes.Turbo.PhysicalSmallBirth22NumericData

namespace Spin.Turbo
open Spin.Structured

def smallBirthOneCode22 (a : Fin 16 × Fin 15) : Nat :=
  9 * (smallBirthClassWeight22 (MapSpectrum.fastWeight 64
    (smallBirthSingleImages22.getD (15*a.1.val+a.2.val) 0))).val +
    MapSpectrum.fastWeight 4 (a.2.val+1)

def smallBirthTwoCode22 (a b : Fin 16) (u v : Fin 15) : Nat :=
  9 * (smallBirthClassWeight22 (MapSpectrum.fastWeight 64
    ((smallBirthSingleImages22.getD (15*a.val+u.val) 0) ^^^
      (smallBirthSingleImages22.getD (15*b.val+v.val) 0)))).val +
    (MapSpectrum.fastWeight 4 (u.val+1) + MapSpectrum.fastWeight 4 (v.val+1))

def smallBirthOneCodes22 : List Nat :=
  (List.ofFn (fun a : Fin 16 => List.ofFn (fun u : Fin 15 =>
    smallBirthOneCode22 (a,u)))).flatten

def smallBirthPairCodes22 (a b : Fin 16) : List Nat :=
  (List.ofFn (fun u : Fin 15 => List.ofFn (fun v : Fin 15 =>
    smallBirthTwoCode22 a b u v))).flatten
def smallBirthTwoShardCodes22 (a : Fin 16) : List Nat :=
  (List.ofFn (fun b : Fin 16 => if a < b then
    (List.ofFn (fun u : Fin 15 => List.ofFn (fun v : Fin 15 =>
      smallBirthTwoCode22 a b u v))).flatten else [])).flatten

def smallBirthCodeHistogram22 (codes : List Nat) : List Nat :=
  codes.foldl (fun h w => MapSpectrum.bump w h) (List.replicate 99 0)


end Spin.Turbo


