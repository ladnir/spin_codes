import SpinCodes.Structured.PolyIdentityDefs

/-! Executable integer formulas for the concrete Fourier certificates. -/

namespace Spin.Structured.FiberNumerics

def kraw (j w : Nat) : Int :=
  ((List.range (j + 1)).map fun h =>
    (-1 : Int) ^ h * (polyChoose w h : Int) * (polyChoose (128 - w) (j - h) : Int)).sum

def transform (spectrum : List Nat) (values : List Int) (f : Int → Int) : Int :=
  ((List.range 129).map fun w => (spectrum.getD w 0 : Int) * f (values.getD w 0)).sum

def signedSum (spectrum : List Nat) (values : List Int) : Int :=
  transform spectrum values id

def squareSum (spectrum : List Nat) (values : List Int) : Int :=
  transform spectrum values (fun k => k * k)

def absSum (spectrum : List Nat) (values : List Int) : Int :=
  transform spectrum values (fun k => (k.natAbs : Int))

end Spin.Structured.FiberNumerics
