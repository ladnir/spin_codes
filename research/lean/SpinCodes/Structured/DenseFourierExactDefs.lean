import SpinCodes.Structured.DenseScalarExactDefs
import SpinCodes.Structured.FiberNumericsData.Tables
import SpinCodes.Structured.SparseWeights

namespace Spin.Structured.DenseFourierExact
open DenseScalarExact

def u0 (qn qd zn zd : Int) : Int := Int.ofNat (((qd-qn)*zd-qn*zn).natAbs)
def u1 (qn qd zn zd : Int) : Int := Int.ofNat ((qn*zd-(qd-qn)*zn).natAbs)

/-- The emission factor clears both Fourier ratio denominators. -/
def monomial (qn qd zn zd : Int) (d w h : Nat) : Int :=
  g0 qn qd zn zd^(128-d-(w-h))*g1 qn qd zn zd^(d-h)*
    u0 qn qd zn zd^(w-h)*u1 qn qd zn zd^h

def overlap (qn qd zn zd : Int) (d w : Nat) : Int :=
  max (monomial qn qd zn zd d w (d+w-128)) (monomial qn qd zn zd d w (min d w))

def spectrum (qn qd zn zd : Int) (d : Nat) : Int :=
  ((List.range 129).map fun w => (FiberNumerics.Data.spectrum.getD w 0:Int)*overlap qn qd zn zd d w).sum

def target (qn qd zn zd : Int) (d : Nat) : Int :=
  524287*spectrum qn qd zn zd d+524288*(g0 qn qd zn zd^(128-d)*g1 qn qd zn zd^d)

def denominator (qd zd : Int) : Int := 2*524287*524288*(qd*zd)^128

end Spin.Structured.DenseFourierExact
