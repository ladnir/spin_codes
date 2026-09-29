import SpinCodes.Structured.DenseFourierExactCheck
import SpinCodes.Structured.DenseFourierExactRate

noncomputable section
namespace Spin.Structured.DenseFourierExact
open Spin.Imt

example {qn qd zn zd rn rd : Int} (hqd : 0 < qd) (hzd : 0 < zd) (hrd : 0 < rd)
    (hA : 0 < DenseScalarExact.g0 qn qd zn zd) (hB : 0 < DenseScalarExact.g1 qn qd zn zd)
    (v : ICoords) (hv : 0 ≤ weightedTotal v)
    (hZ : checkZ qn qd zn zd rn rd v = true)
    (hD : checkD qn qd zn zd rn rd v = true)
    (hS : ∀ i, checkS qn qd zn zd rn rd v i = true) :
    ((ConcreteFourier.matrix ((qn:ℝ)/qd) ((zn:ℝ)/zd)).applyCol v.real).le
      (Coords.smul ((rn:ℝ)/rd) v.real) :=
  checked_collatz hqd hzd hrd hA hB v hv hZ hD hS

#print axioms targetCap_rational
#print axioms zeroColumn_rational
#print axioms checked_collatz
#print axioms routed_counted_rate
#print axioms fourier_pointRate
end Spin.Structured.DenseFourierExact
