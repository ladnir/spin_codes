import SpinCodes.Structured.DenseScalarExactRate
import SpinCodes.Structured.DenseOccupationFamilyRateDefs

noncomputable section
namespace Spin.Structured.DenseScalarExact
open DenseOccupationFixed

/-- Scalar transfer and a local exponent bound supply the common actual-experiment rate. -/
theorem pointRate {α x p y radius z m c η : ℝ}
    (hα0 : 0<α) (hα1 : α≤1) (hx0 : 0<x) (hx1 : x≤1)
    (hp0 : 0<p) (hp1 : p<1) (hy0 : 0<y) (hy1 : y<1)
    (hr : 0<radius) (hz0 : 0<z) (hz1 : z≤1)
    (hscalar : ConcreteScalar.scalarBound (p*y) z ≤ radius)
    (hbox : boxExponent m c p y radius z α x≤-η) :
    PointRate α x m c η 1 := by
  intro L b R d hL hb e rows hw0 hw1 hden hround hd K hK0 hK
  simpa only [mul_one] using routed_counted_rate hL hb e rows hw0 hw1
    hα0 hα1 hx0 hx1 hden hp0 hp1 hy0 hy1 hr hz0 hz1 hscalar hround hd hK0 hK hbox

#print axioms pointRate
end Spin.Structured.DenseScalarExact
