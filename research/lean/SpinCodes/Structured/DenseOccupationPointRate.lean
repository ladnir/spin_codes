import SpinCodes.Structured.DenseOccupationFamilyRateDefs
import SpinCodes.Structured.DenseOccupationFixedRate

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Spin.Imt

/-- A concrete occupation matrix and a checked box exponent give the common family interface. -/
theorem occupation_pointRate {α x p y radius z m c η wmin : ℝ} {w : Coords 5}
    (hα0 : 0 < α) (hα1 : α ≤ 1) (hx0 : 0 < x) (hx1 : x ≤ 1)
    (hp0 : 0 < p) (hp1 : p < 1) (hy0 : 0 < y) (hy1 : y < 1)
    (hr : 0 < radius) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (hmin : 0 < wmin) (hZ : wmin ≤ w.Z) (hD : wmin ≤ w.D) (hS : ∀ i, wmin ≤ w.S i)
    (hcol : ((Occupation.Sparse.numericalMatrix (p*y) z).applyCol w).le (Coords.smul radius w))
    (hbox : boxExponent m c p y radius z α x ≤ -η) :
    PointRate α x m c η (w.Z/wmin) := by
  intro L b R d hL hb e rows hw0 hw1 hden hround hd K hK0 hK
  exact routed_counted_rate hL hb e rows hw0 hw1 hα0 hα1 hx0 hx1 hden hp0 hp1 hy0 hy1
    hr hz0 hz1 hmin hZ hD hS hcol hround hd hK0 hK hbox

#print axioms occupation_pointRate
end Spin.Structured.DenseOccupationFixed
