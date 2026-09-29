import SpinCodes.Structured.ConcreteFourierRouted
import SpinCodes.Structured.DenseOccupationFixedRate
import SpinCodes.Structured.DenseOccupationFamilyRateDefs

noncomputable section
namespace Spin.Structured.DenseFourierExact
open Spin.Imt ConcreteRoute ConcreteRoutedEncoder DenseOccupationFixed

/-- Finite positive-occupation rate; the outer counting cost remains an explicit input. -/
theorem routed_counted_rate {L b R d : Nat} (hL : 0<L) (hb : 0<b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0<totalWeight rows) (hw1 : totalWeight rows<L*b)
    {α x p y radius z m c η C wmin : ℝ} {w : Coords 5}
    (hα0 : 0<α) (hα1 : α≤1) (hx0 : 0<x) (hx1 : x≤1)
    (hden : density rows=α*x) (hp0 : 0<p) (hp1 : p<1) (hy0 : 0<y) (hy1 : y<1)
    (hr : 0<radius) (hz0 : 0<z) (hz1 : z≤1)
    (hmin : 0<wmin) (hZ : wmin≤w.Z) (hD : wmin≤w.D) (hS : ∀ i,wmin≤w.S i)
    (hcol : ((ConcreteFourier.matrix (p*y) z).applyCol w).le (Coords.smul radius w))
    (hround : 128*R=L*b) (hd : (d:ℝ)≤(11/100)*((L:ℝ)*b))
    (hC0 : 0≤C) (hC : C≤Real.exp (((L:ℝ)*b)*α*(m*x+c)))
    (hbox : boxExponent m c p y radius z α x≤-η) :
    C*(experimentLaw L b R).prob (fun ω => weight e rows ω≤d) ≤
      (((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*(w.Z/wmin)*Real.exp (-η*((L:ℝ)*b)) := by
  have hpY0 : 0<p*y := mul_pos hp0 hy0
  have hpY1 : p*y<1 := by nlinarith
  have hprob := ConcreteFourier.routed_probability hL hb e rows hw0 hw1 hpY0 hpY1
    hz0 hz1 hmin hZ hD hS hcol d
  have hpre0 : 0≤(((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*(w.Z/wmin) := by
    have : 0≤w.Z := hmin.le.trans hZ
    positivity
  have hs := finite_exponent_le (N := (L:ℝ)*b) (m := m) (c := c) (R := R) (d := d)
    (by positivity) hα0 hα1 hx0 hx1 hp0 hp1 hy0 hy1 hr hz0 hz1
    (by exact_mod_cast hround) hd
  have hs' : Real.exp (((L:ℝ)*b)*α*(m*x+c))*
      Real.exp (((L:ℝ)*b)*binKL (α*x) (p*y))*radius^R/z^d ≤
      Real.exp (-η*((L:ℝ)*b)) := hs.trans (Real.exp_le_exp.mpr (by
        have := mul_le_mul_of_nonneg_left hbox (show 0≤(L:ℝ)*b by positivity)
        nlinarith))
  calc
    C*(experimentLaw L b R).prob (fun ω => weight e rows ω≤d)
      ≤ C*(((((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*
        Real.exp (((L:ℝ)*b)*binKL (density rows) (p*y)))*(radius^R*w.Z/wmin)/z^d) :=
          mul_le_mul_of_nonneg_left hprob hC0
    _ ≤ Real.exp (((L:ℝ)*b)*α*(m*x+c))*
        (((((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*
        Real.exp (((L:ℝ)*b)*binKL (density rows) (p*y)))*(radius^R*w.Z/wmin)/z^d) := by
          apply mul_le_mul_of_nonneg_right hC
          have : 0≤w.Z := hmin.le.trans hZ
          positivity
    _ = ((((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*(w.Z/wmin))*
        (Real.exp (((L:ℝ)*b)*α*(m*x+c))*Real.exp (((L:ℝ)*b)*binKL (α*x) (p*y))*radius^R/z^d) := by
          rw [hden]
          ring
    _ ≤ _ := mul_le_mul_of_nonneg_left hs' hpre0

/-- A concrete Fourier matrix and a checked box exponent give the common family interface. -/
theorem fourier_pointRate {α x p y radius z m c η wmin : ℝ} {w : Coords 5}
    (hα0 : 0 < α) (hα1 : α ≤ 1) (hx0 : 0 < x) (hx1 : x ≤ 1)
    (hp0 : 0 < p) (hp1 : p < 1) (hy0 : 0 < y) (hy1 : y < 1)
    (hr : 0 < radius) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (hmin : 0 < wmin) (hZ : wmin ≤ w.Z) (hD : wmin ≤ w.D) (hS : ∀ i, wmin ≤ w.S i)
    (hcol : ((ConcreteFourier.matrix (p*y) z).applyCol w).le (Coords.smul radius w))
    (hbox : boxExponent m c p y radius z α x ≤ -η) :
    PointRate α x m c η (w.Z/wmin) := by
  intro L b R d hL hb e rows hw0 hw1 hden hround hd K hK0 hK
  exact routed_counted_rate hL hb e rows hw0 hw1 hα0 hα1 hx0 hx1 hden hp0 hp1 hy0 hy1
    hr hz0 hz1 hmin hZ hD hS hcol hround hd hK0 hK hbox


#print axioms routed_counted_rate
#print axioms fourier_pointRate
end Spin.Structured.DenseFourierExact
