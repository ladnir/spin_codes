import SpinCodes.Structured.DenseScalarExactSound
import SpinCodes.Structured.DenseOccupationFixedRate

noncomputable section
namespace Spin.Structured.DenseScalarExact
open ConcreteRoute ConcreteRoutedEncoder ConcreteScalar DenseOccupationFixed

/-- The scalar alternative retains the exact route and likelihood costs. -/
theorem routed_probability_radius {L b R : ℕ} (hL : 0 < L) (hb : 0 < b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0 < totalWeight rows) (hw1 : totalWeight rows < L*b)
    {p z radius : ℝ} (hp0 : 0 < p) (hp1 : p < 1) (hz : 0 < z) (hz1 : z ≤ 1)
    (hr : scalarBound p z ≤ radius) (d : ℕ) :
    (experimentLaw L b R).prob (fun ω => weight e rows ω ≤ d) ≤
      ((((b : ℝ)+1)^activeRows rows * ((L : ℝ)+1)^b) *
        Real.exp (((L : ℝ)*b)*binKL (density rows) p)) * radius^R / z^d := by
  apply (routed_scalar_probability hL hb e rows hw0 hw1 hp0 hp1 hz hz1 d).trans
  apply div_le_div_of_nonneg_right _ (pow_nonneg hz.le _)
  apply mul_le_mul_of_nonneg_left _ (by positivity)
  exact pow_le_pow_left₀ (scalarBound_nonneg hp0.le hp1.le hz.le) hr R

/-- A certified scalar box gives the counted finite rate with prefactor one. -/
theorem routed_counted_rate {L b R d : ℕ} (hL : 0<L) (hb : 0<b)
    (e : (Fin b × Fin L) ≃ (Fin R × Fin 128)) (rows : Fin L → Finset (Fin b))
    (hw0 : 0<totalWeight rows) (hw1 : totalWeight rows<L*b)
    {α x p y radius z m c η C : ℝ}
    (hα0 : 0<α) (hα1 : α≤1) (hx0 : 0<x) (hx1 : x≤1)
    (hden : density rows=α*x) (hp0 : 0<p) (hp1 : p<1) (hy0 : 0<y) (hy1 : y<1)
    (hr : 0<radius) (hz0 : 0<z) (hz1 : z≤1)
    (hscalar : scalarBound (p*y) z ≤ radius)
    (hround : 128*R=L*b) (hd : (d:ℝ)≤(11/100)*((L:ℝ)*b))
    (hC0 : 0≤C) (hC : C≤Real.exp (((L:ℝ)*b)*α*(m*x+c)))
    (hbox : boxExponent m c p y radius z α x≤-η) :
    C*(experimentLaw L b R).prob (fun ω => weight e rows ω≤d) ≤
      (((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*Real.exp (-η*((L:ℝ)*b)) := by
  have hpY0 : 0<p*y := mul_pos hp0 hy0
  have hpY1 : p*y<1 := by nlinarith
  have hprob := routed_probability_radius hL hb e rows hw0 hw1 hpY0 hpY1 hz0 hz1 hscalar d
  have hs := finite_exponent_le (N := (L:ℝ)*b) (m := m) (c := c) (R := R) (d := d)
    (by positivity) hα0 hα1 hx0 hx1 hp0 hp1 hy0 hy1 hr hz0 hz1
    (by exact_mod_cast hround) hd
  have hs' : Real.exp (((L:ℝ)*b)*α*(m*x+c))*
      Real.exp (((L:ℝ)*b)*binKL (α*x) (p*y))*radius^R/z^d ≤
      Real.exp (-η*((L:ℝ)*b)) := hs.trans (Real.exp_le_exp.mpr (by
        have := mul_le_mul_of_nonneg_left hbox (show 0≤(L:ℝ)*b by positivity)
        nlinarith))
  calc
    _ ≤ C*(((((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*
      Real.exp (((L:ℝ)*b)*binKL (density rows) (p*y)))*radius^R/z^d) :=
        mul_le_mul_of_nonneg_left hprob hC0
    _ ≤ Real.exp (((L:ℝ)*b)*α*(m*x+c))*
      (((((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*
      Real.exp (((L:ℝ)*b)*binKL (density rows) (p*y)))*radius^R/z^d) := by
        apply mul_le_mul_of_nonneg_right hC
        positivity
    _ = (((b:ℝ)+1)^activeRows rows*((L:ℝ)+1)^b)*
      (Real.exp (((L:ℝ)*b)*α*(m*x+c))*Real.exp (((L:ℝ)*b)*binKL (α*x) (p*y))*radius^R/z^d) := by
        rw [hden]
        ring
    _ ≤ _ := mul_le_mul_of_nonneg_left hs' (by positivity)

#print axioms routed_probability_radius
#print axioms routed_counted_rate

end Spin.Structured.DenseScalarExact
