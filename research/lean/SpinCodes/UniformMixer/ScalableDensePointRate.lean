import SpinCodes.UniformMixer.ScalableDenseExponent

/-! One dense occupation point, assuming an explicit positive transfer witness. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteRoute

theorem scalable_dense_counted_rate (m d : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m)))
    (hw0 : 0 < totalWeight rows)
    (hw1 : totalWeight rows < Lsched m * bsched m)
    {α x p y ρ z slope intercept η C scale : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hx0 : 0 < x) (hx1 : x ≤ 1)
    (hp0 : 0 < p) (hp1 : p < 1)
    (hy0 : 0 < y) (hy1 : y < 1)
    (hρ : 0 < ρ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (v : Fin 3 → ℝ) (hscale : 0 ≤ scale)
    (hones : ∀ i, (1 : ℝ) ≤ scale * v i)
    (hv : ∀ i,
      (paperSpectrumTransfer (p * y) z).mulVec v i ≤ ρ * v i)
    (hden : density rows = α * x)
    (hd : (d : ℝ) ≤ (11 / 100) *
      ((Lsched m : ℝ) * bsched m))
    (hC0 : 0 ≤ C)
    (hC : C ≤ Real.exp (((Lsched m : ℝ) * bsched m) *
      α * (slope * x + intercept)))
    (hbox : scalableBoxExponent slope intercept p y ρ z α x ≤ -η) :
    C * ((law rows).prod (scalableStateMapLaw m)).prob
      (fun ω => sampledOutputWeight
        (scalableFlatInput m ω.1) (scalableFlatMaps m ω.2) 0 ≤ d) ≤
      (((bsched m : ℝ) + 1) ^ activeRows rows *
        ((Lsched m : ℝ) + 1) ^ bsched m) *
        (scale * v 0) *
        Real.exp (-η * ((Lsched m : ℝ) * bsched m)) := by
  have hpy0 : 0 < p * y := mul_pos hp0 hy0
  have hpy1 : p * y < 1 := by nlinarith
  have hprob := scalable_dense_tilted_probability m d rows hw0 hw1
    hpy0 hpy1 hz0 hz1 v hscale hρ.le hones hv
  have hN : (0 : ℝ) ≤ (Lsched m : ℝ) * bsched m := by positivity
  have hround : 64 * ((bsched m * scalableStepsPerRegion m : ℕ) : ℝ) =
      (Lsched m : ℝ) * bsched m := by
    exact_mod_cast (show 64 * (bsched m * scalableStepsPerRegion m) =
      Lsched m * bsched m by rw [scalable_region_length]; ring)
  have hs := scalable_dense_finite_exponent_le (m := slope) (c := intercept)
    hN hα0 hα1
    hx0 hx1 hp0 hp1 hy0 hy1 hρ hz0 hz1 hround hd
  have hs' :
      Real.exp (((Lsched m : ℝ) * bsched m) * α *
        (slope * x + intercept)) *
      Real.exp (((Lsched m : ℝ) * bsched m) *
        binKL (α * x) (p * y)) *
      ρ ^ (bsched m * scalableStepsPerRegion m) / z ^ d ≤
      Real.exp (-η * ((Lsched m : ℝ) * bsched m)) :=
    hs.trans (Real.exp_le_exp.mpr (by
      have := mul_le_mul_of_nonneg_left hbox hN
      nlinarith))
  have hfactor :
      0 ≤ (((bsched m : ℝ) + 1) ^ activeRows rows *
        ((Lsched m : ℝ) + 1) ^ bsched m) * (scale * v 0) := by
    have hv0 : 0 ≤ scale * v 0 := (by norm_num : (0 : ℝ) ≤ 1).trans (hones 0)
    positivity
  calc
    _ ≤ C * (((((bsched m : ℝ) + 1) ^ activeRows rows *
          ((Lsched m : ℝ) + 1) ^ bsched m) *
          Real.exp (((Lsched m : ℝ) * bsched m) *
            binKL (density rows) (p * y))) *
          (scale * ρ ^ (bsched m * scalableStepsPerRegion m) * v 0) /
          z ^ d) := mul_le_mul_of_nonneg_left hprob hC0
    _ ≤ Real.exp (((Lsched m : ℝ) * bsched m) * α *
          (slope * x + intercept)) *
          (((((bsched m : ℝ) + 1) ^ activeRows rows *
            ((Lsched m : ℝ) + 1) ^ bsched m) *
            Real.exp (((Lsched m : ℝ) * bsched m) *
              binKL (density rows) (p * y))) *
            (scale * ρ ^ (bsched m * scalableStepsPerRegion m) * v 0) /
            z ^ d) := by
      apply mul_le_mul_of_nonneg_right hC
      have hv0 : 0 ≤ scale * v 0 := (by norm_num : (0 : ℝ) ≤ 1).trans (hones 0)
      have hzpow : 0 < z ^ d := pow_pos hz0 _
      have hρpow : 0 ≤ ρ ^ (bsched m * scalableStepsPerRegion m) := pow_nonneg hρ.le _
      have hroute : 0 ≤ ((bsched m : ℝ) + 1) ^ activeRows rows *
        ((Lsched m : ℝ) + 1) ^ bsched m := by positivity
      have hinner : 0 ≤ scale *
          ρ ^ (bsched m * scalableStepsPerRegion m) * v 0 := by
        nlinarith [mul_nonneg hρpow hv0]
      positivity
    _ = ((((bsched m : ℝ) + 1) ^ activeRows rows *
          ((Lsched m : ℝ) + 1) ^ bsched m) * (scale * v 0)) *
          (Real.exp (((Lsched m : ℝ) * bsched m) * α *
            (slope * x + intercept)) *
            Real.exp (((Lsched m : ℝ) * bsched m) *
              binKL (α * x) (p * y)) *
            ρ ^ (bsched m * scalableStepsPerRegion m) / z ^ d) := by
      rw [hden]
      ring
    _ ≤ _ := by
      exact mul_le_mul_of_nonneg_left hs' hfactor

#print axioms scalable_dense_counted_rate

end Spin.UniformMixer
