import SpinCodes.Structured.DenseOccupationFixedB093
import SpinCodes.Structured.DenseOccupationPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseOccupationFixed.B093
open Spin.Numeric DenseGeometry

def supportLine : ℚ × ℚ := ((m.num:ℚ)/m.den,(c.num:ℚ)/c.den)

lemma supportLine_mem : supportLine ∈ Spin.Majorant.refined.supports := by decide +kernel

lemma geometry_match : boxes.getD 235 zeroRect =
    ⟨(a0.num:ℚ)/a0.den,(a1.num:ℚ)/a1.den,(x0.num:ℚ)/x0.den,(x1.num:ℚ)/x1.den⟩ := by
  change (⟨241/10240,10031/320000,428613422794653/2000000000000000,499776708257287/2000000000000000⟩ : Rect) = _
  congr 1 <;> norm_num [a0,a1,x0,x1]

/-- The occupation witness certifies its original globally indexed rectangle. -/
theorem certified : CertifiedBox 235 (4/10000000) 847 := by
  intro α x hg
  rw [geometry_match] at hg
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [Rect.Contains,a0,a1,x0,x1,QInput.real] at hg ⊢
    exact ⟨hg.1,hg.2.1⟩
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [Rect.Contains,a0,a1,x0,x1,QInput.real] at hg ⊢
    exact ⟨hg.2.2.1,hg.2.2.2⟩
  have hf := W043.witness_floor
  have hrate : PointRate α x m.real c.real (4/10000000) (W043.w.Z/(1/847)) := occupation_pointRate
    (by norm_num [a0,QInput.real] at hα; linarith [hα.1] : 0 < α)
    (by norm_num [a1,QInput.real] at hα; linarith [hα.2] : α ≤ 1)
    (by norm_num [x0,QInput.real] at hx; linarith [hx.1] : 0 < x)
    (by norm_num [x1,QInput.real] at hx; linarith [hx.2] : x ≤ 1)
    (by norm_num [p,QInput.real] : 0 < p.real) (by norm_num [p,QInput.real] : p.real < 1)
    (by norm_num [y,QInput.real] : 0 < y.real) (by norm_num [y,QInput.real] : y.real < 1)
    (by rw [parameters_match.2.1]; exact W043.parameters.2.2.2.2.1)
    (by rw [parameters_match.2.2]; exact W043.parameters.2.2.1)
    (by rw [parameters_match.2.2]; exact W043.parameters.2.2.2.1.le)
    (by norm_num : (0:ℝ) < 1/847) hf.1 hf.2.1 hf.2.2 collatz (exponent_bound hα hx)
  refine ⟨supportLine,supportLine_mem,?_⟩
  norm_num [supportLine,m,c,QInput.real,W043.w,W043.v,Fix.sc,scale] at hrate ⊢
  exact hrate

theorem certified_uniform : CertifiedBox 235 (4/10000000) 7000000 :=
  certified.mono_constant (by norm_num)

#print axioms certified
#print axioms certified_uniform
end Spin.Structured.DenseOccupationFixed.B093
