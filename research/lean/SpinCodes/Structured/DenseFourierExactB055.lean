import SpinCodes.Structured.DenseFourierExactB055Data
import SpinCodes.Structured.DenseFourierExactW018
import SpinCodes.Structured.DenseOccupationFixedVertex
import SpinCodes.Structured.DenseFourierExactRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseFourierExact.B055
open Spin.Numeric Spin.Imt DenseOccupationFixed DenseGeometry Set
lemma vertex00 : boxExponent m.real c.real p.real y.real radius.real z.real a0.real x0.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a0 x0 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_00
  norm_num [upper,scale] at h ⊢
  exact h
lemma vertex01 : boxExponent m.real c.real p.real y.real radius.real z.real a0.real x1.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a0 x1 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_01
  norm_num [upper,scale] at h ⊢
  exact h
lemma vertex10 : boxExponent m.real c.real p.real y.real radius.real z.real a1.real x0.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a1 x0 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_10
  norm_num [upper,scale] at h ⊢
  exact h
lemma vertex11 : boxExponent m.real c.real p.real y.real radius.real z.real a1.real x1.real ≤
    -(4/10000000:ℝ) := by
  have h := vertexCheck_sound m c p y radius z a1 x1 (by decide) (by decide) (by decide) (by decide)
    (by decide) (by decide) (by decide) (by decide) check_11
  norm_num [upper,scale] at h ⊢
  exact h

/-- Every point in the Fourier box 55 has the stated negative exponent. -/
theorem exponent_bound {α x : ℝ} (hα : α∈Icc a0.real a1.real) (hx : x∈Icc x0.real x1.real) :
    boxExponent m.real c.real p.real y.real radius.real z.real α x ≤ -(4/10000000:ℝ) := by
  apply boxExponent_le_of_vertices m.real c.real p.real y.real radius.real z.real _
    (by norm_num [p,QInput.real]) (by norm_num [p,QInput.real])
    (by norm_num [y,QInput.real]) (by norm_num [y,QInput.real])
    (by norm_num [a0,QInput.real]) (by norm_num [a1,QInput.real])
    (by norm_num [x0,QInput.real]) (by norm_num [x1,QInput.real]) hα hx
    vertex00 vertex01 vertex10 vertex11


theorem parameters_match : p.real*y.real = (W018.qn:ℝ)/W018.qd ∧ radius.real = (W018.rn:ℝ)/W018.rd ∧ z.real = (W018.zn:ℝ)/W018.zd := by
  norm_num [p,y,radius,z,QInput.real,W018.qn,W018.qd,W018.rn,W018.rd,W018.zn,W018.zd]
theorem collatz : ((ConcreteFourier.matrix (p.real*y.real) z.real).applyCol W018.v.real).le (Coords.smul radius.real W018.v.real) := by
  rw [parameters_match.1,parameters_match.2.1,parameters_match.2.2]
  exact W018.collatz
def supportLine : ℚ×ℚ := ((m.num:ℚ)/m.den,(c.num:ℚ)/c.den)
lemma supportLine_mem : supportLine∈Spin.Majorant.refined.supports := by decide +kernel
lemma geometry_match : boxes.getD 224 zeroRect = ⟨(a0.num:ℚ)/a0.den,(a1.num:ℚ)/a1.den,(x0.num:ℚ)/x0.den,(x1.num:ℚ)/x1.den⟩ := by decide +kernel
theorem certified : CertifiedBox 224 (4/10000000) 206550720 := by
  intro α x hg
  rw [geometry_match] at hg
  have hα : α∈Icc a0.real a1.real := by
    norm_num [Rect.Contains,a0,a1,x0,x1,QInput.real] at hg ⊢
    exact ⟨hg.1,hg.2.1⟩
  have hx : x∈Icc x0.real x1.real := by
    norm_num [Rect.Contains,a0,a1,x0,x1,QInput.real] at hg ⊢
    exact ⟨hg.2.2.1,hg.2.2.2⟩
  have hf := W018.witness_floor
  have hr : PointRate α x m.real c.real (4/10000000) (W018.v.real.Z/W018.wmin) := fourier_pointRate
    (by norm_num [a0,QInput.real] at hα; linarith [hα.1] : 0<α)
    (by norm_num [a1,QInput.real] at hα; linarith [hα.2] : α≤1)
    (by norm_num [x0,QInput.real] at hx; linarith [hx.1] : 0<x)
    (by norm_num [x1,QInput.real] at hx; linarith [hx.2] : x≤1)
    (by norm_num [p,QInput.real] : 0<p.real) (by norm_num [p,QInput.real] : p.real<1)
    (by norm_num [y,QInput.real] : 0<y.real) (by norm_num [y,QInput.real] : y.real<1)
    (by norm_num [radius,QInput.real] : 0<radius.real)
    (by norm_num [z,QInput.real] : 0<z.real) (by norm_num [z,QInput.real] : z.real≤1)
    hf.1 hf.2.1 hf.2.2.1 hf.2.2.2 collatz (exponent_bound hα hx)
  have hrFinal : PointRate α x m.real c.real (4/10000000) W018.prefactor := hr.mono_constant W018.witness_prefactor
  refine ⟨supportLine,supportLine_mem,?_⟩
  simpa only [PointRate,supportLine,QInput.real,Rat.cast_div,Rat.cast_intCast,W018.prefactor] using @hrFinal
theorem certified_uniform : CertifiedBox 224 (4/10000000) 24000000000000 := certified.mono_constant (by norm_num)
#print axioms exponent_bound
#print axioms collatz
#print axioms certified
#print axioms certified_uniform
end Spin.Structured.DenseFourierExact.B055
