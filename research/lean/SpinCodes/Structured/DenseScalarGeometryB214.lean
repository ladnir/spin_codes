import SpinCodes.Structured.DenseScalarExactB214
import SpinCodes.Structured.DenseScalarPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseScalarExact.B214
open DenseOccupationFixed DenseGeometry
set_option maxRecDepth 100000

/-- The local scalar witness certifies its original indexed box and actual outer support. -/
theorem certified : CertifiedBox 709 (4/10000000) 1 := by
  intro α x hp
  have he : boxes.getD 709 zeroRect = ⟨(630001/640000),(1270001/1280000),(253319723456283/500000000000000),(43385881287743932506819717224474801189/85070591730234615865843651857942052864)⟩ := by rfl
  rw [he] at hp
  norm_num [Rect.Contains] at hp
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [QInput.real,a0,a1,Set.mem_Icc]
    exact And.intro hp.1 hp.2.1
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [QInput.real,x0,x1,Set.mem_Icc]
    exact And.intro hp.2.2.1 hp.2.2.2
  refine ⟨((-3/100),(3076970006775960790657142966285747138917/8507059173023461586584365185794205286400)), Spin.Majorant.member21, ?_⟩
  have hm : (((-3/100):ℚ):ℝ) = m.real := by norm_num [m,QInput.real]
  have hc : (((3076970006775960790657142966285747138917/8507059173023461586584365185794205286400):ℚ):ℝ) = c.real := by norm_num [c,QInput.real]
  change PointRate α x (((-3/100):ℚ):ℝ) (((3076970006775960790657142966285747138917/8507059173023461586584365185794205286400):ℚ):ℝ) (4/10000000) 1
  rw [hm,hc]
  exact pointRate
    (lt_of_lt_of_le (by norm_num [a0,QInput.real]) hα.1)
    (le_trans hα.2 (by norm_num [a1,QInput.real]))
    (lt_of_lt_of_le (by norm_num [x0,QInput.real]) hx.1)
    (le_trans hx.2 (by norm_num [x1,QInput.real]))
    (by norm_num [p,QInput.real]) (by norm_num [p,QInput.real])
    (by norm_num [y,QInput.real]) (by norm_num [y,QInput.real])
    (by norm_num [radius,QInput.real]) (by norm_num [z,QInput.real])
    (by norm_num [z,QInput.real]) scalar_bound (exponent_bound hα hx)

#print axioms certified
end Spin.Structured.DenseScalarExact.B214
