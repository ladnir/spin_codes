import SpinCodes.Structured.DenseScalarExactB251
import SpinCodes.Structured.DenseScalarPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseScalarExact.B251
open DenseOccupationFixed DenseGeometry
set_option maxRecDepth 100000

/-- The local scalar witness certifies its original indexed box and actual outer support. -/
theorem certified : CertifiedBox 788 (4/10000000) 1 := by
  intro α x hp
  have he : boxes.getD 788 zeroRect = ⟨(1270001/1280000),(2550001/2560000),(54471334987363575166937446516365236680666043211/103845937170696552570609926584401920000000000000),(53744075929833/100000000000000)⟩ := by rfl
  rw [he] at hp
  norm_num [Rect.Contains] at hp
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [QInput.real,a0,a1,Set.mem_Icc]
    exact And.intro hp.1 hp.2.1
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [QInput.real,x0,x1,Set.mem_Icc]
    exact And.intro hp.2.2.1 hp.2.2.2
  refine ⟨((-1/10),(397841147503783/1000000000000000)), Spin.Majorant.member25, ?_⟩
  have hm : (((-1/10):ℚ):ℝ) = m.real := by norm_num [m,QInput.real]
  have hc : (((397841147503783/1000000000000000):ℚ):ℝ) = c.real := by norm_num [c,QInput.real]
  change PointRate α x (((-1/10):ℚ):ℝ) (((397841147503783/1000000000000000):ℚ):ℝ) (4/10000000) 1
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
end Spin.Structured.DenseScalarExact.B251
