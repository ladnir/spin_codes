import SpinCodes.Structured.DenseScalarExactB227
import SpinCodes.Structured.DenseScalarPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseScalarExact.B227
open DenseOccupationFixed DenseGeometry
set_option maxRecDepth 100000

/-- The local scalar witness certifies its original indexed box and actual outer support. -/
theorem certified : CertifiedBox 722 (4/10000000) 1 := by
  intro α x hp
  have he : boxes.getD 722 zeroRect = ⟨(20470001/20480000),(40950001/40960000),(87621922880410294620235849331902475503/170141183460469231731687303715884105728),(12967792439091/25000000000000)⟩ := by rfl
  rw [he] at hp
  norm_num [Rect.Contains] at hp
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [QInput.real,a0,a1,Set.mem_Icc]
    exact And.intro hp.1 hp.2.1
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [QInput.real,x0,x1,Set.mem_Icc]
    exact And.intro hp.2.2.1 hp.2.2.2
  refine ⟨((-7/100),(1625681846115929475145509125033299608399/4253529586511730793292182592897102643200)), Spin.Majorant.member23, ?_⟩
  have hm : (((-7/100):ℚ):ℝ) = m.real := by norm_num [m,QInput.real]
  have hc : (((1625681846115929475145509125033299608399/4253529586511730793292182592897102643200):ℚ):ℝ) = c.real := by norm_num [c,QInput.real]
  change PointRate α x (((-7/100):ℚ):ℝ) (((1625681846115929475145509125033299608399/4253529586511730793292182592897102643200):ℚ):ℝ) (4/10000000) 1
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
end Spin.Structured.DenseScalarExact.B227
