import SpinCodes.Structured.DenseScalarExactB050
import SpinCodes.Structured.DenseScalarPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseScalarExact.B050
open DenseOccupationFixed DenseGeometry
set_option maxRecDepth 100000

/-- The local scalar witness certifies its original indexed box and actual outer support. -/
theorem certified : CertifiedBox 52 (4/10000000) 1 := by
  intro α x hp
  have he : boxes.getD 52 zeroRect = ⟨(327670001/327680000),(655350001/655360000),(129223289418938324344194200630890939964085310021/259614842926741381426524816461004800000000000000),(1/2)⟩ := by rfl
  rw [he] at hp
  norm_num [Rect.Contains] at hp
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [QInput.real,a0,a1,Set.mem_Icc]
    exact And.intro hp.1 hp.2.1
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [QInput.real,x0,x1,Set.mem_Icc]
    exact And.intro hp.2.2.1 hp.2.2.2
  refine ⟨((0/1),(3465735902799727/10000000000000000)), Spin.Majorant.member19, ?_⟩
  have hm : (((0/1):ℚ):ℝ) = m.real := by norm_num [m,QInput.real]
  have hc : (((3465735902799727/10000000000000000):ℚ):ℝ) = c.real := by norm_num [c,QInput.real]
  change PointRate α x (((0/1):ℚ):ℝ) (((3465735902799727/10000000000000000):ℚ):ℝ) (4/10000000) 1
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
end Spin.Structured.DenseScalarExact.B050
