import SpinCodes.Structured.DenseScalarExactB183
import SpinCodes.Structured.DenseScalarPointRate
import SpinCodes.Structured.DenseGeometryRate

noncomputable section
namespace Spin.Structured.DenseScalarExact.B183
open DenseOccupationFixed DenseGeometry
set_option maxRecDepth 100000

/-- The local scalar witness certifies its original indexed box and actual outer support. -/
theorem certified : CertifiedBox 636 (4/10000000) 1 := by
  intro α x hp
  have he : boxes.getD 636 zeroRect = ⟨(10230001/10240000),(20470001/20480000),(84219921256861908840067417679371902625/170141183460469231731687303715884105728),(129223289418938324344194200630890939964085310021/259614842926741381426524816461004800000000000000)⟩ := by rfl
  rw [he] at hp
  norm_num [Rect.Contains] at hp
  have hα : α ∈ Set.Icc a0.real a1.real := by
    norm_num [QInput.real,a0,a1,Set.mem_Icc]
    exact And.intro hp.1 hp.2.1
  have hx : x ∈ Set.Icc x0.real x1.real := by
    norm_num [QInput.real,x0,x1,Set.mem_Icc]
    exact And.intro hp.2.2.1 hp.2.2.2
  refine ⟨((1/100),(58119563056842377037993588567825857659/170141183460469231731687303715884105728)), Spin.Majorant.member18, ?_⟩
  have hm : (((1/100):ℚ):ℝ) = m.real := by norm_num [m,QInput.real]
  have hc : (((58119563056842377037993588567825857659/170141183460469231731687303715884105728):ℚ):ℝ) = c.real := by norm_num [c,QInput.real]
  change PointRate α x (((1/100):ℚ):ℝ) (((58119563056842377037993588567825857659/170141183460469231731687303715884105728):ℚ):ℝ) (4/10000000) 1
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
end Spin.Structured.DenseScalarExact.B183
