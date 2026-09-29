import SpinCodes.Structured.DenseOccupationGeometry
import SpinCodes.Structured.DenseOccupationGeometryData

namespace Spin.Structured.DenseGeometry

/-- Exact global geometric coverage, independently of which numerical family certifies each box. -/
theorem global_cover {α x : ℝ}
    (hα : α ∈ Set.Icc (1/10000) 1) (hx : x ∈ Set.Icc (13/125) (112/125)) :
    ∃ i, i < 1023 ∧ (boxes.getD i zeroRect).Contains α x := by
  have hp : whole.Contains α x := by
    norm_num [Rect.Contains,whole]
    exact ⟨hα.1,hα.2,hx.1,hx.2⟩
  simpa only [boxes_length] using check_sound boxes tree checked hp

/-- A pointwise certificate for each indexed box discharges the entire dense rectangle. -/
theorem of_box_certificates (P : ℝ → ℝ → Prop)
    (h : ∀ i, i < 1023 → ∀ α x, (boxes.getD i zeroRect).Contains α x → P α x)
    {α x : ℝ} (hα : α ∈ Set.Icc (1/10000) 1)
    (hx : x ∈ Set.Icc (13/125) (112/125)) : P α x := by
  obtain ⟨i,hi,hp⟩ := global_cover hα hx
  exact h i hi α x hp

#print axioms global_cover
#print axioms of_box_certificates
end Spin.Structured.DenseGeometry
