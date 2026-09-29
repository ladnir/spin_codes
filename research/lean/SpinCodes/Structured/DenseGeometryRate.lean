import SpinCodes.Structured.DenseOccupationGeometryCover
import SpinCodes.Structured.DenseOccupationFamilyRateDefs

noncomputable section
namespace Spin.Structured.DenseGeometry
open DenseOccupationFixed

/-- A numerical family certifies an indexed rectangle with an actual outer support. -/
def CertifiedBox (i : ℕ) (η C : ℝ) : Prop :=
  ∀ α x, (boxes.getD i zeroRect).Contains α x →
    ∃ p : ℚ × ℚ, p ∈ Spin.Majorant.refined.supports ∧ PointRate α x p.1 p.2 η C

theorem CertifiedBox.mono_constant {i : ℕ} {η C C' : ℝ}
    (h : CertifiedBox i η C) (hC : C ≤ C') : CertifiedBox i η C' := by
  intro α x hp
  obtain ⟨p,hp,hr⟩ := h α x hp
  exact ⟨p,hp,hr.mono_constant hC⟩

/-- The independently checked geometry combines all indexed numerical families. -/
theorem denseRates_of_boxes {η C : ℝ} (h : ∀ i, i<1023 → CertifiedBox i η C) :
    DenseRates η C := by
  intro α x hα hx
  obtain ⟨i,hi,hp⟩ := global_cover hα hx
  exact h i hi α x hp

#print axioms CertifiedBox.mono_constant
#print axioms denseRates_of_boxes
end Spin.Structured.DenseGeometry
