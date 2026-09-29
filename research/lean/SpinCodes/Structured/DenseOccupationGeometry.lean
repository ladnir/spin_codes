import SpinCodes.Structured.DenseOccupationGeometryDefs
import Mathlib.Tactic

namespace Spin.Structured.DenseGeometry

def Rect.Contains (r : Rect) (α x : ℝ) : Prop :=
  (r.alo:ℝ) ≤ α ∧ α ≤ (r.ahi:ℝ) ∧ (r.xlo:ℝ) ≤ x ∧ x ≤ (r.xhi:ℝ)

/-- Kernel-checked subdivisions cover every real point of the requested rectangle. -/
theorem check_sound (boxes : List Rect) (t : Tree) {r : Rect}
    (h : check boxes t r = true) {α x : ℝ} (hp : r.Contains α x) :
    ∃ i, i < boxes.length ∧ (boxes.getD i zeroRect).Contains α x := by
  induction t generalizing r with
  | leaf i =>
    have hh := of_decide_eq_true h
    refine ⟨i,hh.1,?_,?_,?_,?_⟩
    · exact le_trans (by exact_mod_cast hh.2.1) hp.1
    · exact le_trans hp.2.1 (by exact_mod_cast hh.2.2.1)
    · exact le_trans (by exact_mod_cast hh.2.2.2.1) hp.2.2.1
    · exact le_trans hp.2.2.2 (by exact_mod_cast hh.2.2.2.2)
  | alpha c l u hl hu =>
    have hh := h
    simp only [check, Bool.and_eq_true] at hh
    by_cases ha : α ≤ (c:ℝ)
    · exact hl hh.1 ⟨hp.1,ha,hp.2.2.1,hp.2.2.2⟩
    · exact hu hh.2 ⟨(not_le.mp ha).le,hp.2.1,hp.2.2.1,hp.2.2.2⟩
  | density c l u hl hu =>
    have hh := h
    simp only [check, Bool.and_eq_true] at hh
    by_cases hx : x ≤ (c:ℝ)
    · exact hl hh.1 ⟨hp.1,hp.2.1,hp.2.2.1,hx⟩
    · exact hu hh.2 ⟨hp.1,hp.2.1,(not_le.mp hx).le,hp.2.2.2⟩

#print axioms check_sound
end Spin.Structured.DenseGeometry
