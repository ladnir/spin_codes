import SpinCodes.Structured.DenseOccupationSound

noncomputable section
namespace Spin.Structured.DenseOccupation
open Spin.Imt

/-- A successful rational check yields the precise concrete matrix hypothesis. -/
theorem check_sound {q z radius floor : ℚ} {v : QCoords}
    (h : check q z radius floor v = true) :
    (0:ℝ) ≤ q ∧ (q:ℝ) ≤ 1 ∧ (0:ℝ) < z ∧ (z:ℝ) ≤ 1 ∧
    (0:ℝ) < floor ∧ (floor:ℝ) ≤ v.real.Z ∧ (floor:ℝ) ≤ v.real.D ∧
    (∀ i, (floor:ℝ) ≤ v.real.S i) ∧
    ((Occupation.Sparse.numericalMatrix q z).applyCol v.real).le
      (Coords.smul radius v.real) := by
  simp only [check, decide_eq_true_eq] at h
  obtain ⟨hq0,hq1,hz0,hz1,hf,hvZ,hvD,hvS,hZ,hD,hS⟩ := h
  refine ⟨by exact_mod_cast hq0, by exact_mod_cast hq1,
    by exact_mod_cast hz0, by exact_mod_cast hz1, by exact_mod_cast hf,
    by change (floor:ℝ) ≤ (v.Z:ℝ); exact_mod_cast hvZ, by change (floor:ℝ) ≤ (v.D:ℝ); exact_mod_cast hvD,
    fun i => by change (floor:ℝ) ≤ (v.S i:ℝ); exact_mod_cast hvS i, ?_⟩
  rw [← cast_column]
  refine ⟨?_,?_,fun i => ?_⟩
  · change ((column q z v).Z:ℝ) ≤ (radius:ℝ)*(v.Z:ℝ)
    exact_mod_cast hZ
  · change ((column q z v).D:ℝ) ≤ (radius:ℝ)*(v.D:ℝ)
    exact_mod_cast hD
  · change ((column q z v).S i:ℝ) ≤ (radius:ℝ)*(v.S i:ℝ)
    exact_mod_cast hS i

/-- A checked witness contracts every concrete finite matrix iterate. -/
theorem checked_iterate {q z radius floor : ℚ} {v : QCoords}
    (h : check q z radius floor v = true) (R : ℕ) :
    (((Occupation.Sparse.numericalMatrix q z).apply)^[R] (Coords.eZ 5)).total ≤
      (radius:ℝ)^R*(v.Z:ℝ)/(floor:ℝ) := by
  obtain ⟨hq0,hq1,hz0,hz1,hf,hvZ,hvD,hvS,hcol⟩ := check_sound h
  exact PositiveFinite.occupation_iterate_le hq0 hq1 hz0.le hf hvZ hvD hvS hcol R

end Spin.Structured.DenseOccupation
