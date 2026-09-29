import SpinCodes.Structured.ConcreteOuterTailDenseShift
import SpinCodes.Structured.ConcreteOuterTailDenseReal
import SpinCodes.Structured.ConcreteOuterMajorantBridge

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Spin.Numeric

/-- The certified gap tolerates the half-bit endpoint displacement forced by parity. -/
theorem dense_near_gap {B a c u : ℝ} (hB : 1 ≤ B)
    (ha : B/125 ≤ a) (haB : a ≤ B)
    (hf₁ : a/2 ≤ c) (hf₂ : c ≤ B-a/2)
    (hf₃ : c/2 ≤ u) (hf₄ : u ≤ B-c/2)
    (hu : u ≤ 13*B/125+1/2) :
    B*(gBA (a/B)+piBA (a/B) (c/B)+piBA (c/B) (u/B)) ≤
      -(768/10^10)*B+36*(Real.log B+1)+250 := by
  have hB0 : 0 < B := by linarith
  let u' := min u (13*B/125)
  let c' := min c (2*u')
  let a' := min a (2*c')
  have ha0 : 0 ≤ a := by linarith
  have hc0 : 0 ≤ c := by linarith
  have hu0 : B/500 ≤ u := by linarith
  have hu'lo : B/500 ≤ u' := le_min hu0 (by linarith)
  have hu'u : u' ≤ u := min_le_left _ _
  have hu'hi : u' ≤ 13*B/125 := min_le_right _ _
  have hud : u-u' ≤ 1/2 := by
    dsimp [u']; rcases le_total u (13*B/125) with h|h
    · rw [min_eq_left h]; linarith
    · rw [min_eq_right h]; linarith
  have hc'lo : B/250 ≤ c' := le_min (by linarith) (by linarith)
  have hc'c : c' ≤ c := min_le_left _ _
  have hc'u : c' ≤ 2*u' := min_le_right _ _
  have hcd : c-c' ≤ 1 := by
    dsimp [c']; rcases le_total c (2*u') with h|h
    · rw [min_eq_left h]; linarith
    · rw [min_eq_right h]; linarith
  have ha'lo : B/125 ≤ a' := le_min ha (by linarith)
  have ha'a : a' ≤ a := min_le_left _ _
  have ha'c : a' ≤ 2*c' := min_le_right _ _
  have had : a-a' ≤ 2 := by
    dsimp [a']; rcases le_total a (2*c') with h|h
    · rw [min_eq_left h]; linarith
    · rw [min_eq_right h]; linarith
  have ha'0 : 0 ≤ a' := by linarith
  have hc'0 : 0 ≤ c' := by linarith
  have hf₁' : a'/2 ≤ c' := by linarith
  have hf₂' : c' ≤ B-a'/2 := by linarith
  have hf₃' : c'/2 ≤ u' := by linarith
  have hf₄' : u' ≤ B-c'/2 := by linarith
  have ha'B : a' ≤ B := ha'a.trans haB
  have hcc : |c-c'| ≤ 1 := abs_le.mpr ⟨by linarith, hcd⟩
  have haa : |a-a'| ≤ 2 := abs_le.mpr ⟨by linarith, had⟩
  have huu : |u-u'| ≤ 1 := abs_le.mpr ⟨by linarith, by linarith⟩
  have hg := gBA_count_step hB0 ha haB ha'0 ha'a had
  have hπ₁ := piBA_count_step hB ha0 hf₁ hf₂ ha'0 hf₁' hf₂' haa hcc
  have hπ₂ := piBA_count_step hB hc0 hf₃ hf₄ hc'0 hf₃' hf₄' (by linarith : |c-c'| ≤ 2) huu
  have hgap := denseTail_closed_real
    (show a'/B ∈ Set.Icc ((1:ℝ)/125) 1 from
      ⟨(le_div_iff₀ hB0).mpr (by linarith), (div_le_one hB0).mpr ha'B⟩)
    (show c'/B ∈ Set.Icc (0:ℝ) 1 from
      ⟨by positivity, (div_le_one hB0).mpr (by linarith)⟩)
    (show u'/B ∈ Set.Icc ((1:ℝ)/500) (13/125) from
      ⟨(le_div_iff₀ hB0).mpr (by linarith), (div_le_iff₀ hB0).mpr (by linarith)⟩)
    ⟨(normalized_feasible hB0 hf₁' hf₂').1, (normalized_feasible hB0 hf₁' hf₂').2,
      (normalized_feasible hB0 hf₃' hf₄').1, (normalized_feasible hB0 hf₃' hf₄').2,
      by positivity, (div_le_one hB0).mpr ha'B⟩
  have hgapB := mul_le_mul_of_nonneg_left hgap.le hB0.le
  linarith [(abs_le.mp hπ₁).2, (abs_le.mp hπ₂).2]

end Spin.Structured.ConcreteOuter
