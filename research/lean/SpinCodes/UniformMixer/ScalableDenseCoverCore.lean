import SpinCodes.Structured.DenseOccupationGeometry
import SpinCodes.UniformMixer.ScalableDenseRateDefs

/-! Chunked geometric assembly for the revised dense certificate. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.DenseGeometry

/-- A property holds throughout one exact rational rectangle. -/
def DenseClaim (P : ℝ → ℝ → Prop) (r : Rect) : Prop :=
  ∀ α x : ℝ, r.Contains α x → P α x

/-- The pointwise target shared by every generated dense chunk. -/
def DenseTarget (η C : ℝ) (α x : ℝ) : Prop :=
  ∃ p : ℚ × ℚ, p ∈ Spin.Majorant.refined.supports ∧
    ScalableDensePointRate α x p.1 p.2 η C

/-- One checked tree and checked leaves establish its rectangle claim. -/
theorem denseClaim_of_checked_tree (P : ℝ → ℝ → Prop)
    (boxes : List Rect) (tree : Tree) (r : Rect)
    (hcheck : check boxes tree r = true)
    (hleaf : ∀ i, i < boxes.length →
      DenseClaim P (boxes.getD i zeroRect)) :
    DenseClaim P r := by
  intro α x hr
  obtain ⟨i, hi, hbox⟩ := check_sound boxes tree hcheck hr
  exact hleaf i hi α x hbox

/-- Combine two independently compiled claims across an alpha cut. -/
theorem denseClaim_alpha (P : ℝ → ℝ → Prop)
    (r : Rect) (cut : ℚ)
    (left : DenseClaim P {r with ahi := cut})
    (right : DenseClaim P {r with alo := cut}) :
    DenseClaim P r := by
  intro α x hr
  by_cases h : α ≤ (cut : ℝ)
  · exact left α x ⟨hr.1, h, hr.2.2.1, hr.2.2.2⟩
  · exact right α x ⟨(not_le.mp h).le, hr.2.1,
      hr.2.2.1, hr.2.2.2⟩

/-- Combine two independently compiled claims across a density cut. -/
theorem denseClaim_density (P : ℝ → ℝ → Prop)
    (r : Rect) (cut : ℚ)
    (left : DenseClaim P {r with xhi := cut})
    (right : DenseClaim P {r with xlo := cut}) :
    DenseClaim P r := by
  intro α x hr
  by_cases h : x ≤ (cut : ℝ)
  · exact left α x ⟨hr.1, hr.2.1, hr.2.2.1, h⟩
  · exact right α x ⟨hr.1, hr.2.1,
      (not_le.mp h).le, hr.2.2.2⟩

/-- The checked global rectangle is exactly the Scalable dense domain. -/
theorem scalableDenseRates_of_claim (η C : ℝ)
    (h : DenseClaim (DenseTarget η C)
      ⟨1 / 10000, 1, 13 / 125, 112 / 125⟩) :
    ScalableDenseRates η C := by
  intro α x hα hx
  have hr : (⟨1 / 10000, 1, 13 / 125, 112 / 125⟩ : Rect).Contains α x := by
    norm_num [Rect.Contains]
    exact ⟨hα.1, hα.2, hx.1, hx.2⟩
  exact h α x hr

#print axioms denseClaim_of_checked_tree
#print axioms denseClaim_alpha
#print axioms denseClaim_density
#print axioms scalableDenseRates_of_claim

end Spin.UniformMixer
