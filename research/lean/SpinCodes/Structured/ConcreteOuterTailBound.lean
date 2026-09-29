import SpinCodes.Structured.ConcreteOuterTailDenseCount
import SpinCodes.Structured.ConcreteOuterTailSparseLimit
import SpinCodes.Structured.ConcreteNativeSelection

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric
open ConcreteNativeFamily

def outerTail (k : ℕ) : ℝ :=
  ∑ w ∈ badWeights (k*24) (weightWindow (k*24)),
    (seedLaw k).expect (fun seed => (spectrum seed w:ℝ))

lemma expected_spectrum_split (k w : ℕ) :
    (seedLaw k).expect (fun seed => (spectrum seed w:ℝ)) =
      (∑ x ∈ sparseMessages k, (seedLaw k).prob (fun seed => wtF (encode seed x)=w))+denseShell k w := by
  simp_rw [spectrum_eq_sum]
  rw [FinPMF.expect_sum]
  simp_rw [←FinPMF.prob_eq_expect_indicator]
  unfold denseShell sparseMessages denseMessages
  simp only [sum_filter]
  rw [←sum_add_distrib]
  apply sum_congr rfl
  intro x hx
  by_cases hx0 : x=0 <;> by_cases hsp : 1000*golayOccupation x ≤ k*24 <;>
    simp [hx0,hsp,Nat.lt_of_not_ge,hsp,show (k*24<1000*golayOccupation x) ↔ ¬1000*golayOccupation x ≤ k*24 from lt_iff_not_ge]

lemma bad_weight_tail {b w : ℕ} (hw : w∈badWeights b (weightWindow b)) :
    w ≤ b ∧ (125*w ≤ 13*b ∨ 112*b ≤ 125*w) := by
  obtain ⟨hwb,hnot⟩ := mem_sdiff.mp hw
  have hwb' : w ≤ b := by have := mem_range.mp hwb; omega
  have hh : ¬(13*b ≤ 125*w ∧ 125*w ≤ 112*b) := by
    intro hh
    exact hnot (mem_filter.mpr ⟨hwb,hh⟩)
  exact ⟨hwb',by omega⟩

lemma sum_bad_le_tails {b : ℕ} (hb : 0<b) (f : ℕ→ℝ) (hf : ∀ w,0 ≤ f w) :
    (∑ w ∈ badWeights b (weightWindow b), f w) ≤
      (∑ w ∈ range (b*13/125+1),f w)+(∑ w ∈ Icc (b-b*13/125) b,f w) := by
  have hsub : badWeights b (weightWindow b) ⊆ range (b*13/125+1) ∪ Icc (b-b*13/125) b := by
    intro w hw
    obtain ⟨hwb,hlo|hhi⟩ := bad_weight_tail hw
    · apply mem_union_left
      apply mem_range.mpr
      omega
    · apply mem_union_right
      apply mem_Icc.mpr
      constructor <;> omega
  have hdis : Disjoint (range (b*13/125+1)) (Icc (b-b*13/125) b) := by
    apply disjoint_left.mpr
    intro w hw₁ hw₂
    have h₁ := mem_range.mp hw₁
    have h₂ := (mem_Icc.mp hw₂).1
    omega
  calc
    (∑ w ∈ badWeights b (weightWindow b),f w) ≤
        ∑ w ∈ range (b*13/125+1) ∪ Icc (b-b*13/125) b,f w :=
      sum_le_sum_of_subset_of_nonneg hsub (fun w _ _ => hf w)
    _ = _ := sum_union hdis

lemma sparse_bad_tail_le {k : ℕ} (hk : 0<k) :
    (∑ w ∈ badWeights (k*24) (weightWindow (k*24)),
      ∑ x ∈ sparseMessages k,(seedLaw k).prob (fun seed => wtF (encode seed x)=w)) ≤
      sparseLowerTail k+sparseUpperTail k := by
  rw [sum_comm]
  unfold sparseLowerTail sparseUpperTail
  rw [←sum_add_distrib]
  exact sum_le_sum (fun x _ => sum_bad_le_tails (by omega) _ (fun w => FinPMF.prob_nonneg _ _))

lemma dense_bad_tail_le {k : ℕ} (hk : 0<k) :
    (∑ w ∈ badWeights (k*24) (weightWindow (k*24)),denseShell k w) ≤
      (((k*24:ℕ):ℝ)+1)^5*Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) := by
  have hh := sum_le_sum (fun w (hw : w∈badWeights (k*24) (weightWindow (k*24))) =>
    denseShell_le hk w (bad_weight_tail hw).2)
  refine hh.trans ?_
  rw [sum_const,nsmul_eq_mul]
  have hc : (badWeights (k*24) (weightWindow (k*24))).card ≤ k*24+1 := by
    exact (card_le_card sdiff_subset).trans_eq (card_range _)
  have hcR : ((badWeights (k*24) (weightWindow (k*24))).card:ℝ) ≤ ((k*24:ℕ):ℝ)+1 := by exact_mod_cast hc
  have hh := mul_le_mul_of_nonneg_right hcR
    (show 0 ≤ (((k*24:ℕ):ℝ)+1)^4*Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) by positivity)
  convert hh using 1 <;> ring

/-- Finite expected number of outer words outside the window, with every term explicit. -/
theorem outerTail_le {k : ℕ} (hk : 0<k) :
    outerTail k ≤ sparseLowerTail k+sparseUpperTail k+
      (((k*24:ℕ):ℝ)+1)^5*Real.exp (-(768/10^10)*((k*24:ℕ):ℝ)+denseLogError (k*24)) := by
  unfold outerTail
  simp_rw [expected_spectrum_split]
  rw [sum_add_distrib]
  exact add_le_add (sparse_bad_tail_le hk) (dense_bad_tail_le hk)

end Spin.Structured.ConcreteOuter

