import SpinCodes.Structured.ConcretePlacementLimitRiemann
import SpinCodes.Structured.ConcreteFixedLargeLaplace

/-! Normalization and elementary laws of the actual ordered-site spacings. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Filter MeasureTheory
open scoped Topology

theorem orderedSite_integral_const_one (Q : ℕ) :
    (∫ _x in orderedSiteDomain Q, (1:ℝ)) = 1/(Q.factorial:ℝ) := by
  have hs := orderedSite_sum_tendsto (a := Q) (fun _ => (1:ℝ)) continuous_const
  simp only [sum_const, card_univ, card_blockSubset, nsmul_eq_mul, mul_one] at hs
  have hf : (Q.factorial:ℝ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero Q
  have ht := (pow_div_choose_tendsto Q).inv₀ hf
  simp only [inv_div, one_div] at ht ⊢
  exact tendsto_nhds_unique hs ht

theorem orderedSite_normalized_one (Q : ℕ) :
    (Q.factorial:ℝ) * (∫ _x in orderedSiteDomain Q, (1:ℝ)) = 1 := by
  rw [orderedSite_integral_const_one]
  exact mul_one_div_cancel (by exact_mod_cast Nat.factorial_ne_zero Q)

theorem positionGaps_sum {Q : ℕ} (x : Fin Q → ℝ) : ∑ i, positionGaps x i = 1 := by
  simp only [positionGaps, sum_sub_distrib, Fin.sum_snoc, Fin.sum_cons]
  ring

theorem positionGaps_nonneg {Q : ℕ} {x : Fin Q → ℝ} (hx : x ∈ orderedSiteDomain Q)
    (i : Fin (Q+1)) : 0 ≤ positionGaps x i := by
  cases Q with
  | zero => fin_cases i; norm_num [positionGaps, Fin.snoc]
  | succ Q =>
    refine Fin.cases ?_ (fun j => ?_) i
    · simpa [positionGaps] using (hx.1 0).1
    · refine Fin.lastCases ?_ (fun k => ?_) j
      · simpa [positionGaps] using sub_nonneg.mpr (hx.1 (Fin.last Q)).2.le
      · simp only [positionGaps, Fin.cons_succ, Fin.succ_castSucc, Fin.snoc_castSucc]
        exact sub_nonneg.mpr (hx.2 (by exact Fin.castSucc_lt_succ)).le

def liveDuration {Q : ℕ} (S : Finset (Fin (Q+1))) (x : Fin Q → ℝ) : ℝ :=
  ∑ i ∈ S, positionGaps x i

theorem liveDuration_mem {Q : ℕ} (S : Finset (Fin (Q+1))) {x : Fin Q → ℝ}
    (hx : x ∈ orderedSiteDomain Q) : liveDuration S x ∈ Set.Icc (0:ℝ) 1 := by
  constructor
  · exact sum_nonneg (fun i _ => positionGaps_nonneg hx i)
  · calc liveDuration S x ≤ ∑ i, positionGaps x i := by
           exact sum_le_sum_of_subset_of_nonneg (subset_univ S) (fun i _ _ => positionGaps_nonneg hx i)
         _ = 1 := positionGaps_sum x

@[simp] theorem liveDuration_empty {Q : ℕ} (x : Fin Q → ℝ) : liveDuration ∅ x = 0 := by
  simp [liveDuration]
@[simp] theorem liveDuration_univ {Q : ℕ} (x : Fin Q → ℝ) : liveDuration univ x = 1 := by
  simpa only [liveDuration] using positionGaps_sum x

theorem continuous_liveDuration {Q : ℕ} (S : Finset (Fin (Q+1))) : Continuous (liveDuration S) := by
  exact continuous_finsetSum _ (fun i _ => continuous_positionGaps i)

theorem spacing_tilt_empty (Q : ℕ) (sig : ℝ) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      1/(1-liveDuration ∅ x+liveDuration ∅ x/sig)^(Q+1)) = 1 := by
  simp only [liveDuration_empty, sub_zero, zero_div, add_zero, one_pow, one_div_one]
  exact orderedSite_normalized_one Q

theorem spacing_tilt_univ (Q : ℕ) (sig : ℝ) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      1/(1-liveDuration univ x+liveDuration univ x/sig)^(Q+1)) = sig^(Q+1) := by
  simp only [liveDuration_univ, sub_self, zero_add, one_div, inv_pow, inv_inv]
  have h : (∫ _x in orderedSiteDomain Q, sig^(Q+1)) = sig^(Q+1) *
      (∫ _x in orderedSiteDomain Q, (1:ℝ)) := by
    rw [← integral_const_mul]
    simp
  rw [h, ← mul_assoc, mul_comm (Q.factorial:ℝ), mul_assoc, orderedSite_normalized_one, mul_one]

#print axioms orderedSite_normalized_one
#print axioms liveDuration_mem
#print axioms spacing_tilt_empty
#print axioms spacing_tilt_univ
end Spin.Structured.Placement
