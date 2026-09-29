import SpinCodes.Structured.ConcretePlacementLimitRiemann

/-! The actual omitted impulse durations give exactly an a/R product perturbation. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder FiniteKernel

theorem normalizedSites_gapRight {R a : Nat} (hR : 0 < R) (B : BlockSubset R a) (i : Fin (a + 1)) :
    (Fin.snoc (normalizedSites B) (1 : ℝ) : Fin (a + 1) → ℝ) i = (gapRight B i : ℝ) / R := by
  refine Fin.lastCases ?_ (fun j => ?_) i
  · simp only [Fin.snoc_last, gapRight]
    rw [div_self (by exact_mod_cast Nat.ne_of_gt hR)]
  · simp only [Fin.snoc_castSucc, gapRight, normalizedSites]

theorem positionGaps_normalizedSites {R a : Nat} (hR : 0 < R) (B : BlockSubset R a) (i : Fin (a + 1)) :
    positionGaps (normalizedSites B) i = (emptyGaps B i : ℝ) / R + if i = 0 then 0 else 1 / (R : ℝ) := by
  unfold positionGaps
  rw [normalizedSites_gapRight hR]
  have he : (emptyGaps B i : ℝ) = (gapRight B i : ℝ) - (gapLeft B i : ℝ) :=
    Nat.cast_sub (gapLeft_le_right B i)
  rw [he]
  refine Fin.cases ?_ (fun j => ?_) i
  · simp [gapLeft]
  · simp only [Fin.cons_succ, Fin.succ_ne_zero, ite_false, gapLeft, Nat.cast_add, Nat.cast_one,
      normalizedSites]
    ring

theorem positionGaps_normalizedSites_nonneg {R a : Nat} (hR : 0 < R) (B : BlockSubset R a)
    (i : Fin (a + 1)) : 0 ≤ positionGaps (normalizedSites B) i := by
  rw [positionGaps_normalizedSites hR]
  split_ifs <;> positivity

theorem positionGaps_normalizedSites_distance {R a : Nat} (hR : 0 < R) (B : BlockSubset R a) :
    ∑ i, |(emptyGaps B i : ℝ) / R - positionGaps (normalizedSites B) i| = (a : ℝ) / R := by
  simp only [positionGaps_normalizedSites hR, sub_add_cancel_left, abs_neg]
  rw [Fin.sum_univ_succ]
  simp only [ite_true, abs_zero, Fin.succ_ne_zero, ite_false,
    abs_of_nonneg (by positivity : (0 : ℝ) ≤ 1 / R), sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul, zero_add]
  ring

theorem placementProduct_site_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ) (B : BlockSubset R a) :
    RowError (placementProduct θ B) (siteProduct (θ * epochMean / 128) (normalizedSites B))
      ((θ * epochMean / 128) * ((a : ℝ) / R)) := by
  rw [placementProduct_normalized hR]
  have hγ : 0 ≤ θ * epochMean / 128 := by unfold epochMean; positivity
  have h := simplexProduct_rowError hγ (fun i => (emptyGaps B i : ℝ) / R)
    (positionGaps (normalizedSites B)) (fun i => by positivity) (positionGaps_normalizedSites_nonneg hR B)
  rw [positionGaps_normalizedSites_distance hR B] at h
  exact h

def stateClass (q : State) : Fin 2 := if q = ∅ then 0 else 1

theorem lifted_entry (K : Matrix (Fin 2) (Fin 2) ℝ) (q r : State) :
    (liveLift * K * liveProjection) q r =
      if r = ∅ then K (stateClass q) 0 else K (stateClass q) 1 / 524287 := by
  by_cases hq : q = ∅ <;> by_cases hr : r = ∅ <;>
    simp [Matrix.mul_apply, Fin.sum_univ_two, liveLift, liveProjection, coarseValue, stateClass, hq, hr, div_eq_mul_inv]

theorem rowError_entry_general {α : Type*} [Fintype α] [DecidableEq α]
    {K L : Matrix α α ℝ} {ε : ℝ} (h : RowError K L ε) (q r : α) : |K q r - L q r| ≤ ε :=
  (Finset.single_le_sum (fun s _ => abs_nonneg ((K - L) q s)) (mem_univ r)).trans (h q)

theorem lifted_rowError_entry {K L : Matrix (Fin 2) (Fin 2) ℝ} {ε : ℝ} (hε : 0 ≤ ε)
    (h : RowError K L ε) (q r : State) :
    |(liveLift * K * liveProjection) q r - (liveLift * L * liveProjection) q r| ≤ ε := by
  simp only [lifted_entry]
  by_cases hr : r = ∅
  · simp only [hr, ite_true]
    exact rowError_entry_general h _ _
  · simp only [hr, ite_false, ← sub_div, abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 524287)]
    have hh := div_le_div_of_nonneg_right (rowError_entry_general h (stateClass q) 1) (by norm_num : (0 : ℝ) ≤ 524287)
    exact hh.trans (div_le_self hε (by norm_num))

end Spin.Structured.Placement

