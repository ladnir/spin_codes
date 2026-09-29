import SpinCodes.Structured.ConcreteStep

/-! Domination of the actual nonzero-state transition, with the paper's
separate zero, diffuse, and uniform-shell budgets. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open Spin.Imt

private theorem uniform_shells (m : ℝ) (r : Finset (Fin 19)) :
    (∑ i, if r ∈ actualShellSystem.shell i then
      ((((actualShellSystem.shell i).card : ℝ) * m / (2 * 524287)) /
        (actualShellSystem.shell i).card) else 0) =
      if r = ∅ then 0 else m / (2 * 524287) := by
  rw [← actualShellSystem.sum_mem r (m / (2 * 524287))]
  apply Finset.sum_congr rfl
  intro i _
  by_cases hr : r ∈ actualShellSystem.shell i
  · simp only [if_pos hr]
    have hc := (actualShellSystem.card_pos i).ne'
    field_simp
  · simp only [if_neg hr]

theorem liveDominates_branches (lazy fresh : Finset (Fin 19) → ℝ) (m ell u : ℝ)
    (hl : ∀ r, 0 ≤ lazy r) (ht : ∑ r, lazy r ≤ m) (hz : lazy ∅ ≤ ell)
    (hf : ∀ r, fresh r ≤ m) (hfz : fresh ∅ ≤ u) :
    LiveDominates actualShellSystem
      ⟨ell / 2 + u / (2 * 524287), m / 2,
        fun i => ((actualShellSystem.shell i).card : ℝ) * m / (2 * 524287)⟩
      (fun r => lazy r / 2 + fresh r / (2 * 524287)) := by
  classical
  rw [liveDominates_iff_witness]
  refine ⟨fun r => if r = ∅ then 0 else lazy r / 2, ?_, by simp, ?_, ?_⟩
  · intro r
    change 0 ≤ if r = ∅ then 0 else lazy r / 2
    split_ifs
    · exact le_rfl
    · exact div_nonneg (hl r) (by norm_num)
  · change (∑ r, if r = ∅ then (0 : ℝ) else lazy r / 2) ≤ m / 2
    calc
      _ ≤ ∑ r, lazy r / 2 := by
        apply Finset.sum_le_sum
        intro r _
        split_ifs
        · exact div_nonneg (hl r) (by norm_num)
        · exact le_rfl
      _ = (∑ r, lazy r) / 2 := (Finset.sum_div _ _ _).symm
      _ ≤ _ := div_le_div_of_nonneg_right ht (by norm_num)
  · intro r
    change lazy r / 2 + fresh r / (2 * 524287) ≤
      (if r = ∅ then ell / 2 + u / (2 * 524287) else 0) +
      (if r = ∅ then 0 else lazy r / 2) + _
    rw [uniform_shells]
    by_cases hr : r = ∅
    · subst r
      simp only [ite_true, add_zero]
      exact add_le_add (div_le_div_of_nonneg_right hz (by norm_num))
        (div_le_div_of_nonneg_right hfz (by norm_num))
    · simp only [if_neg hr, zero_add]
      have hm : fresh r / (2 * 524287) ≤ m / (2 * 524287) :=
        div_le_div_of_nonneg_right (hf r) (by norm_num)
      linarith

def pointRowBudget (j : ℕ) (z : ℝ) (q : Finset (Fin 19)) (live : ℝ) : Coords 5 :=
  ⟨cancellationMoment j z q / 2 + min (emittedMoment j z q) live / (2 * 524287),
    emittedMoment j z q / 2,
    fun i => ((actualShellSystem.shell i).card : ℝ) * emittedMoment j z q / (2 * 524287)⟩

theorem stepRow_dominates_point {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (j : Fin 129)
    {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    LiveDominates actualShellSystem
      (pointRowBudget j z q (Occupation.Sparse.live j (SparsePolynomial.weightN j)))
      (stepRow j z q) := by
  have h := liveDominates_branches (lazyRow j z q) (refreshRow j z q)
    (emittedMoment j z q) (cancellationMoment j z q)
    (min (emittedMoment j z q) (Occupation.Sparse.live j (SparsePolynomial.weightN j)))
    (lazyRow_nonneg hz j q) (lazyRow_total j z q).le (lazyRow_zero j z q).le
    (refreshRow_le_moment hz j q) (refreshRow_zero_le_min hz hz1 j q)
  have he : stepRow j z q = fun r => lazyRow j z q r / 2 + refreshRow j z q r / (2 * 524287) :=
    funext (stepRow_nonzero j z hq)
  rw [he]
  exact h

theorem stepRow_dominates_D {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (j : Fin 129)
    {q : Finset (Fin 19)} (hq : q ≠ ∅)
    (hc : cancellationMoment j z q ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelD) :
    LiveDominates actualShellSystem
      (Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowD (stepRow j z q) := by
  apply (stepRow_dominates_point hz hz1 j hq).mono
  · exact add_le_add (div_le_div_of_nonneg_right hc (by norm_num))
      (div_le_div_of_nonneg_right (min_le_min (emittedMoment_le_maximum j z hq) le_rfl) (by norm_num))
  · exact div_le_div_of_nonneg_right (emittedMoment_le_maximum j z hq) (by norm_num)
  · intro i
    change ((actualShellSystem.shell i).card : ℝ) * emittedMoment j z q / (2 * 524287) ≤
      Occupation.Sparse.count i * Occupation.Sparse.maximumMoment j z / (2 * 524287)
    rw [actualShellSystem_card]
    exact div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left
      (emittedMoment_le_maximum j z hq) (Occupation.Sparse.count_pos i).le) (by norm_num)

end Spin.Structured.ConcreteMaps
