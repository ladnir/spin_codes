import SpinCodes.Structured.SparseCancellation

/-! Nonnegativity of the frozen-data occupation matrix, including the
zero-live-mass endpoint rows. -/

noncomputable section
namespace Spin.Imt.Occupation

theorem matrix_nonneg {k : ℕ} (n : ℕ) (count : Fin k → ℝ) (M β : ℝ)
    (r : Fin (n + 1) → RowData k) (h0 : 0 ≤ β) (h1 : β ≤ 1)
    (hr : ∀ j, (fixed count M (r j)).Nonneg) : (matrix n count M β r).Nonneg := by
  have hh (c : Fin (n + 1) → Coords k) (p : Fin (n + 1) → ℝ)
      (hc : ∀ j, (c j).Nonneg) (hp : ∀ j, 0 ≤ p j) :
      (Coords.sum fun j => Coords.smul (p j) (c j)).Nonneg := by
    refine ⟨?_, ?_, fun i => ?_⟩ <;>
      simp only [Coords.sum, Coords.smul] <;> apply Finset.sum_nonneg <;> intro j _
    · exact mul_nonneg (hp j) (hc j).1
    · exact mul_nonneg (hp j) (hc j).2.1
    · exact mul_nonneg (hp j) ((hc j).2.2 i)
  simp only [Transfer.Nonneg, matrix]
  refine ⟨hh (fun j : Fin (n + 1) => (fixed count M (r j)).rowZ)
      (fun j => probability n j β) (fun j => (hr j).1) (fun j => probability_nonneg n j h0 h1),
    hh (fun j : Fin (n + 1) => (fixed count M (r j)).rowD)
      (fun j => probability n j β) (fun j => (hr j).2.1) (fun j => probability_nonneg n j h0 h1),
    fun i => hh (fun j : Fin (n + 1) => (fixed count M (r j)).rowS i)
      (fun j => probability n j β) (fun j => (hr j).2.2 i) (fun j => probability_nonneg n j h0 h1)⟩

namespace Sparse
open Spin.Structured Spin.Structured.SparsePolynomial

theorem foldr_min_nonneg (xs : List ℝ) (base : ℝ) (hb : 0 ≤ base)
    (hx : ∀ x ∈ xs, 0 ≤ x) : 0 ≤ xs.foldr min base := by
  induction xs with
  | nil => exact hb
  | cons x xs ih => exact le_min (hx x (by simp)) (ih fun y hy => hx y (by simp [hy]))

theorem foldr_max_nonneg (xs : List ℝ) : 0 ≤ xs.foldr max 0 := by
  induction xs with
  | nil => exact le_rfl
  | cons x xs ih => exact ih.trans (le_max_right _ _)

theorem live_bounds {j : ℕ} (hj : j ≤ 128) (d : WeightData)
    (hk : d.kernel ≤ polyChoose 128 j) : 0 ≤ live j d ∧ live j d ≤ 1 := by
  have hn : 0 < (Nat.choose 128 j : ℝ) := by exact_mod_cast Nat.choose_pos hj
  have hk' : (d.kernel : ℝ) ≤ Nat.choose 128 j := by
    exact_mod_cast (show d.kernel ≤ Nat.choose 128 j by simpa [polyChoose_eq] using hk)
  unfold live
  exact ⟨div_nonneg (sub_nonneg.mpr hk') hn.le,
    (div_le_one hn).mpr (sub_le_self _ (Nat.cast_nonneg _))⟩

theorem cancellationBoundsD_nonneg {j : ℕ} (hj : j ≤ 128) (d : WeightData)
    (hk : d.kernel ≤ polyChoose 128 j) {z : ℝ} (hz : 0 ≤ z) :
    ∀ a ∈ cancellationBoundsD j d z, 0 ≤ a := by
  intro a ha
  simp only [cancellationBoundsD, List.mem_append, List.mem_cons,
    List.not_mem_nil, or_false] at ha
  rcases ha with (rfl | rfl | rfl) | ha
  · exact foldr_max_nonneg _
  · exact (live_bounds hj d hk).1
  · positivity
  · cases h : d.lowPatterns with
    | none => simp [h] at ha
    | some ps =>
      simp only [h, List.mem_singleton] at ha
      subst a
      exact foldr_max_nonneg _

theorem cancellationBoundsS_nonneg (j : ℕ) (d : WeightData) (i : Fin 5)
    {z : ℝ} (hz : 0 ≤ z) : ∀ a ∈ cancellationBoundsS j d z i, 0 ≤ a := by
  intro a ha
  simp only [cancellationBoundsS, List.mem_append, List.mem_cons,
    List.not_mem_nil, or_false] at ha
  rcases ha with (rfl | rfl) | ha
  · exact hyperMoment_nonneg _ _ hz
  · positivity
  · cases h : d.lowPatterns with
    | none => simp [h] at ha
    | some ps =>
      simp only [h, List.mem_singleton] at ha
      subst a
      exact div_nonneg (pattern_nonneg _ _ hz) (Nat.cast_nonneg _)

theorem fixed_row_nonneg {j : ℕ} (hj : j ≤ 128) (d : WeightData)
    (hk : d.kernel ≤ polyChoose 128 j) {z : ℝ} (hz : 0 ≤ z) :
    (fixed count 524287 (row j d z)).Nonneg := by
  apply fixed_nonneg _ _ _ (fun i => (count_pos i).le) (by norm_num)
    (live_bounds hj d hk).1 (live_bounds hj d hk).2 (pow_nonneg hz _)
    (foldr_max_nonneg _) (fun i => hyperMoment_nonneg _ _ hz)
  · exact foldr_min_nonneg _ _ (foldr_max_nonneg _) (cancellationBoundsD_nonneg hj d hk hz)
  · intro i
    exact foldr_min_nonneg _ _ (hyperMoment_nonneg _ _ hz) (cancellationBoundsS_nonneg j d i hz)

end Sparse
end Spin.Imt.Occupation
