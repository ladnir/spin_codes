import SpinCodes.Structured.ConcreteRows

/-! The exact zero row and the all-kernel branch of the actual transition. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open Spin.Imt
open scoped symmDiff

theorem kernelLayer_card_le (j : ℕ) : (kernelLayer j).card ≤ Nat.choose 128 j := by
  rw [← Spin.card_layer]
  exact Finset.card_le_card (Finset.filter_subset _ _)

theorem stepRow_zero_zero (j : Fin 129) (z : ℝ) :
    stepRow j z ∅ ∅ = (1 - Occupation.Sparse.live j (SparsePolynomial.weightN j)) * z ^ j.val := by
  have hn : (Nat.choose 128 j : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.choose_pos (Nat.le_of_lt_succ j.isLt)).ne'
  rw [stepRow_from_zero, syndromeFiber_zero, kernel_card_frozen]
  unfold Occupation.Sparse.live
  field_simp
  ring

theorem stepRow_zero_nonzero_total (j : Fin 129) (z : ℝ) :
    ∑ r ∈ nonzeroSyndromes, stepRow j z ∅ r =
      Occupation.Sparse.live j (SparsePolynomial.weightN j) * z ^ j.val := by
  simp only [stepRow_from_zero]
  rw [← Finset.sum_mul, ← Finset.sum_div, ← Nat.cast_sum, sum_nonzero_fiber_card,
    Nat.cast_sub (kernelLayer_card_le j), kernel_card_frozen j]
  rfl

theorem stepRow_dominates_Z {z : ℝ} (hz : 0 ≤ z) (j : Fin 129) :
    LiveDominates actualShellSystem
      (Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowZ (stepRow j z ∅) := by
  classical
  rw [liveDominates_iff_witness]
  refine ⟨fun r => if r = ∅ then 0 else stepRow j z ∅ r, ?_, by simp, ?_, ?_⟩
  · intro r
    change 0 ≤ if r = ∅ then 0 else stepRow j z ∅ r
    split_ifs
    · exact le_rfl
    · exact stepRow_nonneg hz j ∅ r
  · change (∑ r, if r = ∅ then (0 : ℝ) else stepRow j z ∅ r) ≤
      Occupation.Sparse.live j (SparsePolynomial.weightN j) * z ^ j.val
    have he : (univ : Finset (Finset (Fin 19))).filter (fun r => r ≠ ∅) = nonzeroSyndromes := by
      ext r
      simp [nonzeroSyndromes]
    have hs := stepRow_zero_nonzero_total j z
    rw [← he, Finset.sum_filter] at hs
    simpa only [ite_not] using hs.le
  · intro r
    change stepRow j z ∅ r ≤
      (if r = ∅ then (1 - Occupation.Sparse.live j (SparsePolynomial.weightN j)) * z ^ j.val else 0) +
      (if r = ∅ then 0 else stepRow j z ∅ r) +
      ∑ i, if r ∈ actualShellSystem.shell i then (0 : ℝ) / (actualShellSystem.shell i).card else 0
    simp only [zero_div, ite_self, Finset.sum_const_zero, add_zero]
    by_cases hr : r = ∅
    · subst r
      simp only [ite_true, add_zero, stepRow_zero_zero]
      exact le_rfl
    · simp only [if_neg hr, zero_add, le_refl]

theorem live_zero_syndrome (j : Fin 129)
    (hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0)
    {x : Finset (Fin 128)} (hx : x ∈ Spin.layer 128 j) : Cset x = ∅ := by
  have hn : (Nat.choose 128 j : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.choose_pos (Nat.le_of_lt_succ j.isLt)).ne'
  have he : (Nat.choose 128 j : ℝ) = (SparsePolynomial.weightN j).kernel := by
    exact sub_eq_zero.mp ((div_eq_zero_iff.mp hl).resolve_right hn)
  have hc : (kernelLayer j).card = (Spin.layer 128 j).card := by
    rw [kernel_card_frozen, Spin.card_layer]
    exact_mod_cast he.symm
  have hf : kernelLayer j = Spin.layer 128 j :=
    Finset.eq_of_subset_of_card_le (Finset.filter_subset _ _) hc.ge
  rw [← hf] at hx
  exact (Finset.mem_filter.mp hx).2

theorem stepRow_all_kernel (j : Fin 129) (z : ℝ)
    (hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0)
    (q r : Finset (Fin 19)) :
    stepRow j z q r = emittedMoment j z q * Spin.Imt.actLaw q r := by
  unfold stepRow emittedMoment
  have he : ∀ x ∈ Spin.layer 128 j,
      emitted z q x * Spin.Imt.actLaw q (r ∆ Cset x) = emitted z q x * Spin.Imt.actLaw q r := by
    intro x hx
    rw [live_zero_syndrome j hl hx, show r ∆ (∅ : Finset (Fin 19)) = r from symmDiff_bot r]
  rw [Finset.sum_congr rfl he, ← Finset.sum_mul]
  ring

end Spin.Structured.ConcreteMaps
