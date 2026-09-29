import SpinCodes.Structured.ConcreteMoments

/-! General cancellation envelopes for the actual feedback fibers. The special
weight-one and weight-two cancellation tables are separate obligations. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset

theorem exists_weightShell {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    ∃ i : Fin 5, q ∈ weightShell i := by
  have h := Spin.mem_nonzeroStates hq
  rw [← weightShell_covers actual_shell_counts] at h
  obtain ⟨i, _, hi⟩ := Finset.mem_biUnion.mp h
  exact ⟨i, hi⟩

theorem weightShell_card_actual (i : Fin 5) : (weightShell i).card = shellCount i :=
  (weightShell_card i).trans (actual_shell_counts i)

private theorem le_foldr_max_of_mem {a : ℝ} {xs : List ℝ} (ha : a ∈ xs) :
    a ≤ xs.foldr max 0 := by
  induction xs with
  | nil => simp at ha
  | cons b xs ih =>
    rcases List.mem_cons.mp ha with rfl | ha
    · exact le_max_left _ _
    · exact (ih ha).trans (le_max_right _ _)

theorem shell_moment_le_maximum (j : ℕ) (z : ℝ) (i : Fin 5) :
    Spin.Imt.Occupation.hyperMoment (shellWeight i) j z ≤
      Spin.Imt.Occupation.Sparse.maximumMoment j z := by
  apply le_foldr_max_of_mem
  apply List.mem_map.mpr
  refine ⟨shellWeight i, ?_, rfl⟩
  fin_cases i <;> simp [shellWeight, SparsePolynomial.levels]

theorem minDistance_le_shell (j : ℕ) (i : Fin 5) :
    SparsePolynomial.minDistance j ≤ SparsePolynomial.distance (shellWeight i) j := by
  fin_cases i <;> simp [SparsePolynomial.minDistance, shellWeight,
    SparsePolynomial.levels, min_le_iff]

theorem emittedMoment_shell (j : Fin 129) (z : ℝ) {i : Fin 5} {q : Finset (Fin 19)}
    (hq : q ∈ weightShell i) :
    emittedMoment j z q = Spin.Imt.Occupation.hyperMoment (shellWeight i) j z := by
  rw [emittedMoment_eq_hyper, weightShell_weight hq]

theorem emittedMoment_le_maximum (j : Fin 129) (z : ℝ) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    emittedMoment j z q ≤ Spin.Imt.Occupation.Sparse.maximumMoment j z := by
  obtain ⟨i, hi⟩ := exists_weightShell hq
  rw [emittedMoment_shell j z hi]
  exact shell_moment_le_maximum j z i

theorem cancellationMoment_le_maximum {z : ℝ} (hz : 0 ≤ z) (j : Fin 129)
    {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    cancellationMoment j z q ≤ Spin.Imt.Occupation.Sparse.maximumMoment j z :=
  (cancellationMoment_le_moment hz j q).trans (emittedMoment_le_maximum j z hq)

theorem cancellationMoment_le_live {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    cancellationMoment j z q ≤ Spin.Imt.Occupation.Sparse.live j (SparsePolynomial.weightN j) := by
  have hc := nonzero_fiber_le_complement j hq
  rw [kernel_card_frozen j] at hc
  have hc' : ((syndromeFiber j q).card : ℝ) + (SparsePolynomial.weightN j).kernel ≤
      (Nat.choose 128 j : ℝ) := by exact_mod_cast hc
  unfold cancellationMoment Spin.Imt.Occupation.Sparse.live
  apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
  calc
    _ ≤ ∑ _x ∈ syndromeFiber j q, (1 : ℝ) :=
      Finset.sum_le_sum fun x _ => emitted_le_one hz hz1 q x
    _ = ((syndromeFiber j q).card : ℝ) := by simp only [Finset.sum_const, nsmul_eq_mul, mul_one]
    _ ≤ _ := by linarith

theorem cancellationMoment_le_minDistance {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    cancellationMoment j z q ≤ (SparsePolynomial.weightN j).cap / (Nat.choose 128 j : ℝ) *
      z ^ SparsePolynomial.minDistance j := by
  obtain ⟨i, hi⟩ := exists_weightShell hq
  apply (cancellationMoment_le_cap hz hz1 j hq).trans
  apply mul_le_mul_of_nonneg_left _ (div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _))
  apply pow_le_pow_of_le_one hz hz1
  rw [weightShell_weight hi]
  exact minDistance_le_shell j i

def shellCancellationMoment (j : ℕ) (z : ℝ) (i : Fin 5) : ℝ :=
  (∑ q ∈ weightShell i, cancellationMoment j z q) / (shellCount i : ℝ)

theorem shellCancellationMoment_le_moment {z : ℝ} (hz : 0 ≤ z) (j : Fin 129) (i : Fin 5) :
    shellCancellationMoment j z i ≤ Spin.Imt.Occupation.hyperMoment (shellWeight i) j z := by
  have hc : (0 : ℝ) < shellCount i := by exact_mod_cast shellCount_pos i
  unfold shellCancellationMoment
  apply (div_le_iff₀ hc).mpr
  calc
    _ ≤ ∑ _q ∈ weightShell i, Spin.Imt.Occupation.hyperMoment (shellWeight i) j z := by
      apply Finset.sum_le_sum
      intro q hq
      exact (cancellationMoment_le_moment hz j q).trans_eq (emittedMoment_shell j z hq)
    _ = _ := by rw [Finset.sum_const, nsmul_eq_mul, weightShell_card_actual]; ring

theorem shell_fiber_count_le (j : Fin 129) (i : Fin 5) :
    ∑ q ∈ weightShell i, (syndromeFiber j q).card ≤
      min (Nat.choose 128 j - (SparsePolynomial.weightN j).kernel)
        (shellCount i * (SparsePolynomial.weightN j).cap) := by
  apply le_min
  · rw [← kernel_card_frozen j, ← sum_nonzero_fiber_card]
    apply Finset.sum_le_sum_of_subset_of_nonneg
    · intro q hq
      exact Finset.mem_erase.mpr ⟨weightShell_nonzero hq, Finset.mem_univ _⟩
    · intro _ _ _; exact Nat.zero_le _
  · calc
      _ ≤ ∑ _q ∈ weightShell i, (SparsePolynomial.weightN j).cap :=
        Finset.sum_le_sum fun q hq => syndromeFiber_cap_frozen j (weightShell_nonzero hq)
      _ = _ := by rw [Finset.sum_const, smul_eq_mul, weightShell_card_actual]

theorem shellCancellationMoment_le_cap {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (i : Fin 5) :
    shellCancellationMoment j z i ≤
      (min (Nat.choose 128 j - (SparsePolynomial.weightN j).kernel)
        (shellCount i * (SparsePolynomial.weightN j).cap) : ℕ) /
        ((shellCount i : ℝ) * Nat.choose 128 j) * z ^ SparsePolynomial.distance (shellWeight i) j := by
  let b := z ^ SparsePolynomial.distance (shellWeight i) j
  have hn : (0 : ℝ) ≤ Nat.choose 128 j := Nat.cast_nonneg _
  have hs : (0 : ℝ) ≤ shellCount i := Nat.cast_nonneg _
  have hb : 0 ≤ b := pow_nonneg hz _
  have he : ∑ q ∈ weightShell i, cancellationMoment j z q ≤
      (∑ q ∈ weightShell i, ((syndromeFiber j q).card : ℝ)) / Nat.choose 128 j * b := by
    calc
      _ ≤ ∑ q ∈ weightShell i, ((syndromeFiber j q).card : ℝ) / Nat.choose 128 j * b := by
        apply Finset.sum_le_sum
        intro q hq
        simpa only [weightShell_weight hq] using cancellationMoment_le_fiber hz hz1 j q
      _ = _ := by rw [← Finset.sum_mul, ← Finset.sum_div]
  have hc : (∑ q ∈ weightShell i, ((syndromeFiber j q).card : ℝ)) ≤
      (min (Nat.choose 128 j - (SparsePolynomial.weightN j).kernel)
        (shellCount i * (SparsePolynomial.weightN j).cap) : ℕ) := by
    exact_mod_cast shell_fiber_count_le j i
  unfold shellCancellationMoment
  apply (div_le_div_of_nonneg_right he hs).trans
  calc
    _ ≤ ((min (Nat.choose 128 j - (SparsePolynomial.weightN j).kernel)
        (shellCount i * (SparsePolynomial.weightN j).cap) : ℕ) /
        (Nat.choose 128 j : ℝ) * b) / shellCount i :=
      div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_right
        (div_le_div_of_nonneg_right hc hn) hb) hs
    _ = _ := by dsimp [b]; ring

end Spin.Structured.ConcreteMaps
