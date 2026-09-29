import SpinCodes.Structured.ConcreteZeroRow
import SpinCodes.Structured.ConcreteCancellationRows

/-! Averaging the actual transition over an entering weight shell. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open Spin.Imt

def shellStepRow (j : ℕ) (z : ℝ) (i : Fin 5) (r : Finset (Fin 19)) : ℝ :=
  (∑ q ∈ weightShell i, stepRow j z q r) / (shellCount i : ℝ)

def shellLazyRow (j : ℕ) (z : ℝ) (i : Fin 5) (r : Finset (Fin 19)) : ℝ :=
  (∑ q ∈ weightShell i, lazyRow j z q r) / (shellCount i : ℝ)

def shellRefreshRow (j : ℕ) (z : ℝ) (i : Fin 5) (r : Finset (Fin 19)) : ℝ :=
  (∑ q ∈ weightShell i, refreshRow j z q r) / (shellCount i : ℝ)

theorem shell_average_le (i : Fin 5) (f : Finset (Fin 19) → ℝ) (a : ℝ)
    (h : ∀ q ∈ weightShell i, f q ≤ a) :
    (∑ q ∈ weightShell i, f q) / (shellCount i : ℝ) ≤ a := by
  have hc : (0 : ℝ) < shellCount i := by exact_mod_cast shellCount_pos i
  apply (div_le_iff₀ hc).mpr
  calc
    _ ≤ ∑ _q ∈ weightShell i, a := Finset.sum_le_sum h
    _ = _ := by rw [Finset.sum_const, nsmul_eq_mul, weightShell_card_actual]; ring

theorem shell_average_eq (i : Fin 5) (f : Finset (Fin 19) → ℝ) (a : ℝ)
    (h : ∀ q ∈ weightShell i, f q = a) :
    (∑ q ∈ weightShell i, f q) / (shellCount i : ℝ) = a := by
  have hc : (shellCount i : ℝ) ≠ 0 := by exact_mod_cast (shellCount_pos i).ne'
  rw [Finset.sum_congr rfl h, Finset.sum_const, nsmul_eq_mul, weightShell_card_actual]
  exact mul_div_cancel_left₀ _ hc

theorem shellStepRow_branches (j : ℕ) (z : ℝ) (i : Fin 5) (r : Finset (Fin 19)) :
    shellStepRow j z i r = shellLazyRow j z i r / 2 + shellRefreshRow j z i r / (2 * 524287) := by
  unfold shellStepRow shellLazyRow shellRefreshRow
  rw [Finset.sum_congr rfl (fun q hq => stepRow_nonzero j z (weightShell_nonzero hq) r),
    Finset.sum_add_distrib]
  simp only [← Finset.sum_div]
  ring

theorem shellLazyRow_nonneg {z : ℝ} (hz : 0 ≤ z) (j i r) : 0 ≤ shellLazyRow j z i r :=
  div_nonneg (Finset.sum_nonneg fun q _ => lazyRow_nonneg hz j q r) (Nat.cast_nonneg _)

theorem shellLazyRow_total (j : Fin 129) (z : ℝ) (i : Fin 5) :
    ∑ r, shellLazyRow j z i r = Occupation.hyperMoment (shellWeight i) j z := by
  unfold shellLazyRow
  rw [← Finset.sum_div, Finset.sum_comm]
  simp only [lazyRow_total]
  exact shell_average_eq i _ _ (fun q hq => emittedMoment_shell j z hq)

theorem shellLazyRow_zero (j : ℕ) (z : ℝ) (i : Fin 5) :
    shellLazyRow j z i ∅ = shellCancellationMoment j z i := by
  simp only [shellLazyRow, shellCancellationMoment, lazyRow_zero]

theorem shellRefreshRow_le_moment {z : ℝ} (hz : 0 ≤ z) (j : Fin 129) (i : Fin 5) (r) :
    shellRefreshRow j z i r ≤ Occupation.hyperMoment (shellWeight i) j z := by
  apply shell_average_le
  intro q hq
  exact (refreshRow_le_moment hz j q r).trans_eq (emittedMoment_shell j z hq)

theorem shellRefreshRow_zero_le_min {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (i : Fin 5) :
    shellRefreshRow j z i ∅ ≤ min (Occupation.hyperMoment (shellWeight i) j z)
      (Occupation.Sparse.live j (SparsePolynomial.weightN j)) := by
  apply shell_average_le
  intro q hq
  have h := refreshRow_zero_le_min hz hz1 j q
  rwa [emittedMoment_shell j z hq] at h

theorem shellStepRow_dominates_live {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (i : Fin 5)
    (hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) ≠ 0)
    (hc : shellCancellationMoment j z i ≤ (Occupation.Sparse.row j (SparsePolynomial.weightN j) z).cancelS i) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowS i) (shellStepRow j z i) := by
  have h := liveDominates_branches (shellLazyRow j z i) (shellRefreshRow j z i)
    (Occupation.hyperMoment (shellWeight i) j z) (shellCancellationMoment j z i)
    (min (Occupation.hyperMoment (shellWeight i) j z) (Occupation.Sparse.live j (SparsePolynomial.weightN j)))
    (shellLazyRow_nonneg hz j i) (shellLazyRow_total j z i).le (shellLazyRow_zero j z i).le
    (shellRefreshRow_le_moment hz j i) (shellRefreshRow_zero_le_min hz hz1 j i)
  have he : shellStepRow j z i = fun r => shellLazyRow j z i r / 2 + shellRefreshRow j z i r / (2 * 524287) :=
    funext (shellStepRow_branches j z i)
  rw [← he] at h
  apply h.mono
  · exact add_le_add (div_le_div_of_nonneg_right hc (by norm_num : (0 : ℝ) ≤ 2)) le_rfl
  · change Occupation.hyperMoment (shellWeight i) j z / 2 ≤
      if Occupation.Sparse.live j (SparsePolynomial.weightN j) = 0 then 0
        else Occupation.hyperMoment (shellWeight i) j z / 2
    rw [if_neg hl]
  · intro h
    simp only [Occupation.fixed, Occupation.Sparse.row, hl, false_and, ite_false, add_zero,
      actualShellSystem_card]
    exact le_rfl

theorem shellStepRow_dominates_live_generic {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (hj1 : j.val ≠ 1) (hj2 : j.val ≠ 2) (i : Fin 5)
    (hl : Occupation.Sparse.live j (SparsePolynomial.weightN j) ≠ 0) :
    LiveDominates actualShellSystem
      ((Occupation.fixed Occupation.Sparse.count 524287
        (Occupation.Sparse.row j (SparsePolynomial.weightN j) z)).rowS i) (shellStepRow j z i) :=
  shellStepRow_dominates_live hz hz1 j i hl (shellCancellationMoment_le_rowS hz hz1 j hj1 hj2 i)

end Spin.Structured.ConcreteMaps
