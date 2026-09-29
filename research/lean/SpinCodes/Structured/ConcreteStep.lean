import SpinCodes.Structured.ConcreteCancellation
import SpinCodes.Structured.LiveInduction

/-! The weighted transition for uniform fixed-weight input and the actual
transvection law, decomposed into its lazy and refresh branches. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open scoped symmDiff

def stepRow (j : ℕ) (z : ℝ) (q r : Finset (Fin 19)) : ℝ :=
  (∑ x ∈ Spin.layer 128 j, emitted z q x * Spin.Imt.actLaw q (r ∆ Cset x)) /
    (Nat.choose 128 j : ℝ)

def lazyRow (j : ℕ) (z : ℝ) (q r : Finset (Fin 19)) : ℝ :=
  (∑ x ∈ Spin.layer 128 j, if q ∆ Cset x = r then emitted z q x else 0) /
    (Nat.choose 128 j : ℝ)

def refreshRow (j : ℕ) (z : ℝ) (q r : Finset (Fin 19)) : ℝ :=
  (∑ x ∈ Spin.layer 128 j, if Cset x = r then 0 else emitted z q x) /
    (Nat.choose 128 j : ℝ)

theorem nonzeroStates_card_actual : (Spin.nonzeroStates 19).card = 524287 := by
  simp [Spin.nonzeroStates_eq_erase]

theorem stepRow_nonzero (j : ℕ) (z : ℝ) {q : Finset (Fin 19)} (hq : q ≠ ∅)
    (r : Finset (Fin 19)) :
    stepRow j z q r = lazyRow j z q r / 2 + refreshRow j z q r / (2 * 524287) := by
  classical
  have he (x : Finset (Fin 128)) :
      emitted z q x * Spin.Imt.actLaw q (r ∆ Cset x) =
        (if q ∆ Cset x = r then emitted z q x else 0) / 2 +
        (if Cset x = r then 0 else emitted z q x) / (2 * 524287) := by
    have hs : r ∆ Cset x = q ↔ q ∆ Cset x = r := by
      constructor <;> intro h
      · have h' := congrArg (fun w => w ∆ Cset x) h
        simpa only [Spin.symmDiff_cancel_tail] using h'.symm
      · have h' := congrArg (fun w => w ∆ Cset x) h
        simpa only [Spin.symmDiff_cancel_tail] using h'.symm
    have hzero : r ∆ Cset x = ∅ ↔ Cset x = r := by
      rw [Finset.symmDiff_eq_empty]
      exact eq_comm
    rw [Spin.Imt.actLaw, Spin.prob_act hq, nonzeroStates_card_actual]
    simp only [hs, hzero, Nat.cast_ofNat]
    split_ifs <;> ring
  unfold stepRow lazyRow refreshRow
  rw [Finset.sum_congr rfl (fun x _ => he x), Finset.sum_add_distrib]
  simp only [← Finset.sum_div]
  ring

theorem lazyRow_nonneg {z : ℝ} (hz : 0 ≤ z) (j q r) : 0 ≤ lazyRow j z q r := by
  apply div_nonneg _ (Nat.cast_nonneg _)
  apply Finset.sum_nonneg
  intro x _
  split_ifs <;> first | exact emitted_nonneg hz q x | exact le_rfl

theorem refreshRow_nonneg {z : ℝ} (hz : 0 ≤ z) (j q r) : 0 ≤ refreshRow j z q r := by
  apply div_nonneg _ (Nat.cast_nonneg _)
  apply Finset.sum_nonneg
  intro x _
  split_ifs <;> first | exact emitted_nonneg hz q x | exact le_rfl

theorem stepRow_nonneg {z : ℝ} (hz : 0 ≤ z) (j q r) : 0 ≤ stepRow j z q r := by
  apply div_nonneg _ (Nat.cast_nonneg _)
  apply Finset.sum_nonneg
  intro x _
  exact mul_nonneg (emitted_nonneg hz q x)
    (div_nonneg (Nat.cast_nonneg _) (Nat.cast_nonneg _))

theorem lazyRow_zero (j : ℕ) (z : ℝ) (q : Finset (Fin 19)) :
    lazyRow j z q ∅ = cancellationMoment j z q := by
  unfold lazyRow cancellationMoment syndromeFiber
  rw [Finset.sum_filter]
  simp only [Finset.symmDiff_eq_empty, eq_comm (a := q)]

theorem lazyRow_total (j : ℕ) (z : ℝ) (q : Finset (Fin 19)) :
    ∑ r, lazyRow j z q r = emittedMoment j z q := by
  classical
  unfold lazyRow emittedMoment
  rw [← Finset.sum_div, Finset.sum_comm]
  congr 1
  apply Finset.sum_congr rfl
  intro x _
  simp only [Finset.sum_ite_eq, Finset.mem_univ, if_true]

theorem refreshRow_le_moment {z : ℝ} (hz : 0 ≤ z) (j q r) :
    refreshRow j z q r ≤ emittedMoment j z q := by
  apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
  apply Finset.sum_le_sum
  intro x _
  split_ifs <;> first | exact emitted_nonneg hz q x | exact le_rfl

theorem stepRow_zero_column (j : ℕ) (z : ℝ) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    stepRow j z q ∅ = cancellationMoment j z q / 2 + refreshRow j z q ∅ / (2 * 524287) := by
  rw [stepRow_nonzero j z hq, lazyRow_zero]

theorem refreshRow_zero_le_live {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (q : Finset (Fin 19)) :
    refreshRow j z q ∅ ≤ Spin.Imt.Occupation.Sparse.live j (SparsePolynomial.weightN j) := by
  classical
  have hc := Finset.card_filter_add_card_filter_not (s := Spin.layer 128 j)
    (fun x => Cset x = ∅)
  have hkernel : ((Spin.layer 128 j).filter (fun x => Cset x = ∅)).card =
      (SparsePolynomial.weightN j).kernel := kernel_card_frozen j
  rw [hkernel, Spin.card_layer] at hc
  have hc' : (((Spin.layer 128 j).filter (fun x => Cset x ≠ ∅)).card : ℝ) =
      (Nat.choose 128 j : ℝ) - (SparsePolynomial.weightN j).kernel := by
    have ht : ((SparsePolynomial.weightN j).kernel : ℝ) +
        (((Spin.layer 128 j).filter (fun x => Cset x ≠ ∅)).card : ℝ) =
        (Nat.choose 128 j : ℝ) := by exact_mod_cast hc
    linarith
  have h := (Spin.sum_refresh_le_min (Spin.layer 128 j) (emitted z q) Cset
    (fun x _ => emitted_nonneg hz q x) (fun x _ => emitted_le_one hz hz1 q x)).trans
      (min_le_right _ _)
  have he : (∑ x ∈ Spin.layer 128 j, emitted z q x * (if Cset x = ∅ then (0 : ℝ) else 1)) =
      ∑ x ∈ Spin.layer 128 j, if Cset x = ∅ then 0 else emitted z q x := by
    apply Finset.sum_congr rfl
    intro x _
    split_ifs <;> simp
  rw [he, hc'] at h
  exact div_le_div_of_nonneg_right h (Nat.cast_nonneg _)

theorem refreshRow_zero_le_min {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) (q : Finset (Fin 19)) :
    refreshRow j z q ∅ ≤ min (emittedMoment j z q)
      (Spin.Imt.Occupation.Sparse.live j (SparsePolynomial.weightN j)) :=
  le_min (refreshRow_le_moment hz j q ∅) (refreshRow_zero_le_live hz hz1 j q)

private theorem transPairs_card_pos_actual : 0 < (Spin.transPairs 19).card := by
  have hq : ({0} : Finset (Fin 19)) ≠ ∅ := Finset.singleton_ne_empty _
  rw [Spin.card_transPairs hq]
  exact Nat.mul_pos (Finset.card_pos.mpr ⟨{0}, Spin.mem_nonzeroStates hq⟩)
    (Spin.card_perp_pos hq)

theorem actLaw_empty_actual (r : Finset (Fin 19)) :
    Spin.Imt.actLaw ∅ r = if r = ∅ then 1 else 0 := by
  classical
  have hn : ((Spin.transPairs 19).card : ℝ) ≠ 0 := by
    exact_mod_cast (Nat.ne_of_gt transPairs_card_pos_actual)
  unfold Spin.Imt.actLaw
  by_cases hr : r = ∅
  · subst r
    simp only [Spin.act_empty, filter_true, ite_true, div_self hn]
  · have hh : (Spin.transPairs 19).filter (fun p => Spin.act p.1 p.2 ∅ = r) = ∅ := by
      apply Finset.filter_eq_empty_iff.mpr
      intro p _
      simpa only [Spin.act_empty] using Ne.symm hr
    rw [hh]
    simp only [Finset.card_empty, Nat.cast_zero, zero_div, if_neg hr]

theorem stepRow_from_zero (j : ℕ) (z : ℝ) (r : Finset (Fin 19)) :
    stepRow j z ∅ r = ((syndromeFiber j r).card : ℝ) / Nat.choose 128 j * z ^ j := by
  classical
  have he (x : Finset (Fin 128)) (hx : x ∈ Spin.layer 128 j) :
      emitted z ∅ x * Spin.Imt.actLaw ∅ (r ∆ Cset x) =
        if Cset x = r then z ^ j else 0 := by
    have hc := (Finset.mem_powersetCard.mp hx).2
    rw [emitted, Aset_empty, show x ∆ (∅ : Finset (Fin 128)) = x from symmDiff_bot x,
      hc, actLaw_empty_actual]
    simp only [Finset.symmDiff_eq_empty, eq_comm (a := r)]
    split_ifs <;> simp
  unfold stepRow
  rw [Finset.sum_congr rfl he, ← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul]
  change ((syndromeFiber j r).card : ℝ) * z ^ j / Nat.choose 128 j = _
  ring

end Spin.Structured.ConcreteMaps
