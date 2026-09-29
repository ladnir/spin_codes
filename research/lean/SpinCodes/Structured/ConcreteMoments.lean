import SpinCodes.Structured.ConcreteCounts
import SpinCodes.Structured.SparseModel

/-! Emitted moments and elementary cancellation bounds for the actual maps. -/

noncomputable section
namespace Spin.Structured.ConcreteMaps
open Finset
open scoped symmDiff

private theorem real_list_sum_range (f : ℕ → ℝ) (n : ℕ) :
    ((List.range n).map f).sum = ∑ i ∈ Finset.range n, f i := by
  induction n with
  | zero => simp
  | succ n ih => simp [List.range_succ, ih, Finset.sum_range_succ]

def emitted (z : ℝ) (q : Finset (Fin 19)) (x : Finset (Fin 128)) : ℝ :=
  z ^ (x ∆ Aset q).card

def emittedMoment (j : ℕ) (z : ℝ) (q : Finset (Fin 19)) : ℝ :=
  (∑ x ∈ Spin.layer 128 j, emitted z q x) / (Nat.choose 128 j : ℝ)

def cancellationMoment (j : ℕ) (z : ℝ) (q : Finset (Fin 19)) : ℝ :=
  (∑ x ∈ syndromeFiber j q, emitted z q x) / (Nat.choose 128 j : ℝ)

theorem emitted_nonneg {z : ℝ} (hz : 0 ≤ z) (q) (x) : 0 ≤ emitted z q x :=
  pow_nonneg hz _

theorem emitted_le_one {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (q) (x) :
    emitted z q x ≤ 1 := pow_le_one₀ hz hz1

theorem symmDiff_distance_le {n : ℕ} (x a : Finset (Fin n)) :
    SparsePolynomial.distance a.card x.card ≤ (x ∆ a).card := by
  have h := Spin.card_symmDiff_add_two_mul_card_inter x a
  have h₁ := Finset.card_le_card (Finset.inter_subset_left (s₁ := x) (s₂ := a))
  have h₂ := Finset.card_le_card (Finset.inter_subset_right (s₁ := x) (s₂ := a))
  unfold SparsePolynomial.distance
  omega

theorem emitted_le_distance {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (q : Finset (Fin 19)) {x : Finset (Fin 128)} {j : ℕ} (hx : x.card = j) :
    emitted z q x ≤ z ^ SparsePolynomial.distance (Aset q).card j := by
  apply pow_le_pow_of_le_one hz hz1
  rw [← hx]
  exact symmDiff_distance_le x (Aset q)

theorem emittedMoment_nonneg {z : ℝ} (hz : 0 ≤ z) (j q) :
    0 ≤ emittedMoment j z q :=
  div_nonneg (Finset.sum_nonneg fun x _ => emitted_nonneg hz q x) (Nat.cast_nonneg _)

theorem cancellationMoment_nonneg {z : ℝ} (hz : 0 ≤ z) (j q) :
    0 ≤ cancellationMoment j z q :=
  div_nonneg (Finset.sum_nonneg fun x _ => emitted_nonneg hz q x) (Nat.cast_nonneg _)

theorem cancellationMoment_le_moment {z : ℝ} (hz : 0 ≤ z) (j q) :
    cancellationMoment j z q ≤ emittedMoment j z q := by
  apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
  exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
    (fun x _ _ => emitted_nonneg hz q x)

theorem cancellationMoment_le_fiber {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1) (j q) :
    cancellationMoment j z q ≤
      ((syndromeFiber j q).card : ℝ) / Nat.choose 128 j *
        z ^ SparsePolynomial.distance (Aset q).card j := by
  unfold cancellationMoment
  rw [div_mul_eq_mul_div]
  apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
  calc
    _ ≤ ∑ _x ∈ syndromeFiber j q, z ^ SparsePolynomial.distance (Aset q).card j := by
      apply Finset.sum_le_sum
      intro x hx
      apply emitted_le_distance hz hz1
      exact (Finset.mem_powersetCard.mp (Finset.mem_filter.mp hx).1).2
    _ = _ := by rw [Finset.sum_const, nsmul_eq_mul]

theorem cancellationMoment_le_cap {z : ℝ} (hz : 0 ≤ z) (hz1 : z ≤ 1)
    (j : Fin 129) {q : Finset (Fin 19)} (hq : q ≠ ∅) :
    cancellationMoment j z q ≤ (SparsePolynomial.weightN j).cap / (Nat.choose 128 j : ℝ) *
      z ^ SparsePolynomial.distance (Aset q).card j := by
  apply (cancellationMoment_le_fiber hz hz1 j q).trans
  apply mul_le_mul_of_nonneg_right _ (pow_nonneg hz _)
  apply div_le_div_of_nonneg_right _ (Nat.cast_nonneg _)
  exact_mod_cast syndromeFiber_cap_frozen j hq

theorem emittedMoment_eq_hyper (j : Fin 129) (z : ℝ) (q : Finset (Fin 19)) :
    emittedMoment j z q = Spin.Imt.Occupation.hyperMoment (Aset q).card j z := by
  classical
  unfold emittedMoment Spin.Imt.Occupation.hyperMoment
  congr 1
  rw [real_list_sum_range]
  have hm : ∀ x ∈ Spin.layer 128 j, (x ∩ Aset q).card ∈ range 129 := by
    intro x hx
    have hc := (Finset.mem_powersetCard.mp hx).2
    have hi := Finset.card_le_card (Finset.inter_subset_left (s₁ := x) (s₂ := Aset q))
    simp only [Finset.mem_range]
    omega
  rw [← Finset.sum_fiberwise_of_maps_to hm]
  apply Finset.sum_congr rfl
  intro h _
  by_cases hh : h ≤ j
  · rw [if_pos hh]
    have he : ∀ x ∈ (Spin.layer 128 j).filter (fun x => (x ∩ Aset q).card = h),
        emitted z q x = z ^ ((Aset q).card + j - 2 * h) := by
      intro x hx
      obtain ⟨hx, hi⟩ := Finset.mem_filter.mp hx
      have hc := (Finset.mem_powersetCard.mp hx).2
      have hs := Spin.card_symmDiff_add_two_mul_card_inter x (Aset q)
      unfold emitted
      congr 1
      omega
    rw [Finset.sum_congr rfl he, Finset.sum_const, nsmul_eq_mul,
      Spin.card_meet_fiber (Aset q) hh]
    push_cast
    ring
  · rw [if_neg hh]
    apply Finset.sum_eq_zero
    intro x hx
    obtain ⟨hx, hi⟩ := Finset.mem_filter.mp hx
    have hc := (Finset.mem_powersetCard.mp hx).2
    have hb := Finset.card_le_card (Finset.inter_subset_left (s₁ := x) (s₂ := Aset q))
    omega

end Spin.Structured.ConcreteMaps
