import SpinCodes.Structured.DenseOccupationOuterAffine
import SpinCodes.Structured.ConcreteOuterCountingProduct

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset ConcreteOuter

/-- Actual messages with fixed active positions and a prescribed encoded weight at each. -/
def profileMessages {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (weights : Fin L → ℕ) : Finset (Fin L → LocalMessage k) :=
  univ.filter (fun x => ∀ i, if i ∈ S then
    x i ≠ 0 ∧ wtF (encode seed (x i)) = weights i else x i = 0)

theorem profileMessages_card {L k : ℕ} (seed : Seed k) (S : Finset (Fin L))
    (weights : Fin L → ℕ) :
    ((profileMessages seed S weights).card:ℝ) =
      ∏ i ∈ S, (spectrum seed (weights i):ℝ) := by
  classical
  let f : Fin L → LocalMessage k → ℝ := fun i x =>
    if (if i ∈ S then x ≠ 0 ∧ wtF (encode seed x) = weights i else x = 0) then 1 else 0
  have hf : ∀ i, ∑ x, f i x = if i ∈ S then (spectrum seed (weights i):ℝ) else 1 := by
    intro i
    by_cases hi : i ∈ S
    · simp only [f, hi, ite_true, ConcreteOuter.spectrum, Finset.card_filter, Nat.cast_sum,
        Nat.cast_ite, Nat.cast_one, Nat.cast_zero]
    · simp [f, hi]
  have hp : ∀ x : Fin L → LocalMessage k, ∏ i, f i (x i) =
      if (∀ i, if i ∈ S then x i ≠ 0 ∧ wtF (encode seed (x i)) = weights i else x i = 0)
      then 1 else 0 := by
    intro x
    by_cases hx : ∀ i, if i ∈ S then x i ≠ 0 ∧ wtF (encode seed (x i)) = weights i else x i = 0
    · simp [f, hx]
    · rw [if_neg hx]
      obtain ⟨i,hi⟩ := not_forall.mp hx
      exact Finset.prod_eq_zero (Finset.mem_univ i) (if_neg hi)
  have h := Spin.prod_sum_eq_sum_prod f
  simp only [hf, hp] at h
  rw [Finset.prod_ite] at h
  simp only [Finset.prod_const_one, mul_one] at h
  simpa only [profileMessages, Finset.card_filter, Nat.cast_sum, Nat.cast_ite,
    Nat.cast_one, Nat.cast_zero, Finset.filter_mem_eq_inter, Finset.univ_inter] using h.symm

#print axioms profileMessages_card
end Spin.Structured.DenseOccupationFixed
