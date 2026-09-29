import SpinCodes.FiniteLaw

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin
open Finset

theorem prod_sum_eq_sum_prod {ι : Type*} [Fintype ι] [DecidableEq ι]
    {α : ι → Type*} [∀ i, Fintype (α i)] (f : ∀ i, α i → ℝ) :
    (∏ i, ∑ x, f i x) = ∑ x : ∀ i, α i, ∏ i, f i (x i) := by
  have h := Finset.prod_univ_sum (fun i => (univ : Finset (α i))) f
  rwa [Fintype.piFinset_univ] at h

theorem piPMF_bind {ι : Type*} [Fintype ι] [DecidableEq ι]
    {α β : ι → Type*} [∀ i, Fintype (α i)] [∀ i, Fintype (β i)]
    (P : ∀ i, FinPMF (α i)) (K : ∀ i, α i → FinPMF (β i)) :
    (piPMF P).bind (fun x => piPMF (fun i => K i (x i))) =
      piPMF (fun i => (P i).bind (K i)) := by
  ext y
  simp only [FinPMF.bind_p, piPMF_apply, ← Finset.prod_mul_distrib]
  exact (prod_sum_eq_sum_prod (fun i x => (P i).p x * (K i x).p (y i))).symm

theorem piPMF_map {ι : Type*} [Fintype ι] [DecidableEq ι]
    {α β : ι → Type*} [∀ i, Fintype (α i)] [∀ i, Fintype (β i)]
    (P : ∀ i, FinPMF (α i)) (f : ∀ i, α i → β i) :
    (piPMF P).map (fun x i => f i (x i)) = piPMF (fun i => (P i).map (f i)) := by
  classical
  ext y
  simp only [FinPMF.map_p, piPMF_apply]
  rw [prod_sum_eq_sum_prod]
  apply sum_congr rfl
  intro x _
  by_cases he : (fun i => f i (x i)) = y
  · subst y
    simp
  · rw [if_neg he]
    obtain ⟨i, hi⟩ : ∃ i, f i (x i) ≠ y i := by
      by_contra hc
      simp only [not_exists, not_not] at hc
      exact he (funext hc)
    symm
    exact Finset.prod_eq_zero (Finset.mem_univ i) (by simp [hi])

namespace FinPMF
variable {α β : Type*} [Fintype α] [Fintype β]

theorem map_equiv_apply (P : FinPMF α) (e : α ≃ β) (y : β) :
    (P.map e).p y = P.p (e.symm y) := by
  classical
  rw [map_p]
  have h (x : α) : e x = y ↔ x = e.symm y := e.apply_eq_iff_eq_symm_apply
  simp only [h]
  simp

theorem expect_prod_mul (P : FinPMF α) (Q : FinPMF β) (f : α → ℝ) (g : β → ℝ) :
    (P.prod Q).expect (fun xy => f xy.1 * g xy.2) = P.expect f * Q.expect g := by
  rw [expect_prod]
  simp only [expect, mul_sum, sum_mul]
  rw [sum_comm]
  apply sum_congr rfl
  intro x _
  apply sum_congr rfl
  intro y _
  ring

end FinPMF
end Spin
