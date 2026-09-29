import SpinCodes.Structured.Domination

/-! Pushforward and finite stochastic composition for the concrete encoder
and routing experiments. -/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin
open Finset
namespace FinPMF
variable {α β γ : Type*} [Fintype α] [Fintype β] [Fintype γ]

def map (P : FinPMF α) (f : α → β) : FinPMF β := by
  classical
  exact {
    p := fun y => ∑ x, if f x = y then P.p x else 0
    nonneg := fun y => sum_nonneg fun x _ => by split_ifs <;> first | exact P.nonneg x | exact le_rfl
    total := by rw [sum_comm]; simpa using P.total }

def bind (P : FinPMF α) (K : α → FinPMF β) : FinPMF β where
  p y := ∑ x, P.p x * (K x).p y
  nonneg y := sum_nonneg fun x _ => mul_nonneg (P.nonneg x) ((K x).nonneg y)
  total := by
    rw [sum_comm]
    simp only [← mul_sum, FinPMF.total, mul_one]

theorem map_p (P : FinPMF α) (f : α → β) (y : β) :
    (P.map f).p y = ∑ x, if f x = y then P.p x else 0 := by
  rfl

theorem bind_p (P : FinPMF α) (K : α → FinPMF β) (y : β) :
    (P.bind K).p y = ∑ x, P.p x * (K x).p y := rfl

theorem expect_map (P : FinPMF α) (f : α → β) (g : β → ℝ) :
    (P.map f).expect g = P.expect (fun x => g (f x)) := by
  classical
  simp only [expect, map_p, sum_mul]
  rw [sum_comm]
  apply sum_congr rfl
  intro x _
  simp [ite_mul]

theorem expect_bind (P : FinPMF α) (K : α → FinPMF β) (f : β → ℝ) :
    (P.bind K).expect f = P.expect (fun x => (K x).expect f) := by
  simp only [expect, bind_p, sum_mul, mul_sum, mul_assoc]
  rw [sum_comm]

theorem map_prob (P : FinPMF α) (f : α → β) (s : β → Prop) [DecidablePred s] :
    (P.map f).prob s = P.prob (fun x => s (f x)) := by
  rw [prob_eq_expect_indicator, expect_map, prob_eq_expect_indicator]

theorem bind_prob (P : FinPMF α) (K : α → FinPMF β) (s : β → Prop) [DecidablePred s] :
    (P.bind K).prob s = P.expect (fun x => (K x).prob s) := by
  rw [prob_eq_expect_indicator, expect_bind]
  simp only [prob_eq_expect_indicator]

theorem eq_of_expect_eq {P Q : FinPMF α} (h : ∀ f : α → ℝ, P.expect f = Q.expect f) : P = Q := by
  classical
  ext x
  simpa [expect] using h (fun y => if y = x then 1 else 0)

theorem map_id (P : FinPMF α) : P.map id = P := by
  apply eq_of_expect_eq
  intro f
  exact expect_map _ _ _

theorem map_comp (P : FinPMF α) (f : α → β) (g : β → γ) :
    (P.map f).map g = P.map (g ∘ f) := by
  apply eq_of_expect_eq
  intro h
  simp only [expect_map, Function.comp_def]

theorem bind_assoc (P : FinPMF α) (K : α → FinPMF β) (J : β → FinPMF γ) :
    (P.bind K).bind J = P.bind (fun x => (K x).bind J) := by
  apply eq_of_expect_eq
  intro f
  simp only [expect_bind]

theorem map_bind (P : FinPMF α) (K : α → FinPMF β) (f : β → γ) :
    (P.bind K).map f = P.bind (fun x => (K x).map f) := by
  apply eq_of_expect_eq
  intro g
  simp only [expect_map, expect_bind]

theorem bind_map (P : FinPMF α) (f : α → β) (K : β → FinPMF γ) :
    (P.map f).bind K = P.bind (fun x => K (f x)) := by
  apply eq_of_expect_eq
  intro g
  simp only [expect_map, expect_bind]

end FinPMF

namespace Dominates
variable {α β : Type*} [Fintype α] [Fintype β] {P Q : FinPMF α} {c : ℝ}

theorem map (h : Dominates P Q c) (f : α → β) : Dominates (P.map f) (Q.map f) c := by
  classical
  intro y
  simp only [FinPMF.map_p, mul_sum]
  apply sum_le_sum
  intro x _
  by_cases he : f x = y <;> simp [he, h x]

theorem bind (h : Dominates P Q c) (K : α → FinPMF β) :
    Dominates (P.bind K) (Q.bind K) c := by
  intro y
  simp only [FinPMF.bind_p, mul_sum]
  apply sum_le_sum
  intro x _
  simpa only [mul_assoc] using mul_le_mul_of_nonneg_right (h x) ((K x).nonneg y)

end Dominates
end Spin
