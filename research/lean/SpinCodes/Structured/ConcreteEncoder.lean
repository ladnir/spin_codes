import SpinCodes.Structured.ConcreteOccupation
import SpinCodes.Structured.Domination

/-! The finite encoder experiment. Each round consumes one 128-bit input,
emits its XOR with the expanded state, and updates the state using an
independently sampled valid transvection. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open scoped symmDiff
open ConcreteMaps

abbrev Input := Finset (Fin 128)
abbrev State := Finset (Fin 19)
abbrev Transvection := {p : State × State // p ∈ Spin.transPairs 19}

instance : Nonempty Transvection := ⟨⟨({0}, ∅), by simp [Spin.transPairs]⟩⟩

def transvectionLaw : Spin.FinPMF Transvection := Spin.FinPMF.uniform _

def nextState (q : State) (x : Input) (t : Transvection) : State :=
  Spin.act t.val.1 t.val.2 q ∆ Cset x

/-- Total emitted Hamming weight, including all rounds and no terminal tail. -/
def outputWeight : {R : Nat} → (Fin R → Input) → (Fin R → Transvection) → State → Nat
  | 0, _, _, _ => 0
  | R + 1, xs, ts, q => (xs 0 ∆ Aset q).card +
      outputWeight (Fin.tail xs) (Fin.tail ts) (nextState q (xs 0) (ts 0))

/-- Independent transvections averaged for arbitrary deterministic input blocks. -/
def inputMoment {R : Nat} (z : ℝ) (xs : Fin R → Input) (q : State := ∅) : ℝ :=
  (Spin.piPMF (fun _ : Fin R => transvectionLaw)).expect
    (fun ts => z ^ outputWeight xs ts q)

theorem inputMoment_nonneg {R : Nat} {z : ℝ} (hz : 0 ≤ z)
    (xs : Fin R → Input) (q : State) : 0 ≤ inputMoment z xs q :=
  Spin.FinPMF.expect_nonneg _ (fun _ => pow_nonneg hz _)

theorem sum_tuple_succ {α : Type*} [Fintype α] {R : Nat}
    (f : (Fin (R + 1) → α) → ℝ) :
    (∑ xs, f xs) = ∑ x, ∑ xs : Fin R → α, f (Fin.cons x xs) := by
  rw [← Equiv.sum_comp (Fin.consEquiv (fun _ : Fin (R + 1) => α))]
  rw [Fintype.sum_prod_type]
  rfl

theorem expect_iid_succ {α : Type*} [Fintype α] (P : Spin.FinPMF α) {R : Nat}
    (f : (Fin (R + 1) → α) → ℝ) :
    (Spin.piPMF (fun _ : Fin (R + 1) => P)).expect f =
      P.expect (fun x => (Spin.piPMF (fun _ : Fin R => P)).expect
        (fun xs => f (Fin.cons x xs))) := by
  classical
  simp only [Spin.FinPMF.expect, Spin.piPMF_apply, sum_tuple_succ,
    Fin.prod_univ_succ, Fin.cons_zero, Fin.cons_succ, Finset.mul_sum, mul_assoc]

theorem inputMoment_zero (z : ℝ) (xs : Fin 0 → Input) (q : State) :
    inputMoment z xs q = 1 := by
  exact Spin.FinPMF.expect_const _ 1

theorem inputMoment_succ {R : Nat} (z : ℝ) (xs : Fin (R + 1) → Input) (q : State) :
    inputMoment z xs q = emitted z q (xs 0) *
      transvectionLaw.expect (fun t =>
        inputMoment z (Fin.tail xs) (nextState q (xs 0) t)) := by
  unfold inputMoment
  rw [expect_iid_succ]
  simp only [outputWeight, Fin.cons_zero, Fin.tail_cons, pow_add, emitted,
    Spin.FinPMF.expect, Finset.mul_sum]
  congr 1
  funext t
  apply Finset.sum_congr rfl
  intro ts _
  ring

theorem sum_fiber_count {α β : Type*} [Fintype β] [DecidableEq β]
    (s : Finset α) (f : α → β) (g : β → ℝ) :
    (∑ a ∈ s, g (f a)) = ∑ b, ((s.filter (fun a => f a = b)).card : ℝ) * g b := by
  rw [← Finset.sum_fiberwise s f (fun a => g (f a))]
  apply Finset.sum_congr rfl
  intro b _
  rw [Finset.sum_congr rfl (fun a ha => by rw [(Finset.mem_filter.mp ha).2])]
  simp

theorem transvection_expect (q : State) (g : State → ℝ) :
    transvectionLaw.expect (fun t => g (Spin.act t.val.1 t.val.2 q)) =
      ∑ r, Spin.Imt.actLaw q r * g r := by
  classical
  unfold transvectionLaw Spin.FinPMF.expect Spin.FinPMF.uniform Spin.Imt.actLaw
  simp only [Fintype.card_coe]
  simp only [div_mul_eq_mul_div, one_mul, ← Finset.sum_div]
  apply congrArg (fun a : ℝ => a / ((Spin.transPairs 19).card : ℝ))
  rw [Finset.sum_coe_sort (Spin.transPairs 19)
    (fun p : State × State => g (Spin.act p.1 p.2 q))]
  exact sum_fiber_count _ _ _

theorem nextState_expect (q : State) (x : Input) (g : State → ℝ) :
    transvectionLaw.expect (fun t => g (nextState q x t)) =
      ∑ r, Spin.Imt.actLaw q (r ∆ Cset x) * g r := by
  change transvectionLaw.expect (fun t =>
    (fun r => g (r ∆ Cset x)) (Spin.act t.val.1 t.val.2 q)) = _
  rw [transvection_expect q (fun r => g (r ∆ Cset x))]
  apply Fintype.sum_equiv
    (Function.Involutive.toPerm (fun r : State => r ∆ Cset x)
      (fun r => Spin.symmDiff_cancel_tail r (Cset x)))
  intro r
  change Spin.Imt.actLaw q r * g (r ∆ Cset x) =
    Spin.Imt.actLaw q ((r ∆ Cset x) ∆ Cset x) * g (r ∆ Cset x)
  rw [Spin.symmDiff_cancel_tail]

end Spin.Structured.ConcreteEncoder
