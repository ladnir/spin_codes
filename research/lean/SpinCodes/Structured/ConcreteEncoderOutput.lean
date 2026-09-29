import SpinCodes.Structured.ConcreteEncoder
import SpinCodes.Structured.ConcreteBinaryEncoding

noncomputable section
namespace Spin.Structured.ConcreteEncoder
open Finset ConcreteMaps ConcreteBinaryEncoding
open scoped symmDiff

/-- The actual emitted blocks, in time order and without a terminal tail. -/
def outputBlocks : {R : ℕ} → (Fin R → Input) → (Fin R → Transvection) → State → (Fin R → Input)
  | 0, _, _, _ => Fin.elim0
  | R+1, xs, ts, q => Fin.cons (xs 0 ∆ Aset q)
      (outputBlocks (Fin.tail xs) (Fin.tail ts) (nextState q (xs 0) (ts 0)))

theorem outputBlocks_weight {R : ℕ} (xs : Fin R → Input) (ts : Fin R → Transvection) (q : State) :
    (∑ i, (outputBlocks xs ts q i).card) = outputWeight xs ts q := by
  induction R generalizing q with
  | zero => simp [outputBlocks, outputWeight]
  | succ R ih =>
    simp only [outputBlocks, Fin.sum_univ_succ, Fin.cons_zero, Fin.cons_succ, outputWeight]
    exact congrArg (_ + ·) (ih _ _ _)

/-- Each input block is recovered from its output once the entering state is known. -/
theorem outputBlocks_injective {R : ℕ} (ts : Fin R → Transvection) (q : State) :
    Function.Injective (fun xs : Fin R → Input => outputBlocks xs ts q) := by
  induction R generalizing q with
  | zero => intro xs ys _; exact Subsingleton.elim _ _
  | succ R ih =>
    intro xs ys he
    have h0 := congrFun he 0
    change xs 0 ∆ Aset q = ys 0 ∆ Aset q at h0
    have hx : xs 0 = ys 0 := by
      have h := congrArg (fun s => s ∆ Aset q) h0
      simpa only [symmDiff_symmDiff_cancel_right] using h
    have ht := congrArg Fin.tail he
    change outputBlocks (Fin.tail xs) (Fin.tail ts) (nextState q (xs 0) (ts 0)) =
      outputBlocks (Fin.tail ys) (Fin.tail ts) (nextState q (ys 0) (ts 0)) at ht
    rw [hx] at ht
    have htail := ih (Fin.tail ts) _ ht
    funext i
    exact Fin.cases hx (fun j => congrFun htail j) i

theorem act_xor {n : ℕ} (u v q r : Finset (Fin n)) :
    Spin.act u v (q ∆ r) = Spin.act u v q ∆ Spin.act u v r := by
  have hp : Even ((q ∆ r) ∩ v).card ↔ (Even (q ∩ v).card ↔ Even (r ∩ v).card) := by
    rw [Spin.symmDiff_inter_left]
    have hc := Spin.card_symmDiff_add_two_mul_card_inter (q ∩ v) (r ∩ v)
    simp only [even_iff_two_dvd, Nat.dvd_iff_mod_eq_zero]
    omega
  unfold Spin.act
  by_cases hq : Even (q ∩ v).card <;> by_cases hr : Even (r ∩ v).card <;>
    simp only [hp, hq, hr, iff_self, true_iff, false_iff, not_true_eq_false,
      ite_true, ite_false]
  all_goals ext i; simp only [mem_symmDiff]; tauto

theorem nextState_xor (q r : State) (x y : Input) (t : Transvection) :
    nextState (q ∆ r) (x ∆ y) t = nextState q x t ∆ nextState r y t := by
  simp only [nextState, act_xor, Cset_xor]
  ac_rfl

/-- For each fixed transvection sequence, the emitted block map is binary linear. -/
theorem outputBlocks_xor {R : ℕ} (xs ys : Fin R → Input) (ts : Fin R → Transvection) (q r : State) :
    outputBlocks (fun i => xs i ∆ ys i) ts (q ∆ r) =
      fun i => outputBlocks xs ts q i ∆ outputBlocks ys ts r i := by
  induction R generalizing q r with
  | zero => exact Subsingleton.elim _ _
  | succ R ih =>
    funext i
    refine Fin.cases ?_ (fun j => ?_) i
    · simp only [outputBlocks, Fin.cons_zero, Aset_xor]
      ac_rfl
    · simp only [outputBlocks, Fin.cons_succ, nextState_xor]
      exact congrFun (ih (Fin.tail xs) (Fin.tail ys) (Fin.tail ts)
        (nextState q (xs 0) (ts 0)) (nextState r (ys 0) (ts 0))) j

def outputWord {R : ℕ} (xs : Fin R → Input) (ts : Fin R → Transvection) (q : State := ∅) :
    Fin (R*128) → ZMod 2 := rowsToWord (outputBlocks xs ts q)

def binaryWeight {n : ℕ} (word : Fin n → ZMod 2) : ℕ := (univ.filter (fun i => word i ≠ 0)).card

theorem rowsToWord_weight {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    binaryWeight (rowsToWord rows) = ∑ i, (rows i).card := by
  unfold binaryWeight
  rw [card_eq_sum_ones, sum_filter,
    ← Equiv.sum_comp (finProdFinEquiv : Fin L × Fin b ≃ Fin (L*b))]
  simp only [Fintype.sum_prod_type, rowsToWord, Equiv.symm_apply_apply]
  apply sum_congr rfl
  intro i _
  simp [bitEquiv]

theorem outputWord_weight {R : ℕ} (xs : Fin R → Input) (ts : Fin R → Transvection) (q : State) :
    binaryWeight (outputWord xs ts q) = outputWeight xs ts q := by
  rw [outputWord, rowsToWord_weight, outputBlocks_weight]

theorem outputWord_injective {R : ℕ} (ts : Fin R → Transvection) (q : State) :
    Function.Injective (fun xs : Fin R → Input => outputWord xs ts q) := by
  intro xs ys h
  exact outputBlocks_injective ts q ((rowsEquiv R 128).injective h)

end Spin.Structured.ConcreteEncoder
