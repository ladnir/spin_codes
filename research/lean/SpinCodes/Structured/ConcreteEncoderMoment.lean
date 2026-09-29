import SpinCodes.Structured.ConcreteEncoder

/-! Averaging the explicit encoder agrees with the checked occupation kernel. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps Spin.Imt
open scoped symmDiff

def averagedMoment (P : Spin.FinPMF Input) (R : Nat) (z : ℝ) (q : State) : ℝ :=
  (Spin.piPMF (fun _ : Fin R => P)).expect (fun xs => inputMoment z xs q)

/-- The explicit sample space: independent inputs and independent valid transvections. -/
def experimentLaw (P : Spin.FinPMF Input) (R : Nat) :
    Spin.FinPMF ((Fin R → Input) × (Fin R → Transvection)) :=
  (Spin.piPMF (fun _ : Fin R => P)).prod
    (Spin.piPMF (fun _ : Fin R => transvectionLaw))

theorem experiment_expect (P : Spin.FinPMF Input) (R : Nat) (z : ℝ) (q : State) :
    (experimentLaw P R).expect (fun ω => z ^ outputWeight ω.1 ω.2 q) =
      averagedMoment P R z q := by
  exact Spin.FinPMF.expect_prod _ _ _

theorem expect_scale_sum {α ι : Type*} [Fintype α] [Fintype ι]
    (P : Spin.FinPMF α) (a : ℝ) (c : ι → ℝ) (f : α → ι → ℝ) :
    P.expect (fun x => a * ∑ i, c i * f x i) =
      a * ∑ i, c i * P.expect (fun x => f x i) := by
  simp only [Spin.FinPMF.expect, Finset.mul_sum]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro i _
  apply Finset.sum_congr rfl
  intro x _
  ring

theorem averagedMoment_zero (P : Spin.FinPMF Input) (z : ℝ) (q : State) :
    averagedMoment P 0 z q = 1 := by
  unfold averagedMoment
  simp only [inputMoment_zero]
  exact Spin.FinPMF.expect_const _ 1

theorem averagedMoment_succ (P : Spin.FinPMF Input) (R : Nat) (z : ℝ) (q : State) :
    averagedMoment P (R + 1) z q =
      ∑ x, P.p x * emitted z q x *
        ∑ r, actLaw q (r ∆ Cset x) * averagedMoment P R z r := by
  unfold averagedMoment
  rw [expect_iid_succ]
  simp only [inputMoment_succ, Fin.cons_zero, Fin.tail_cons]
  have ht (x : Input) (xs : Fin R → Input) :
      transvectionLaw.expect (fun t => inputMoment z xs (nextState q x t)) =
        ∑ r, actLaw q (r ∆ Cset x) * inputMoment z xs r :=
    nextState_expect q x (fun r => inputMoment z xs r)
  simp_rw [ht, expect_scale_sum]
  simp only [Spin.FinPMF.expect]
  apply Finset.sum_congr rfl
  intro x _
  ring

theorem averagedMoment_bernoulli_succ (P : Spin.FinPMF Input) (β z : ℝ)
    (hP : ∀ x, P.p x = β ^ x.card * (1 - β) ^ (128 - x.card))
    (R : Nat) (q : State) :
    averagedMoment P (R + 1) z q =
      ∑ r, bernoulliRow β z q r * averagedMoment P R z r := by
  rw [averagedMoment_succ]
  simp only [hP, bernoulliRow, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro r _
  apply Finset.sum_congr rfl
  intro x _
  ring

theorem averagedMoment_kernel (P : Spin.FinPMF Input) (β z : ℝ)
    (hP : ∀ x, P.p x = β ^ x.card * (1 - β) ^ (128 - x.card))
    (R : Nat) (μ : State → ℝ) :
    (∑ q, μ q * averagedMoment P R z q) =
      ∑ q, ((bernoulliStep β z)^[R] μ) q := by
  induction R generalizing μ with
  | zero => simp only [averagedMoment_zero, mul_one, Function.iterate_zero, id_eq]
  | succ R ih =>
    simp only [averagedMoment_bernoulli_succ P β z hP R, Finset.mul_sum]
    rw [Finset.sum_comm]
    have he : (∑ r, ∑ q, μ q * (bernoulliRow β z q r * averagedMoment P R z r)) =
        ∑ r, (bernoulliStep β z μ) r * averagedMoment P R z r := by
      simp only [bernoulliStep, kernelApply, Finset.sum_mul, mul_assoc]
    rw [he, ih]
    rw [Function.iterate_succ_apply]

theorem averagedMoment_eq_iterate (P : Spin.FinPMF Input) (β z : ℝ)
    (hP : ∀ x, P.p x = β ^ x.card * (1 - β) ^ (128 - x.card)) (R : Nat) :
    (Spin.piPMF (fun _ : Fin R => P)).expect (fun xs => inputMoment z xs ∅) =
      ∑ q, ((bernoulliStep β z)^[R] (fun q => if q = ∅ then (1 : ℝ) else 0)) q := by
  simpa only [averagedMoment, ite_mul, one_mul, zero_mul, Finset.sum_ite_eq',
    Finset.mem_univ, ite_true]
    using averagedMoment_kernel P β z hP R (fun q => if q = ∅ then (1 : ℝ) else 0)

theorem encoder_moment_bound (P : Spin.FinPMF Input) {β z : ℝ}
    (hP : ∀ x, P.p x = β ^ x.card * (1 - β) ^ (128 - x.card))
    (hβ0 : 0 ≤ β) (hβ1 : β ≤ 1) (hz : 0 ≤ z) (hz1 : z ≤ 1) (R : Nat) :
    (Spin.piPMF (fun _ : Fin R => P)).expect (fun xs => inputMoment z xs ∅) ≤
      (((Occupation.Sparse.numericalMatrix β z).apply)^[R] (Coords.eZ 5)).total := by
  rw [averagedMoment_eq_iterate P β z hP]
  exact bernoulli_moment_bound hβ0 hβ1 hz hz1 R

theorem sparse_encoder_moment (P : Spin.FinPMF Input) {α : ℝ}
    (hP : ∀ x, P.p x = ((4 / 5) * α) ^ x.card *
      (1 - (4 / 5) * α) ^ (128 - x.card))
    (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : Nat) :
    (Spin.piPMF (fun _ : Fin R => P)).expect
        (fun xs => inputMoment (1 - (8 / 5) * α) xs ∅) ≤
      2048 * (1 - 96 * α) ^ R := by
  rw [averagedMoment_eq_iterate P ((4 / 5) * α) (1 - (8 / 5) * α) hP]
  exact sparse_bernoulli_moment h0 h1 R

theorem sparse_encoder_experiment (P : Spin.FinPMF Input) {α : ℝ}
    (hP : ∀ x, P.p x = ((4 / 5) * α) ^ x.card *
      (1 - (4 / 5) * α) ^ (128 - x.card))
    (h0 : 0 < α) (h1 : α ≤ 1 / 10000) (R : Nat) :
    (experimentLaw P R).expect
        (fun ω => (1 - (8 / 5) * α) ^ outputWeight ω.1 ω.2 ∅) ≤
      2048 * (1 - 96 * α) ^ R := by
  rw [experiment_expect]
  exact sparse_encoder_moment P hP h0 h1 R

end Spin.Structured.ConcreteEncoder
