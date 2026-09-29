import SpinCodes.Structured.ConcreteEpoch

/-! Actual endpoint kernels for an empty epoch, with exact structural zeros
and the geometric endpoint mixing estimate. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder
open ConcreteMaps
open scoped symmDiff

def emptyKernel (z : ℝ) (g : Nat) (q r : State) : ℝ :=
  (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect (fun ts =>
    z ^ outputWeight (fun _ => ∅) ts q * (if emptyState ts q = r then 1 else 0))

theorem liveMean_point (r : State) :
    liveMean (fun q => if q = r then 1 else 0) = if r = ∅ then 0 else 1 / 524287 := by
  unfold liveMean
  rw [Finset.sum_ite_eq']
  simp only [Spin.nonzeroStates, Finset.mem_filter, Finset.mem_univ, true_and]
  by_cases hr : r = ∅ <;> simp [hr]

theorem emptyKernel_one (g : Nat) (q r : State) : emptyKernel 1 g q r =
    emptyExpect g (fun w => if w = r then 1 else 0) q := by
  simp only [emptyKernel, emptyExpect, one_pow, one_mul]

/-- The full endpoint law, with both the point mass and uniform-live parts. -/
theorem emptyKernel_one_live (g : Nat) {q : State} (hq : q ≠ ∅) (r : State) :
    emptyKernel 1 g q r = (1 / 2 : ℝ) ^ g * (if q = r then 1 else 0) +
      (1 - (1 / 2 : ℝ) ^ g) * (if r = ∅ then 0 else 1 / 524287) := by
  rw [emptyKernel_one, emptyExpect_live _ _ hq, liveMean_point]
  ring

theorem emptyKernel_one_error (g : Nat) {q r : State} (hq : q ≠ ∅) (hr : r ≠ ∅) :
    |emptyKernel 1 g q r - 1 / 524287| ≤ (1 / 2 : ℝ) ^ g := by
  rw [emptyKernel_one_live g hq r, if_neg hr, abs_le]
  have hp : 0 ≤ (1 / 2 : ℝ) ^ g := by positivity
  by_cases he : q = r <;> simp only [he, ite_true, ite_false] <;>
    constructor <;> nlinarith

theorem emptyKernel_nonneg {z : ℝ} (hz : 0 ≤ z) (g : Nat) (q r : State) :
    0 ≤ emptyKernel z g q r := by
  apply Spin.FinPMF.expect_nonneg
  intro ts
  split_ifs <;> positivity

theorem emptyKernel_live_zero (z : ℝ) (g : Nat) {q : State} (hq : q ≠ ∅) :
    emptyKernel z g q ∅ = 0 := by
  unfold emptyKernel
  simp only [if_neg (emptyState_nonzero _ hq), mul_zero, Spin.FinPMF.expect_const]

theorem empty_output_zero {g : Nat} (ts : Fin g → Transvection) :
    outputWeight (fun _ => ∅) ts ∅ = 0 := by
  induction g with
  | zero => rfl
  | succ g ih =>
    simp only [outputWeight, Aset_empty, nextState_empty, Spin.act_empty]
    rw [show (∅ : Input) ∆ ∅ = ∅ from symmDiff_self ∅, Finset.card_empty, zero_add]
    exact ih _

theorem emptyKernel_zero (z : ℝ) (g : Nat) (r : State) :
    emptyKernel z g ∅ r = if r = ∅ then 1 else 0 := by
  unfold emptyKernel
  simp only [empty_output_zero, pow_zero, one_mul, emptyState_zero]
  rw [Spin.FinPMF.expect_const]
  exact if_congr eq_comm rfl rfl

theorem emptyKernel_total (z : ℝ) (g : Nat) (q : State) :
    ∑ r, emptyKernel z g q r = inputMoment z (fun _ : Fin g => ∅) q := by
  unfold emptyKernel inputMoment
  rw [← Spin.FinPMF.expect_sum]
  apply congrArg (Spin.FinPMF.expect _)
  funext ts
  rw [← Finset.mul_sum]
  simp

end Spin.Structured.ConcreteEncoder
