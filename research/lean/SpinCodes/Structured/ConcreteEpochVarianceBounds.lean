import SpinCodes.Structured.ConcreteEpochVariance

/-! Variance and absolute deviation for the actual empty-epoch emitted sum. -/
noncomputable section
namespace Spin.Structured.ConcreteEncoder

theorem expect_sq_sub {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → ℝ) (c : ℝ) :
    P.expect (fun x => (f x - c) ^ 2) =
      P.expect (fun x => f x ^ 2) - 2 * c * P.expect f + c ^ 2 := by
  have he (x) : f x - c = -c + f x := by ring
  simp only [he, expect_sq_shift]
  ring

def finiteVariance {α : Type*} [Fintype α] (P : Spin.FinPMF α) (f : α → ℝ) : ℝ :=
  P.expect (fun x => (f x - P.expect f) ^ 2)

theorem finiteVariance_le_centered {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → ℝ) (c : ℝ) :
    finiteVariance P f ≤ P.expect (fun x => (f x - c) ^ 2) := by
  rw [finiteVariance, expect_sq_sub, expect_sq_sub]
  nlinarith [sq_nonneg (P.expect f - c)]

theorem expect_abs_le_sqrt_square {α : Type*} [Fintype α] (P : Spin.FinPMF α)
    (f : α → ℝ) :
    P.expect (fun x => |f x|) ≤ Real.sqrt (P.expect (fun x => f x ^ 2)) := by
  have hn := P.expect_nonneg (fun x => sq_nonneg (|f x| - P.expect (fun x => |f x|)))
  rw [expect_sq_sub] at hn
  simp only [sq_abs] at hn
  apply Real.le_sqrt_of_sq_le
  nlinarith

def emptyVariance (g : Nat) (q : State) : ℝ :=
  finiteVariance (Spin.piPMF (fun _ : Fin g => transvectionLaw))
    (fun ts => (outputWeight (fun _ => ∅) ts q : ℝ))

theorem emptyVariance_bound (g : Nat) {q : State} (hq : q ≠ ∅) :
    emptyVariance g q ≤ 3 * 128 ^ 2 * g := by
  exact (finiteVariance_le_centered _ _ (epochMean * g)).trans (emptyCenteredSecond_bound g hq)

def emptyAbsDeviation (g : Nat) (q : State) : ℝ :=
  (Spin.piPMF (fun _ : Fin g => transvectionLaw)).expect
    (fun ts => |(outputWeight (fun _ => ∅) ts q : ℝ) - epochMean * g|)

theorem emptyAbsDeviation_bound (g : Nat) {q : State} (hq : q ≠ ∅) :
    emptyAbsDeviation g q ≤ 128 * Real.sqrt (3 * g) := by
  calc emptyAbsDeviation g q ≤ Real.sqrt (emptyCenteredSecond g q) :=
      expect_abs_le_sqrt_square _ _
    _ ≤ Real.sqrt (3 * 128 ^ 2 * (g : ℝ)) :=
      Real.sqrt_le_sqrt (emptyCenteredSecond_bound g hq)
    _ = _ := by
      rw [show (3 : ℝ) * 128 ^ 2 * g = 128 ^ 2 * (3 * g) by ring,
        Real.sqrt_mul (sq_nonneg (128 : ℝ)), Real.sqrt_sq_eq_abs]
      norm_num

end Spin.Structured.ConcreteEncoder
