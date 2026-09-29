import SpinCodes.Structured.DenseNativeRemainder

noncomputable section
namespace Spin.Structured.DenseOccupationFixed
open Finset

/-- All finite counting/routing prefactors fit the vanishing dense remainder. -/
theorem native_occupation_rate (m Q : ℕ) (hQL : Q ≤ Lsched m) {C : ℝ} (hC : 0<C) (η : ℝ) :
    densePrefactor (Lsched m) (bsched m) Q C * Real.exp (-η*(Nsched m : ℝ)) ≤
      Real.exp (-η*(Nsched m : ℝ)+denseRemainder C m*(Nsched m : ℝ)) := by
  calc
    _ ≤ Real.exp (denseRemainder C m*(Nsched m : ℝ))*Real.exp (-η*(Nsched m : ℝ)) :=
      mul_le_mul_of_nonneg_right (native_densePrefactor_le_exp m hQL hC) (Real.exp_pos _).le
    _ = _ := by rw [← Real.exp_add]; congr 1; ring

/-- Summing the positive-occupation range pays at most one factor L. -/
theorem native_dense_sum (m cut : ℕ) (EZ : ℕ → ℝ) {C : ℝ} (hC : 0<C) (η : ℝ)
    (h : ∀ Q ∈ Ico (cut+1) (Lsched m+1), EZ Q ≤
      densePrefactor (Lsched m) (bsched m) Q C * Real.exp (-η*(Nsched m : ℝ))) :
    ∑ Q ∈ Ico (cut+1) (Lsched m+1), EZ Q ≤
      (Lsched m : ℝ)*Real.exp (-η*(Nsched m : ℝ)+denseRemainder C m*(Nsched m : ℝ)) := by
  have he : ∀ Q ∈ Ico (cut+1) (Lsched m+1), EZ Q ≤
      Real.exp (-η*(Nsched m : ℝ)+denseRemainder C m*(Nsched m : ℝ)) := by
    intro Q hQ
    have hQL : Q ≤ Lsched m := by have := (mem_Ico.mp hQ).2; omega
    exact (h Q hQ).trans (native_occupation_rate m Q hQL hC η)
  calc
    _ ≤ ∑ _Q ∈ Ico (cut+1) (Lsched m+1),
        Real.exp (-η*(Nsched m : ℝ)+denseRemainder C m*(Nsched m : ℝ)) := sum_le_sum he
    _ = ((Ico (cut+1) (Lsched m+1)).card : ℝ)*
        Real.exp (-η*(Nsched m : ℝ)+denseRemainder C m*(Nsched m : ℝ)) := by
      rw [sum_const, nsmul_eq_mul]
    _ ≤ _ := by
      apply mul_le_mul_of_nonneg_right _ (Real.exp_pos _).le
      exact_mod_cast (show (Ico (cut+1) (Lsched m+1)).card ≤ Lsched m by rw [Nat.card_Ico]; omega)

#print axioms native_occupation_rate
#print axioms native_dense_sum

end Spin.Structured.DenseOccupationFixed
