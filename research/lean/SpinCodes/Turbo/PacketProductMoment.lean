import SpinCodes.Turbo.PacketSupportCount
import SpinCodes.FiniteProductLaw

/-!
# Product moments for independent Turbo columns

The old finite proof uses `FiniteProductLaw.prod_sum_eq_sum_prod` to keep
independent setup choices as a product law. The same identity turns the
weight generating function of eight independently randomized columns
into a product of their local generating functions.
-/

noncomputable section

namespace Spin.Turbo

open Finset

theorem piPMF_expect_prod
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    {α : ι → Type*} [∀ i, Fintype (α i)]
    (P : ∀ i, FinPMF (α i)) (f : ∀ i, α i → ℝ) :
    (piPMF P).expect (fun x => ∏ i, f i (x i)) =
      ∏ i, (P i).expect (f i) := by
  calc
    _ = ∑ x : ∀ i, α i, ∏ i, (P i).p (x i) * f i (x i) := by
      simp [FinPMF.expect, piPMF_apply, Finset.prod_mul_distrib]
    _ = ∏ i, ∑ a : α i, (P i).p a * f i a :=
      (prod_sum_eq_sum_prod (fun i a => (P i).p a * f i a)).symm
    _ = _ := rfl

#print axioms piPMF_expect_prod

end Spin.Turbo
