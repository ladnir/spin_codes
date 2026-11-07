import Mathlib.Algebra.Polynomial.Roots

/-!
# The RS[8,4] row's information-set property

The old Structured proof first checks that a nonzero base-code message
cannot produce the zero outer word. For Turbo, the corresponding local
fact is stronger: any four distinct evaluation positions determine a
degree-below-four row polynomial. This is the algebraic starting point
for the MDS column-support enumerator.
-/

namespace Spin.Turbo

variable {F : Type*} [Field F]

/-- Four distinct evaluations determine a polynomial of degree below four. -/
theorem rsRow_eq_of_four_evals
    (a : Fin 4 → F) (ha : Function.Injective a)
    (p q : Polynomial F) (hp : p.natDegree < 4) (hq : q.natDegree < 4)
    (heval : ∀ i : Fin 4, p.eval (a i) = q.eval (a i)) : p = q := by
  apply Polynomial.eq_of_natDegree_lt_card_of_eval_eq p q ha heval
  simpa using (max_lt hp hq)

/-- A degree-below-four row with four distinct zero evaluations is zero. -/
theorem rsRow_zero_of_four_zeros
    (a : Fin 4 → F) (ha : Function.Injective a)
    (p : Polynomial F) (hp : p.natDegree < 4)
    (heval : ∀ i : Fin 4, p.eval (a i) = 0) : p = 0 := by
  apply Polynomial.eq_zero_of_natDegree_lt_card_of_eval_eq_zero p ha heval
  simpa using hp

#print axioms rsRow_eq_of_four_evals
#print axioms rsRow_zero_of_four_zeros

end Spin.Turbo
