import SpinCodes.Turbo.MatrixCoefficient22
import Mathlib.Algebra.Polynomial.Coeff

/-!
# Polynomial identity for the ordered regional recurrence

The external producer uses truncated binary polynomial powering. Its
coefficients live in a noncommutative matrix semiring. This bridge identifies
the finite ordered-choice recurrence with polynomial coefficients, so that
the associativity used by binary powering does not reorder physical steps.
-/
noncomputable section
namespace Spin.Turbo
open Finset Polynomial

variable {A I : Type*} [Fintype A] [Fintype I] [DecidableEq I]

def orderedMatrixPolynomial (degree : A → ℕ) (W : A → Matrix I I ℝ) :
    Polynomial (Matrix I I ℝ) := ∑ a : A, Polynomial.monomial (degree a) (W a)

private theorem coeff_monomial_mul_all (P : Polynomial (Matrix I I ℝ))
    (d q : ℕ) (M : Matrix I I ℝ) :
    (Polynomial.monomial d M * P).coeff q =
      if d ≤ q then M * P.coeff (q - d) else 0 := by
  rw [← Polynomial.C_mul_X_pow_eq_monomial, mul_assoc,
    Polynomial.coeff_C_mul, Polynomial.coeff_X_pow_mul']
  simp only [mul_ite, mul_zero]

theorem orderedMatrixPolynomial_pow_coeff (degree : A → ℕ)
    (W : A → Matrix I I ℝ) (n q : ℕ) :
    ((orderedMatrixPolynomial degree W) ^ n).coeff q =
      orderedMatrixDP degree W n q := by
  induction n generalizing q with
  | zero => simp only [pow_zero, Polynomial.coeff_one, orderedMatrixDP]
  | succ n ih =>
    rw [pow_succ']
    conv_lhs => arg 1; arg 1; unfold orderedMatrixPolynomial
    rw [Finset.sum_mul, Polynomial.finsetSum_coeff]
    change (∑ a : A, (Polynomial.monomial (degree a) (W a) *
      (orderedMatrixPolynomial degree W) ^ n).coeff q) = _
    simp only [coeff_monomial_mul_all, ih, orderedMatrixDP]

/-- Coefficients of ordered binary products: this is the convolution used
by the producer. Terms above q cannot affect coefficient q. -/
theorem orderedMatrixDP_add (degree : A → ℕ)
    (W : A → Matrix I I ℝ) (m n q : ℕ) :
    orderedMatrixDP degree W (m + n) q =
      ∑ ab ∈ Finset.antidiagonal q,
        orderedMatrixDP degree W m ab.1 * orderedMatrixDP degree W n ab.2 := by
  simp only [← orderedMatrixPolynomial_pow_coeff]
  rw [pow_add, Polynomial.coeff_mul]

#print axioms orderedMatrixPolynomial_pow_coeff
#print axioms orderedMatrixDP_add
end Spin.Turbo
