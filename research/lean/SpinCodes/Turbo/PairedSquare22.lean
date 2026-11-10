import SpinCodes.Turbo.SparseRegionalPower22
import SpinCodes.Turbo.ConditionalFractionalMatrix

/-! At fractional exponent α = 1 the paired two-step profile polynomial is the
square of the single-step polynomial, so the paired coefficient is the
coefficient of the single-step power with twice as many steps. -/
noncomputable section
namespace Spin.Turbo
open Finset Polynomial

namespace PairedSquare22

theorem pairedProfileMatrix22_one (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (ab : Fin 17 × Fin 17) :
    pairedProfileMatrix22 B 1 ab =
      ((Nat.choose 16 ab.1.val : ℝ) • B ab.1) * ((Nat.choose 16 ab.2.val : ℝ) • B ab.2) := by
  unfold pairedProfileMatrix22
  have h : matrixEntrywiseRpow (B ab.1 * B ab.2) 1 = B ab.1 * B ab.2 := by
    funext i j
    simp [matrixEntrywiseRpow, Real.rpow_one]
  rw [h, smul_mul_smul_comm]

theorem pairedPoly_eq_sq (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) :
    orderedMatrixPolynomial pairedProfileDegree22 (pairedProfileMatrix22 B 1) =
      (orderedMatrixPolynomial Fin.val (fun j : Fin 17 => (Nat.choose 16 j.val : ℝ) • B j)) ^ 2 := by
  unfold orderedMatrixPolynomial
  rw [sq, Finset.sum_mul_sum, ← Finset.univ_product_univ, Finset.sum_product]
  apply Finset.sum_congr rfl
  intro a _
  apply Finset.sum_congr rfl
  intro b _
  rw [Polynomial.monomial_mul_monomial, pairedProfileMatrix22_one]
  rfl

theorem pairedDP_one_eq (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (n q : ℕ) :
    orderedMatrixDP pairedProfileDegree22 (pairedProfileMatrix22 B 1) n q =
      regionalCountDP22 B (2 * n) q := by
  unfold regionalCountDP22
  rw [← orderedMatrixPolynomial_pow_coeff, ← orderedMatrixPolynomial_pow_coeff, pairedPoly_eq_sq,
    ← pow_mul]

theorem pairedRegionDP22_one_eq (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ) (n q : ℕ) :
    pairedRegionDP22 B 1 n q =
      (1 / (Nat.choose (32 * n) q : ℝ)) • regionalCountDP22 B (2 * n) q := by
  unfold pairedRegionDP22
  rw [pairedDP_one_eq]

end PairedSquare22

#print axioms PairedSquare22.pairedRegionDP22_one_eq
end Spin.Turbo
