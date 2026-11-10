import SpinCodes.Turbo.SparseSpectrum22

/-! Audit: RS8SpectrumNumeric authenticates the finite integer MDS counts;
SparseSpectrum22 authenticates actual Cv coefficients. Clear one common
column denominator here, so the final checker can use integers only. -/
noncomputable section
namespace Spin.Turbo
open Finset Polynomial

/-- Numerator polynomial with common denominator (2^32-1)^8. -/
def sparseRS8IntegerPolynomial22 : Polynomial ℤ :=
  ∑ h ∈ Finset.Icc 5 8,
    C (sparseRS8A22 h * (4294967295 : ℤ) ^ (8 - h)) *
      ((1 + 15 * X) ^ 8 - 1) ^ h

theorem sparseRS8IntegerPolynomial22_map :
    sparseRS8IntegerPolynomial22.map (Int.castRingHom ℝ) =
      C ((4294967295 : ℝ) ^ 8) * sparseRS8CountPolynomial22 := by
  apply Polynomial.funext
  intro z
  let a : ℝ := (1 + 15 * z) ^ 8 - 1
  have hi (h : ℕ) :
      ((C (sparseRS8A22 h * (4294967295 : ℤ) ^ (8 - h)) *
        ((1 + 15 * X) ^ 8 - 1) ^ h).map (Int.castRingHom ℝ)).eval z =
        (sparseRS8A22 h : ℝ) * (4294967295 : ℝ) ^ (8 - h) * a ^ h := by
    simp [a]
  have hp : sparseRS8PacketPolynomial22.eval z = a / 4294967295 := by
    norm_num [sparseRS8PacketPolynomial22, a, div_eq_mul_inv]
  unfold sparseRS8IntegerPolynomial22 sparseRS8CountPolynomial22
  rw [Polynomial.map_sum, Polynomial.eval_mul, Polynomial.eval_C]
  simp only [Polynomial.eval_finsetSum]
  simp_rw [hi, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow, hp]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro h hh
  have hh8 : h ≤ 8 := (Finset.mem_Icc.mp hh).2
  have hc : (4294967295 : ℝ) ^ 8 * (a / 4294967295) ^ h =
      (4294967295 : ℝ) ^ (8 - h) * a ^ h := by
    rw [div_pow]
    have he : (4294967295 : ℝ) ^ 8 =
        (4294967295 : ℝ) ^ (8 - h) * (4294967295 : ℝ) ^ h := by
      rw [← pow_add, Nat.sub_add_cancel hh8]
    rw [he]
    field_simp
  calc
    _ = (sparseRS8A22 h : ℝ) * ((4294967295 : ℝ) ^ (8 - h) * a ^ h) := by ring
    _ = (sparseRS8A22 h : ℝ) * ((4294967295 : ℝ) ^ 8 * (a / 4294967295) ^ h) := by rw [hc]
    _ = _ := by ring

theorem sparseRS8ExpectedPacketCount22_eq_integer_coeff
    (points : Fin 8 → Field16) (hp : Function.Injective points) (v : ℕ) :
    rs8ExpectedPacketCount points v =
      (sparseRS8IntegerPolynomial22.coeff v : ℝ) / (4294967295 : ℝ) ^ 8 := by
  have h := congrArg (fun P : Polynomial ℝ => P.coeff v) sparseRS8IntegerPolynomial22_map
  simp only [Polynomial.coeff_map, Int.coe_castRingHom, Polynomial.coeff_C_mul,
    sparseRS8CountPolynomial22_coeff points hp v] at h
  apply (eq_div_iff (by norm_num : (4294967295 : ℝ) ^ 8 ≠ 0)).mpr
  simpa only [mul_comm] using h.symm

#print axioms sparseRS8IntegerPolynomial22_map
#print axioms sparseRS8ExpectedPacketCount22_eq_integer_coeff
end Spin.Turbo
