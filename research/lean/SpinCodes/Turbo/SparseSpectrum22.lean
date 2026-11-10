import SpinCodes.Turbo.RS8FieldMoment
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Algebra.Polynomial.Coeff

/-! Exact current RS8 Cv coefficient interface. The historical spectrum
proof authenticates the MDS Ah counts and the generating function of the
actual shared setup's expected nonzero-message measure. Polynomial
extensionality converts that checked real identity into coefficientwise
arithmetic targets, including every shell0..64. -/
noncomputable section
namespace Spin.Turbo
open Finset Polynomial

def sparseRS8A22 (h : ℕ) : ℤ :=
  (Nat.choose 8 h : ℤ) * ∑ j ∈ Finset.range (h - 4),
    (-1 : ℤ) ^ j * (h.choose j : ℤ) *
      (((2 ^ 32 : ℕ) ^ (h - 4 - j) - 1 : ℕ) : ℤ)

def sparseRS8PacketPolynomial22 : Polynomial ℝ :=
  ((1 + 15 * X) ^ 8 - 1) * C (1 / (2 ^ 32 - 1 : ℝ))

def sparseRS8CountPolynomial22 : Polynomial ℝ :=
  ∑ h ∈ Finset.Icc 5 8, C (sparseRS8A22 h : ℝ) * sparseRS8PacketPolynomial22 ^ h

def sparseRS8ExpectedPolynomial22 (points : Fin 8 → Field16) : Polynomial ℝ :=
  ∑ v ∈ Finset.Icc 5 64, Polynomial.monomial v (rs8ExpectedPacketCount points v)

theorem sparseRS8CountPolynomial22_eq_expected
    (points : Fin 8 → Field16) (hp : Function.Injective points) :
    sparseRS8CountPolynomial22 = sparseRS8ExpectedPolynomial22 points := by
  apply Polynomial.funext
  intro z
  have h := rs8_expected_packet_count_gf points hp z
  simpa only [sparseRS8CountPolynomial22, sparseRS8ExpectedPolynomial22,
    sparseRS8PacketPolynomial22, sparseRS8A22, Polynomial.eval_finsetSum,
    Polynomial.eval_monomial, Polynomial.eval_mul, Polynomial.eval_C,
    Polynomial.eval_pow, Polynomial.eval_sub, Polynomial.eval_add,
    Polynomial.eval_one, Polynomial.eval_ofNat, Polynomial.eval_X,
    div_eq_mul_inv, one_mul] using h.symm

theorem rs8ExpectedPacketCount_eq_zero_of_not_range22
    (points : Fin 8 → Field16) (hp : Function.Injective points)
    (v : ℕ) (hv : v ∉ Finset.Icc 5 64) : rs8ExpectedPacketCount points v = 0 := by
  unfold rs8ExpectedPacketCount
  apply Finset.sum_eq_zero
  intro msg hm
  have hm0 : msg ≠ 0 := (Finset.mem_filter.mp hm).2
  unfold FinPMF.prob
  apply Finset.sum_eq_zero
  intro coeff hc
  have he := (Finset.mem_filter.mp hc).2
  have hr := rs8_sampled_packet_weight_range points hp msg hm0 coeff
  have hmem : rs8SampledPacketWeight points msg coeff ∈ Finset.Icc 5 64 :=
    Finset.mem_Icc.mpr hr
  exact False.elim (hv (he ▸ hmem))

theorem sparseRS8CountPolynomial22_coeff
    (points : Fin 8 → Field16) (hp : Function.Injective points) (v : ℕ) :
    sparseRS8CountPolynomial22.coeff v = rs8ExpectedPacketCount points v := by
  rw [sparseRS8CountPolynomial22_eq_expected points hp]
  unfold sparseRS8ExpectedPolynomial22
  rw [Polynomial.finsetSum_coeff]
  simp only [Polynomial.coeff_monomial]
  by_cases hv : v ∈ Finset.Icc 5 64
  · simp only [Finset.sum_ite_eq', hv, ite_true]
  · rw [rs8ExpectedPacketCount_eq_zero_of_not_range22 points hp v hv]
    simp only [Finset.sum_ite_eq', hv, ite_false]

/-- Rational polynomial used by the exact certificate arithmetic. -/
def sparseRS8RationalPolynomial22 : Polynomial ℚ :=
  ∑ h ∈ Finset.Icc 5 8, C (sparseRS8A22 h : ℚ) *
    (((1 + 15 * X) ^ 8 - 1) * C (1 / (2 ^ 32 - 1 : ℚ))) ^ h

theorem sparseRS8RationalPolynomial22_map :
    sparseRS8RationalPolynomial22.map (Rat.castHom ℝ) = sparseRS8CountPolynomial22 := by
  simp [Polynomial.map_sum, sparseRS8RationalPolynomial22, sparseRS8CountPolynomial22,
    sparseRS8PacketPolynomial22]

theorem sparseRS8ExpectedPacketCount22_eq_rational_coeff
    (points : Fin 8 → Field16) (hp : Function.Injective points) (v : ℕ) :
    rs8ExpectedPacketCount points v = (sparseRS8RationalPolynomial22.coeff v : ℝ) := by
  rw [← sparseRS8CountPolynomial22_coeff points hp v, ← sparseRS8RationalPolynomial22_map]
  simp only [Polynomial.coeff_map, Rat.coe_castHom]
#print axioms sparseRS8CountPolynomial22_eq_expected
#print axioms sparseRS8ExpectedPacketCount22_eq_rational_coeff
#print axioms sparseRS8CountPolynomial22_coeff
end Spin.Turbo



