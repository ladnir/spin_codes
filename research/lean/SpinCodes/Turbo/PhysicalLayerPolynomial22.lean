import SpinCodes.Turbo.PhysicalStep22
import Mathlib.Algebra.Polynomial.Coeff

/-!
# Exact-layer normalization for arbitrary packet observables

This is the coefficient form of `packetStepLaw22_apply`. It is independent
of the inner output map and also applies to signed character observables.
-/
noncomputable section
namespace Spin.Turbo
open Finset Polynomial

def degreePolynomial22 {Ω : Type*} [Fintype Ω] (degree : Ω → ℕ)
    (w : Ω → ℝ) : Polynomial ℝ :=
  ∑ x, Polynomial.monomial (degree x) (w x)

private theorem finiteLayerCoefficient {Ω : Type*} (s : Finset Ω)
    (p w : Ω → ℝ) (degree : Ω → ℕ) (j : ℕ) (D : ℝ)
    (hp : ∀ x, p x = if degree x = j then 1/D else 0) :
    (∑ x ∈ s, p x * w x) =
      (∑ x ∈ s, Polynomial.monomial (degree x) (w x)).coeff j / D := by
  simp only [hp, Polynomial.finsetSum_coeff, Polynomial.coeff_monomial]
  rw [Finset.sum_div]
  apply Finset.sum_congr rfl
  intro x _
  by_cases hx : degree x = j
  · simp only [hx, ite_true, div_eq_mul_inv, one_mul]
    exact mul_comm (D⁻¹) (w x)
  · simp only [hx, ite_false, zero_mul, zero_div]

private theorem finiteLawCoefficient {Ω : Type*} [Fintype Ω]
    (P : FinPMF Ω) (w : Ω → ℝ) (degree : Ω → ℕ) (j : ℕ) (D : ℝ)
    (hp : ∀ x, P.p x = if degree x = j then 1/D else 0) :
    P.expect w = (degreePolynomial22 degree w).coeff j / D :=
  finiteLayerCoefficient Finset.univ P.p w degree j D hp

theorem packetStepLaw22_expect_layer (j : Fin 17) (w : PacketStep16 → ℝ) :
    (packetStepLaw22 j).expect w =
      ∑ x : PacketStep16, if (packetStepSupport22 x).card = j.val then
        w x / ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val) else 0 := by
  unfold FinPMF.expect
  apply Finset.sum_congr rfl
  intro x _
  rw [packetStepLaw22_apply]
  by_cases hx : (packetStepSupport22 x).card = j.val
  · simp only [hx, ite_true, div_eq_mul_inv, one_mul]
    exact mul_comm (((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val)⁻¹) (w x)
  · simp only [hx, ite_false, zero_mul]

theorem degreePolynomial22_coeff_div {Ω : Type*} [Fintype Ω]
    (degree : Ω → ℕ) (w : Ω → ℝ) (j : ℕ) (D : ℝ) :
    (degreePolynomial22 degree w).coeff j / D =
      ∑ x, if degree x = j then w x / D else 0 := by
  unfold degreePolynomial22
  rw [Polynomial.finsetSum_coeff, Finset.sum_div]
  apply Finset.sum_congr rfl
  intro x _
  rw [Polynomial.coeff_monomial]
  by_cases hx : degree x = j
  · simp only [hx, ite_true]
  · simp only [hx, ite_false, zero_div]

theorem packetStepLaw22_expect_coefficient (j : Fin 17) (w : PacketStep16 → ℝ) :
    (packetStepLaw22 j).expect w =
      (degreePolynomial22 (fun x => (packetStepSupport22 x).card) w).coeff j.val /
        ((Nat.choose 16 j.val : ℝ) * (15 : ℝ)^j.val) := by
  rw [packetStepLaw22_expect_layer, degreePolynomial22_coeff_div]

#print axioms packetStepLaw22_expect_layer
#print axioms packetStepLaw22_expect_coefficient
end Spin.Turbo
