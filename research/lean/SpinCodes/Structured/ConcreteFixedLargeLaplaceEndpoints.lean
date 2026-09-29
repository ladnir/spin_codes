import SpinCodes.Structured.ConcreteFixedLargeGapSymmetry
import SpinCodes.Structured.ConcreteFixedLargeExponent

/-! Normalized beta and deterministic endpoint bounds for the actual spacing law. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Set MeasureTheory

theorem orderedSite_normalized_const (Q : ℕ) (c : ℝ) :
    (Q.factorial:ℝ) * (∫ _x in orderedSiteDomain Q, c) = c := by
  have h : (∫ _x in orderedSiteDomain Q, c) = c * (∫ _x in orderedSiteDomain Q, (1:ℝ)) := by
    rw [← integral_const_mul]
    simp
  rw [h, ← mul_assoc, mul_comm (Q.factorial:ℝ), mul_assoc, orderedSite_normalized_one, mul_one]

theorem beta_laplace_normalized {sig tau G A : ℝ} (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (a b : ℕ) (hA : 0 ≤ A)
    (hNorm : A * (∫ y in (0:ℝ)..1, y^a*(1-y)^b) = 1)
    (hG : ∀ y ∈ Icc (0:ℝ) 1,
      -tau*y + (1+1/((a+b+1:ℕ):ℝ))*Real.log (1-y+y/sig) ≤ G) :
    A * (∫ y in (0:ℝ)..1, y^a*(1-y)^b*Real.exp (-((a+b+1:ℕ):ℝ)*tau*y)) ≤
      Real.exp (((a+b+1:ℕ):ℝ)*G) * sig^(a+1) := by
  have h := mul_le_mul_of_nonneg_left (beta_laplace_le hsig0 hsig1 a b hG) hA
  calc _ ≤ A*(Real.exp (((a+b+1:ℕ):ℝ)*G)*sig^(a+1)*
             (∫ y in (0:ℝ)..1, y^a*(1-y)^b)) := h
    _ = (Real.exp (((a+b+1:ℕ):ℝ)*G)*sig^(a+1)) *
           (A*(∫ y in (0:ℝ)..1, y^a*(1-y)^b)) := by ring
    _ = _ := by rw [hNorm,mul_one]

theorem liveDuration_laplace_empty_le {Q : ℕ} {sig tau G : ℝ}
    (hG : ∀ y ∈ Icc (0:ℝ) 1, -tau*y+(1+1/(Q:ℝ))*Real.log (1-y+y/sig) ≤ G) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      Real.exp (-(Q:ℝ)*tau*liveDuration ∅ x)) ≤ Real.exp ((Q:ℝ)*G)*sig^(∅ : Finset (Fin (Q+1))).card := by
  have hG0 : 0 ≤ G := by simpa using hG 0 (by constructor <;> norm_num)
  simp only [liveDuration_empty,mul_zero,Real.exp_zero,Finset.card_empty,pow_zero,mul_one]
  rw [orderedSite_normalized_one]
  exact Real.one_le_exp (mul_nonneg (Nat.cast_nonneg Q) hG0)

theorem liveDuration_laplace_univ_le {Q : ℕ} {sig tau G : ℝ}
    (hQ : 0 < Q) (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (hG : ∀ y ∈ Icc (0:ℝ) 1, -tau*y+(1+1/(Q:ℝ))*Real.log (1-y+y/sig) ≤ G) :
    (Q.factorial:ℝ) * (∫ x in orderedSiteDomain Q,
      Real.exp (-(Q:ℝ)*tau*liveDuration univ x)) ≤ Real.exp ((Q:ℝ)*G)*sig^(univ : Finset (Fin (Q+1))).card := by
  simp only [liveDuration_univ,mul_one,card_univ,Fintype.card_fin]
  rw [orderedSite_normalized_const]
  have h := laplace_le_tilt hQ hsig0 hsig1 (by norm_num : (0:ℝ) ≤ 1) (hG 1 (by constructor <;> norm_num))
  simpa only [mul_one,sub_self,zero_add,one_div,inv_pow,div_inv_eq_mul] using h

#print axioms beta_laplace_normalized
#print axioms liveDuration_laplace_empty_le
#print axioms liveDuration_laplace_univ_le
end Spin.Structured.Placement
