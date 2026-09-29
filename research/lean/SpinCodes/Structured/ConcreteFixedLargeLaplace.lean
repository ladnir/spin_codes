import SpinCodes.Structured.BetaSub
import SpinCodes.Structured.EmptyEpoch

/-! Scalar Laplace domination used for the actual ordered-spacing integral.
No distributional identification is assumed here: the beta integral is explicit. -/
noncomputable section
namespace Spin.Structured.Placement
open Set MeasureTheory intervalIntegral

/-- A pointwise bound on the tilted exponent gives the rational-power envelope. -/
theorem laplace_le_tilt {Q : ℕ} {sig tau G y : ℝ}
    (hQ : 0 < Q) (hsig0 : 0 < sig) (hsig1 : sig < 1) (hy : 0 ≤ y)
    (hG : -tau*y + (1 + 1/(Q:ℝ))*Real.log (1-y+y/sig) ≤ G) :
    Real.exp (-(Q:ℝ)*tau*y) ≤
      Real.exp ((Q:ℝ)*G) / (1-y+y/sig)^(Q+1) := by
  have hQR : (0:ℝ) < Q := by exact_mod_cast hQ
  have ht : 0 < 1-y+y/sig := lt_of_lt_of_le (by norm_num) (Spin.log_arg_ge_one hsig0 hsig1 hy)
  have he := mul_le_mul_of_nonneg_left hG hQR.le
  have hcast : ((Q+1:ℕ):ℝ) = (Q:ℝ)+1 := by norm_cast
  have he' : -(Q:ℝ)*tau*y ≤ (Q:ℝ)*G - ((Q+1:ℕ):ℝ)*Real.log (1-y+y/sig) := by
    have hc : (Q:ℝ)*(1+1/(Q:ℝ)) = (Q:ℝ)+1 := by field_simp
    rw [hcast]
    rw [mul_add, ← mul_assoc (Q:ℝ) (1+1/(Q:ℝ)), hc] at he
    nlinarith [he]
  have hp : (1-y+y/sig)^(Q+1) = Real.exp (((Q+1:ℕ):ℝ)*Real.log (1-y+y/sig)) := by
    rw [Real.exp_nat_mul, Real.exp_log ht]
  rw [hp, ← Real.exp_sub]
  exact Real.exp_le_exp.mpr he'

/-- All positive beta shapes with total spacing count `a+b+2`. -/
theorem beta_laplace_le {sig tau G : ℝ} (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (a b : ℕ)
    (hG : ∀ y ∈ Icc (0:ℝ) 1,
      -tau*y + (1+1/((a+b+1:ℕ):ℝ))*Real.log (1-y+y/sig) ≤ G) :
    (∫ y in (0:ℝ)..1, y^a*(1-y)^b*Real.exp (-((a+b+1:ℕ):ℝ)*tau*y)) ≤
      Real.exp (((a+b+1:ℕ):ℝ)*G) * sig^(a+1) *
        ∫ y in (0:ℝ)..1, y^a*(1-y)^b := by
  let E := Real.exp (((a+b+1:ℕ):ℝ)*G)
  have hleft : Continuous (fun y : ℝ => y^a*(1-y)^b*Real.exp (-((a+b+1:ℕ):ℝ)*tau*y)) := by fun_prop
  have hright : ContinuousOn (fun y : ℝ => E*(y^a*(1-y)^b/(1-y+y/sig)^(a+b+2))) (Icc 0 1) := by
    apply ContinuousOn.mul continuousOn_const
    apply ContinuousOn.div (by fun_prop) (by fun_prop)
    intro y hy
    exact pow_ne_zero _ (ne_of_gt (lt_of_lt_of_le (by norm_num) (Spin.log_arg_ge_one hsig0 hsig1 hy.1)))
  have hmono := intervalIntegral.integral_mono_on (μ := volume) (by norm_num : (0:ℝ) ≤ 1)
    (hleft.intervalIntegrable 0 1) (hright.intervalIntegrable_of_Icc (by norm_num))
    (fun y hy => ?_)
  · rw [intervalIntegral.integral_const_mul, Spin.beta_substitution hsig0 hsig1.le] at hmono
    simpa only [E, mul_assoc] using hmono
  · have hd : 0 ≤ y^a*(1-y)^b := mul_nonneg (pow_nonneg hy.1 _) (pow_nonneg (sub_nonneg.mpr hy.2) _)
    have h := mul_le_mul_of_nonneg_left (laplace_le_tilt (by omega : 0 < a+b+1) hsig0 hsig1 hy.1 (hG y hy)) hd
    convert h using 1 <;> ring

/-- A single three-site exponent envelope controls every larger occupation. -/
theorem exponent_le_three {Q : ℕ} {sig tau G : ℝ} (hQ : 3 ≤ Q)
    (hsig0 : 0 < sig) (hsig1 : sig < 1)
    (hG : ∀ y ∈ Icc (0:ℝ) 1, -tau*y + (1+1/(3:ℝ))*Real.log (1-y+y/sig) ≤ G)
    (y : ℝ) (hy : y ∈ Icc (0:ℝ) 1) :
    -tau*y + (1+1/(Q:ℝ))*Real.log (1-y+y/sig) ≤ G :=
  (Spin.G_antitone hsig0 hsig1 hy.1 (by decide : 0 < 3) hQ).trans (hG y hy)

#print axioms laplace_le_tilt
#print axioms beta_laplace_le
#print axioms exponent_le_three
end Spin.Structured.Placement
