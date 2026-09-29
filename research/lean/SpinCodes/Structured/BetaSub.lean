/-
The Beta identity of `app:imt-fixed`.

    "If a fixed class path has `l` live spacings, their total duration `Y` has
     distribution `Beta(l, Q+1-l)` ... The identity
     `E[(1-Y+Y/σ)^{-(Q+1)}] = σ^l` follows by the substitution
     `u = y/[σ+(1-σ)y]` in the beta integral."

Written with `a = l-1` and `b = Q-l` as naturals — so `l = a+1`, `Q+1 = a+b+2`
and no truncated subtraction appears — the claim is

    ∫₀¹ yᵃ(1-y)ᵇ / (1-y+y/σ)^{a+b+2} dy  =  σ^{a+1} ∫₀¹ uᵃ(1-u)ᵇ du.

Dividing by the Beta function gives the paper's expectation form, and the
**Beta value is never needed**: the same `∫₀¹ uᵃ(1-u)ᵇ du` appears on both
sides and cancels.  So no Beta or Gamma machinery is used here, only the
change of variables.

Under `y = uσ/D` with `D = 1 - u(1-σ)` the substitution is exact: every power
of `D` cancels, and the transformed integrand is *literally* `σ^{a+1}uᵃ(1-u)ᵇ`
with nothing left over.
-/
import SpinCodes.Structured.WeightedNorm

namespace Spin

open Set intervalIntegral

variable {sig : ℝ}

/-! ## The substitution -/

/-- `y = uσ / (1 - u(1-σ))`, the inverse of `u = y/[σ+(1-σ)y]`. -/
noncomputable def betaSub (sig u : ℝ) : ℝ := u * sig / (1 - u * (1 - sig))

/-- Its derivative, `σ / (1 - u(1-σ))²`. -/
noncomputable def betaSub' (sig u : ℝ) : ℝ := sig / (1 - u * (1 - sig)) ^ 2

/-- The denominator is bounded below by `σ` on `[0,1]`. -/
theorem den_pos (hsig0 : 0 < sig) (hsig1 : sig ≤ 1) {u : ℝ} (hu : u ∈ Icc (0 : ℝ) 1) :
    0 < 1 - u * (1 - sig) := by
  obtain ⟨hu0, hu1⟩ := hu
  nlinarith

theorem betaSub_zero : betaSub sig 0 = 0 := by simp [betaSub]

theorem betaSub_one (hsig0 : 0 < sig) : betaSub sig 1 = 1 := by
  simp only [betaSub, one_mul]
  rw [show (1 : ℝ) - (1 - sig) = sig by ring]
  exact div_self (ne_of_gt hsig0)

theorem hasDerivAt_betaSub (hsig0 : 0 < sig) (hsig1 : sig ≤ 1) {u : ℝ}
    (hu : u ∈ Icc (0 : ℝ) 1) : HasDerivAt (betaSub sig) (betaSub' sig u) u := by
  have hD := den_pos hsig0 hsig1 hu
  have hc : HasDerivAt (fun v : ℝ => v * sig) sig u := by
    simpa using (hasDerivAt_id u).mul_const sig
  have hd : HasDerivAt (fun v : ℝ => 1 - v * (1 - sig)) (-(1 - sig)) u := by
    have h0 : HasDerivAt (fun v : ℝ => v * (1 - sig)) (1 - sig) u := by
      simpa using (hasDerivAt_id u).mul_const (1 - sig)
    simpa using h0.const_sub 1
  have h := hc.div hd (ne_of_gt hD)
  refine h.congr_deriv ?_
  unfold betaSub'
  field_simp
  ring

/-- The image of `[0,1]` lands in `[0,1]`. -/
theorem betaSub_mem (hsig0 : 0 < sig) (hsig1 : sig ≤ 1) {u : ℝ}
    (hu : u ∈ Icc (0 : ℝ) 1) : betaSub sig u ∈ Icc (0 : ℝ) 1 := by
  obtain ⟨hu0, hu1⟩ := hu
  have hD := den_pos hsig0 hsig1 ⟨hu0, hu1⟩
  constructor
  · unfold betaSub
    positivity
  · rw [betaSub, div_le_one hD]
    nlinarith

/-! ## The three algebraic identities -/

theorem one_sub_betaSub (hsig0 : 0 < sig) (hsig1 : sig ≤ 1) {u : ℝ}
    (hu : u ∈ Icc (0 : ℝ) 1) :
    1 - betaSub sig u = (1 - u) / (1 - u * (1 - sig)) := by
  have hD := den_pos hsig0 hsig1 hu
  rw [betaSub]
  field_simp
  ring

theorem tilt_betaSub (hsig0 : 0 < sig) (hsig1 : sig ≤ 1) {u : ℝ}
    (hu : u ∈ Icc (0 : ℝ) 1) :
    1 - betaSub sig u + betaSub sig u / sig = 1 / (1 - u * (1 - sig)) := by
  have hD := den_pos hsig0 hsig1 hu
  rw [betaSub]
  field_simp
  ring

/-! ## The identity -/

/-- **The Beta identity**, with the Beta value left uncancelled on both sides. -/
theorem beta_substitution (hsig0 : 0 < sig) (hsig1 : sig ≤ 1) (a b : ℕ) :
    (∫ y in (0 : ℝ)..1, y ^ a * (1 - y) ^ b / (1 - y + y / sig) ^ (a + b + 2))
      = sig ^ (a + 1) * ∫ u in (0 : ℝ)..1, u ^ a * (1 - u) ^ b := by
  set g : ℝ → ℝ :=
    fun y => y ^ a * (1 - y) ^ b / (1 - y + y / sig) ^ (a + b + 2) with hg
  -- `g` is continuous on `[0,1]`, hence on the image of the substitution
  have htilt : ∀ y ∈ Icc (0 : ℝ) 1, (1 : ℝ) ≤ 1 - y + y / sig := by
    intro y hy
    obtain ⟨hy0, hy1⟩ := hy
    have : y ≤ y / sig := by
      rw [le_div_iff₀ hsig0]
      nlinarith
    linarith
  have hgcont : ContinuousOn g (Icc (0 : ℝ) 1) := by
    refine ContinuousOn.div (by fun_prop) (by fun_prop) ?_
    intro y hy
    have h1 := htilt y hy
    positivity
  have himg : betaSub sig '' uIcc (0 : ℝ) 1 ⊆ Icc (0 : ℝ) 1 := by
    rintro y ⟨u, hu, rfl⟩
    rw [uIcc_of_le (by norm_num : (0:ℝ) ≤ 1)] at hu
    exact betaSub_mem hsig0 hsig1 hu
  -- change of variables
  have hcv := integral_comp_mul_deriv' (a := (0 : ℝ)) (b := 1)
    (f := betaSub sig) (f' := betaSub' sig) (g := g)
    (fun u hu => hasDerivAt_betaSub hsig0 hsig1 (by
      rwa [uIcc_of_le (by norm_num : (0:ℝ) ≤ 1)] at hu))
    (by
      intro u hu
      rw [uIcc_of_le (by norm_num : (0:ℝ) ≤ 1)] at hu
      have hD := den_pos hsig0 hsig1 hu
      refine ContinuousWithinAt.div continuousWithinAt_const ?_ (by positivity)
      · exact (continuousWithinAt_const.sub
          (continuousWithinAt_id.mul continuousWithinAt_const)).pow 2)
    (hgcont.mono himg)
  rw [betaSub_zero, betaSub_one hsig0] at hcv
  rw [← hcv, ← intervalIntegral.integral_const_mul]
  refine intervalIntegral.integral_congr ?_
  intro u hu
  rw [uIcc_of_le (by norm_num : (0:ℝ) ≤ 1)] at hu
  have hD := den_pos hsig0 hsig1 hu
  have h1 : 1 - u * sig / (1 - u * (1 - sig)) = (1 - u) / (1 - u * (1 - sig)) := by
    field_simp; ring
  have h2 : 1 - u * sig / (1 - u * (1 - sig)) + u * sig / (1 - u * (1 - sig)) / sig
      = 1 / (1 - u * (1 - sig)) := by
    field_simp; ring
  have hDne : (1 - u * (1 - sig)) ≠ 0 := ne_of_gt hD
  have e3' : ((1 - u) / (1 - u * (1 - sig))
      + u * sig / (1 - u * (1 - sig)) / sig) ^ (a + b + 2)
      = ((1 - u * (1 - sig)) ^ (a + b + 2))⁻¹ := by
    have hone : (1 - u) / (1 - u * (1 - sig))
        + u * sig / (1 - u * (1 - sig)) / sig = 1 / (1 - u * (1 - sig)) := by
      field_simp
      ring
    rw [hone, one_div, inv_pow]
  simp only [Function.comp_apply, hg, betaSub, betaSub', h1, e3']
  rw [div_eq_mul_inv ((u * sig / (1 - u * (1 - sig))) ^ a
      * ((1 - u) / (1 - u * (1 - sig))) ^ b)
    (((1 - u * (1 - sig)) ^ (a + b + 2))⁻¹), inv_inv, div_pow, div_pow]
  field_simp
  ring

end Spin
