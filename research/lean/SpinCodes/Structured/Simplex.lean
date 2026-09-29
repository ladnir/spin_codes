/-
The continuum region kernel `K_a(θ)` of `app:imt-fixed`.

    K_a(θ) := a! ∫_{u ≥ 0, ∑ u_i = 1} D_γ(u_0) P_0 D_γ(u_1) ⋯ P_0 D_γ(u_a) du

with `D_γ(u) = diag(1, e^{-γu})` and, in the enlarged form, `P_+ = [[0,1],[r,1]]`.

Two separable halves:

*  the **integrand**, an exact matrix identity, checked against Mathlib's own
   matrix product rather than against a transcription of it — this is where a
   sign or a transposition in the paper's `K_a^+` would show up;
*  the **scalar integrals** over `[0,1]`, each proved by exhibiting an
   antiderivative, so neither substitution nor integration by parts is needed.

For `a = 1` the two combine directly into the paper's

    K_1^+ = [[0, f], [r f, a_0]],   f = (1 - e^{-γ})/γ,  a_0 = e^{-γ}.

The `a = 2` integrand is also verified here; assembling it needs a genuine
2-simplex integral and is separate.
-/
import SpinCodes.Structured.Coarse

namespace Spin

open Real intervalIntegral

/-! ## The two factors -/

/-- `D_γ(u) = diag(1, e^{-γu})`. -/
noncomputable def Dg (gam u : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![1, 0; 0, Real.exp (-gam * u)]

/-- The enlarged impulse matrix `P_+ = [[0,1],[r,1]]`. -/
def Pplus (r : ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  !![0, 1; r, 1]

/-! ## The integrands

Checked against Mathlib's matrix multiplication. -/

/-- The `a = 1` integrand, exactly as the matrix product yields it. -/
theorem integrand_one (gam r u : ℝ) :
    Dg gam u * Pplus r * Dg gam (1 - u)
      = !![0, Real.exp (-gam * (1 - u));
           r * Real.exp (-gam * u),
           Real.exp (-gam * u) * Real.exp (-gam * (1 - u))] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Dg, Pplus, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-- The two empty durations of an `a = 1` region total `1`. -/
theorem exp_prod_one (gam u : ℝ) :
    Real.exp (-gam * u) * Real.exp (-gam * (1 - u)) = Real.exp (-gam) := by
  rw [← Real.exp_add, Real.exp_eq_exp]
  ring

/-- The `a = 2` integrand, exactly as the matrix product yields it. -/
theorem integrand_two (gam r u0 u1 u2 : ℝ) :
    Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam u2
      = !![r * Real.exp (-gam * u1), Real.exp (-gam * u1) * Real.exp (-gam * u2);
           r * (Real.exp (-gam * u0) * Real.exp (-gam * u1)),
           r * (Real.exp (-gam * u0) * Real.exp (-gam * u2))
             + Real.exp (-gam * u0) * Real.exp (-gam * u1) * Real.exp (-gam * u2)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Dg, Pplus, Matrix.mul_apply, Fin.sum_univ_two] <;> ring

/-- On the simplex the three durations total `1`, so the `(2,2)` entry picks up
the constant `e^{-γ}`. -/
theorem exp_prod_two {gam u0 u1 u2 : ℝ} (h : u0 + u1 + u2 = 1) :
    Real.exp (-gam * u0) * Real.exp (-gam * u1) * Real.exp (-gam * u2)
      = Real.exp (-gam) := by
  rw [← Real.exp_add, ← Real.exp_add, Real.exp_eq_exp]
  linear_combination (-gam) * h

/-! ## The scalar integrals

Each by exhibiting an antiderivative, so no substitution and no parts. -/

private lemma hasDerivAt_expAffine (c d u : ℝ) :
    HasDerivAt (fun v : ℝ => Real.exp (c * v + d)) (Real.exp (c * u + d) * c) u := by
  have h1 : HasDerivAt (fun v : ℝ => c * v + d) c u := by
    simpa using ((hasDerivAt_id u).const_mul c).add_const d
  simpa using h1.exp

private lemma hasDerivAt_expLin (c u : ℝ) :
    HasDerivAt (fun v : ℝ => Real.exp (c * v)) (Real.exp (c * u) * c) u := by
  have h1 : HasDerivAt (fun v : ℝ => c * v) c u := by
    simpa using (hasDerivAt_id u).const_mul c
  simpa using h1.exp

/-- `∫₀¹ e^{-γu} du = (1 - e^{-γ})/γ`. -/
theorem integral_expNeg {gam : ℝ} (hg : gam ≠ 0) :
    (∫ u in (0 : ℝ)..1, Real.exp (-gam * u)) = (1 - Real.exp (-gam)) / gam := by
  have hF : ∀ u ∈ Set.uIcc (0 : ℝ) 1,
      HasDerivAt (fun v : ℝ => -(Real.exp (-gam * v) / gam)) (Real.exp (-gam * u)) u := by
    intro u _
    have h := ((hasDerivAt_expLin (-gam) u).div_const gam).neg
    convert h using 1
    field_simp
  rw [integral_eq_sub_of_hasDerivAt hF
    (Continuous.intervalIntegrable (by fun_prop) 0 1)]
  simp only [mul_one, mul_zero, Real.exp_zero]
  ring

private lemma exp_rev (gam u : ℝ) :
    Real.exp (-gam * (1 - u)) = Real.exp (gam * u + -gam) := by
  congr 1
  ring

/-- `∫₀¹ e^{-γ(1-u)} du = (1 - e^{-γ})/γ`. -/
theorem integral_expNeg_rev {gam : ℝ} (hg : gam ≠ 0) :
    (∫ u in (0 : ℝ)..1, Real.exp (-gam * (1 - u))) = (1 - Real.exp (-gam)) / gam := by
  simp only [exp_rev]
  have hF : ∀ u ∈ Set.uIcc (0 : ℝ) 1,
      HasDerivAt (fun v : ℝ => Real.exp (gam * v + -gam) / gam)
        (Real.exp (gam * u + -gam)) u := by
    intro u _
    have h := (hasDerivAt_expAffine gam (-gam) u).div_const gam
    convert h using 1
    field_simp
  rw [integral_eq_sub_of_hasDerivAt hF
    (Continuous.intervalIntegrable (by fun_prop) 0 1)]
  simp only [mul_one, mul_zero, zero_add, Real.exp_zero, add_neg_cancel]
  ring

/-- `∫₀¹ u e^{-γu} du = (1 - (1+γ)e^{-γ})/γ²`. -/
theorem integral_mul_expNeg {gam : ℝ} (hg : gam ≠ 0) :
    (∫ u in (0 : ℝ)..1, u * Real.exp (-gam * u))
      = (1 - (1 + gam) * Real.exp (-gam)) / gam ^ 2 := by
  have hF : ∀ u ∈ Set.uIcc (0 : ℝ) 1,
      HasDerivAt (fun v : ℝ => -((v / gam + 1 / gam ^ 2) * Real.exp (-gam * v)))
        (u * Real.exp (-gam * u)) u := by
    intro u _
    have hlin : HasDerivAt (fun v : ℝ => v / gam + 1 / gam ^ 2) (1 / gam) u := by
      simpa using ((hasDerivAt_id u).div_const gam).add_const (1 / gam ^ 2)
    have h := (hlin.mul (hasDerivAt_expLin (-gam) u)).neg
    convert h using 1
    field_simp
    ring
  rw [integral_eq_sub_of_hasDerivAt hF
    (Continuous.intervalIntegrable (by fun_prop) 0 1)]
  simp only [mul_one, mul_zero, Real.exp_zero, zero_div, zero_add]
  field_simp
  ring

/-! ## `K_1^+`

The `a = 1` kernel, entry by entry. -/

theorem K1_entry_00 (gam r : ℝ) :
    (∫ u in (0 : ℝ)..1, (Dg gam u * Pplus r * Dg gam (1 - u)) 0 0) = 0 := by
  simp [integrand_one]

theorem K1_entry_01 {gam : ℝ} (r : ℝ) (hg : gam ≠ 0) :
    (∫ u in (0 : ℝ)..1, (Dg gam u * Pplus r * Dg gam (1 - u)) 0 1)
      = (1 - Real.exp (-gam)) / gam := by
  have h : ∀ u : ℝ, (Dg gam u * Pplus r * Dg gam (1 - u)) 0 1
      = Real.exp (-gam * (1 - u)) := by
    intro u; rw [integrand_one]; rfl
  simp only [h]
  exact integral_expNeg_rev hg

theorem K1_entry_10 {gam : ℝ} (r : ℝ) (hg : gam ≠ 0) :
    (∫ u in (0 : ℝ)..1, (Dg gam u * Pplus r * Dg gam (1 - u)) 1 0)
      = r * ((1 - Real.exp (-gam)) / gam) := by
  have h : ∀ u : ℝ, (Dg gam u * Pplus r * Dg gam (1 - u)) 1 0
      = r * Real.exp (-gam * u) := by
    intro u; rw [integrand_one]; rfl
  simp only [h]
  rw [intervalIntegral.integral_const_mul, integral_expNeg hg]

theorem K1_entry_11 (gam r : ℝ) :
    (∫ u in (0 : ℝ)..1, (Dg gam u * Pplus r * Dg gam (1 - u)) 1 1)
      = Real.exp (-gam) := by
  have h : ∀ u : ℝ, (Dg gam u * Pplus r * Dg gam (1 - u)) 1 1 = Real.exp (-gam) := by
    intro u; rw [integrand_one]; exact exp_prod_one gam u
  simp only [h, intervalIntegral.integral_const]
  norm_num
/-! ## `K_2^+`

The `a = 2` kernel is an iterated integral over the triangle
`{u₀, u₁ ≥ 0, u₀ + u₁ ≤ 1}` with `u₂ = 1 - u₀ - u₁`, written as a nested
`intervalIntegral` with a variable inner endpoint.  No Fubini and no measure
theory: the inner integral is evaluated by an antiderivative for each fixed
`u₀`, and the outer one again. -/

/-- Variable-endpoint exponential integral. -/
theorem integral_expLin {c : ℝ} (hc : c ≠ 0) (b : ℝ) :
    (∫ u in (0 : ℝ)..b, Real.exp (c * u)) = (Real.exp (c * b) - 1) / c := by
  have hF : ∀ u ∈ Set.uIcc (0 : ℝ) b,
      HasDerivAt (fun v : ℝ => Real.exp (c * v) / c) (Real.exp (c * u)) u := by
    intro u _
    have h := (hasDerivAt_expLin c u).div_const c
    convert h using 1
    field_simp
  rw [integral_eq_sub_of_hasDerivAt hF
    (Continuous.intervalIntegrable (by fun_prop) 0 b)]
  simp only [mul_zero, Real.exp_zero]
  ring

/-- The reversed linear-times-exponential integral. -/
theorem integral_one_sub_mul_expNeg {gam : ℝ} (hg : gam ≠ 0) :
    (∫ u in (0 : ℝ)..1, (1 - u) * Real.exp (-gam * (1 - u)))
      = (1 - (1 + gam) * Real.exp (-gam)) / gam ^ 2 := by
  have hF : ∀ u ∈ Set.uIcc (0 : ℝ) 1,
      HasDerivAt (fun v : ℝ => ((1 - v) / gam + 1 / gam ^ 2) * Real.exp (-gam * (1 - v)))
        ((1 - u) * Real.exp (-gam * (1 - u))) u := by
    intro u _
    have hlin : HasDerivAt (fun v : ℝ => (1 - v) / gam + 1 / gam ^ 2) (-(1 / gam)) u := by
      have h0 : HasDerivAt (fun v : ℝ => 1 - v) (-1 : ℝ) u := by
        simpa using (hasDerivAt_id u).const_sub 1
      have h1 : HasDerivAt (fun v : ℝ => (1 - v) / gam) (-(1 / gam)) u := by
        have := h0.div_const gam
        convert this using 1
        field_simp
      simpa using h1.add_const (1 / gam ^ 2)
    have hexp : HasDerivAt (fun v : ℝ => Real.exp (-gam * (1 - v)))
        (Real.exp (-gam * (1 - u)) * gam) u := by
      have h2 : HasDerivAt (fun v : ℝ => -gam * (1 - v)) gam u := by
        have hh : HasDerivAt (fun v : ℝ => 1 - v) (-1 : ℝ) u := by
          simpa using (hasDerivAt_id u).const_sub 1
        have := hh.const_mul (-gam)
        convert this using 1
        ring
      simpa using h2.exp
    have h := hlin.mul hexp
    convert h using 1
    field_simp
    ring
  rw [integral_eq_sub_of_hasDerivAt hF
    (Continuous.intervalIntegrable (by fun_prop) 0 1)]
  simp only [sub_self, sub_zero, mul_zero, mul_one, Real.exp_zero, zero_div, zero_add]
  field_simp
  ring

/-! ### The four entries -/

theorem K2_entry_00 {gam : ℝ} (r : ℝ) (hg : gam ≠ 0) :
    (2 : ℝ) * ∫ u0 in (0 : ℝ)..1, ∫ u1 in (0 : ℝ)..(1 - u0),
        (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 0 0
      = r * (2 * (gam - 1 + Real.exp (-gam)) / gam ^ 2) := by
  have hentry : ∀ u0 u1 : ℝ,
      (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 0 0
        = r * Real.exp (-gam * u1) := by
    intro u0 u1; rw [integrand_two]; rfl
  have hinner : ∀ u0 : ℝ,
      (∫ u1 in (0 : ℝ)..(1 - u0),
          (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 0 0)
        = r * ((Real.exp (-gam * (1 - u0)) - 1) / (-gam)) := by
    intro u0
    simp only [hentry]
    rw [intervalIntegral.integral_const_mul, integral_expLin (by simpa using hg)]
  simp only [hinner]
  rw [intervalIntegral.integral_const_mul]
  have hsplit : (∫ u0 in (0 : ℝ)..1, (Real.exp (-gam * (1 - u0)) - 1) / (-gam))
      = ((1 - Real.exp (-gam)) / gam - 1) / (-gam) := by
    rw [intervalIntegral.integral_div, intervalIntegral.integral_sub
      (Continuous.intervalIntegrable (by fun_prop) 0 1)
      (Continuous.intervalIntegrable (by fun_prop) 0 1),
      integral_expNeg_rev hg]
    simp
  rw [hsplit]
  field_simp
  ring

theorem K2_entry_01 {gam : ℝ} (r : ℝ) (hg : gam ≠ 0) :
    (2 : ℝ) * ∫ u0 in (0 : ℝ)..1, ∫ u1 in (0 : ℝ)..(1 - u0),
        (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 0 1
      = 2 * (1 - (1 + gam) * Real.exp (-gam)) / gam ^ 2 := by
  have hentry : ∀ u0 u1 : ℝ,
      (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 0 1
        = Real.exp (-gam * (1 - u0)) := by
    intro u0 u1
    rw [integrand_two]
    show Real.exp (-gam * u1) * Real.exp (-gam * (1 - u0 - u1))
      = Real.exp (-gam * (1 - u0))
    rw [← Real.exp_add, Real.exp_eq_exp]
    ring
  have hinner : ∀ u0 : ℝ,
      (∫ u1 in (0 : ℝ)..(1 - u0),
          (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 0 1)
        = (1 - u0) * Real.exp (-gam * (1 - u0)) := by
    intro u0
    simp only [hentry, intervalIntegral.integral_const, smul_eq_mul, sub_zero]
  simp only [hinner]
  rw [integral_one_sub_mul_expNeg hg]
  ring
private lemma hasDerivAt_expRev (gam u : ℝ) :
    HasDerivAt (fun v : ℝ => Real.exp (-gam * (1 - v)))
      (Real.exp (-gam * (1 - u)) * gam) u := by
  have hh : HasDerivAt (fun v : ℝ => 1 - v) (-1 : ℝ) u := by
    simpa using (hasDerivAt_id u).const_sub 1
  have h2 : HasDerivAt (fun v : ℝ => -gam * (1 - v)) gam u := by
    have := hh.const_mul (-gam)
    convert this using 1
    ring
  simpa using h2.exp

/-- Variable-endpoint reversed exponential integral. -/
theorem integral_expRev {gam : ℝ} (hg : gam ≠ 0) (b : ℝ) :
    (∫ u in (0 : ℝ)..b, Real.exp (-gam * (1 - u)))
      = (Real.exp (-gam * (1 - b)) - Real.exp (-gam)) / gam := by
  have hF : ∀ u ∈ Set.uIcc (0 : ℝ) b,
      HasDerivAt (fun v : ℝ => Real.exp (-gam * (1 - v)) / gam)
        (Real.exp (-gam * (1 - u))) u := by
    intro u _
    have h := (hasDerivAt_expRev gam u).div_const gam
    convert h using 1
    field_simp
  rw [integral_eq_sub_of_hasDerivAt hF
    (Continuous.intervalIntegrable (by fun_prop) 0 b)]
  simp only [sub_zero, mul_one]
  ring

theorem K2_entry_10 {gam : ℝ} (r : ℝ) (hg : gam ≠ 0) :
    (2 : ℝ) * ∫ u0 in (0 : ℝ)..1, ∫ u1 in (0 : ℝ)..(1 - u0),
        (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 1 0
      = r * (2 * (1 - (1 + gam) * Real.exp (-gam)) / gam ^ 2) := by
  have hentry : ∀ u0 u1 : ℝ,
      (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 1 0
        = r * Real.exp (-gam * u0) * Real.exp (-gam * u1) := by
    intro u0 u1
    rw [integrand_two]
    show r * (Real.exp (-gam * u0) * Real.exp (-gam * u1))
      = r * Real.exp (-gam * u0) * Real.exp (-gam * u1)
    ring
  have hinner : ∀ u0 : ℝ,
      (∫ u1 in (0 : ℝ)..(1 - u0),
          (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 1 0)
        = r * (Real.exp (-gam * u0) - Real.exp (-gam)) / gam := by
    intro u0
    simp only [hentry]
    rw [intervalIntegral.integral_const_mul, integral_expLin (by simpa using hg)]
    have hAB := exp_prod_one gam u0
    simp only [neg_mul] at hAB ⊢
    field_simp
    linear_combination (-r) * hAB
  simp only [hinner]
  rw [intervalIntegral.integral_div, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_sub
      (Continuous.intervalIntegrable (by fun_prop) 0 1)
      (Continuous.intervalIntegrable (by fun_prop) 0 1),
    integral_expNeg hg]
  simp only [intervalIntegral.integral_const, smul_eq_mul, sub_zero, one_mul]
  field_simp
  ring

theorem K2_entry_11 {gam : ℝ} (r : ℝ) (hg : gam ≠ 0) :
    (2 : ℝ) * ∫ u0 in (0 : ℝ)..1, ∫ u1 in (0 : ℝ)..(1 - u0),
        (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 1 1
      = r * (2 * (1 - (1 + gam) * Real.exp (-gam)) / gam ^ 2) + Real.exp (-gam) := by
  have hentry : ∀ u0 u1 : ℝ,
      (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 1 1
        = r * Real.exp (-gam * (1 - u1)) + Real.exp (-gam) := by
    intro u0 u1
    rw [integrand_two]
    show r * (Real.exp (-gam * u0) * Real.exp (-gam * (1 - u0 - u1)))
        + Real.exp (-gam * u0) * Real.exp (-gam * u1) * Real.exp (-gam * (1 - u0 - u1))
      = r * Real.exp (-gam * (1 - u1)) + Real.exp (-gam)
    have h1 : Real.exp (-gam * u0) * Real.exp (-gam * (1 - u0 - u1))
        = Real.exp (-gam * (1 - u1)) := by
      rw [← Real.exp_add, Real.exp_eq_exp]; ring
    have h2 : Real.exp (-gam * u0) * Real.exp (-gam * u1) * Real.exp (-gam * (1 - u0 - u1))
        = Real.exp (-gam) := by
      rw [← Real.exp_add, ← Real.exp_add, Real.exp_eq_exp]; ring
    rw [h1, h2]
  have hinner : ∀ u0 : ℝ,
      (∫ u1 in (0 : ℝ)..(1 - u0),
          (Dg gam u0 * Pplus r * Dg gam u1 * Pplus r * Dg gam (1 - u0 - u1)) 1 1)
        = r * ((Real.exp (-gam * u0) - Real.exp (-gam)) / gam)
          + Real.exp (-gam) * (1 - u0) := by
    intro u0
    simp only [hentry]
    rw [intervalIntegral.integral_add
        (Continuous.intervalIntegrable (by fun_prop) 0 (1 - u0))
        (Continuous.intervalIntegrable (by fun_prop) 0 (1 - u0)),
      intervalIntegral.integral_const_mul, integral_expRev hg,
      intervalIntegral.integral_const]
    simp only [smul_eq_mul, sub_zero, sub_sub_cancel]
    ring
  simp only [hinner]
  rw [intervalIntegral.integral_add
      (Continuous.intervalIntegrable (by fun_prop) 0 1)
      (Continuous.intervalIntegrable (by fun_prop) 0 1),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_div,
    intervalIntegral.integral_sub
      (Continuous.intervalIntegrable (by fun_prop) 0 1)
      (Continuous.intervalIntegrable (by fun_prop) 0 1),
    integral_expNeg hg]
  simp only [intervalIntegral.integral_const, smul_eq_mul, sub_zero, one_mul]
  rw [intervalIntegral.integral_sub
      (Continuous.intervalIntegrable (by fun_prop) 0 1)
      (Continuous.intervalIntegrable (by fun_prop) 0 1),
    intervalIntegral.integral_const, integral_id]
  simp only [smul_eq_mul, sub_zero, one_mul]
  field_simp
  ring

end Spin
