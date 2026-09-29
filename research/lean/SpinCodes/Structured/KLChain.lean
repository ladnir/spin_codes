/-
The chain-rule bound behind `eq:imt-dense-exponent`.

    "Apply route domination to the moment restricted to that weight slice,
     then change its iid reference probability from αx to py.  The likelihood
     cost is exp(N D_KL(αx ‖ py)).  The chain rule bounds its exponent by
     D_KL(α‖p) + α D_KL(x‖y)."

Read as a two-stage experiment this is exact plus data processing: a
Bernoulli(α) "active" flag followed, when active, by a Bernoulli(x) bit has
joint divergence exactly `D(α‖p) + α D(x‖y)` against the `(p,y)` reference,
and the observed bit is a function of that pair.  The inequality here is the
concrete form of that argument, and the only analytic input is the two-term
log-sum inequality, itself one application of `log t ≤ t - 1`.

The bound holds on the paper's full range, `α` and `x` up to and including
`1`; the degenerate terms vanish there under the usual `0 log 0 = 0`
convention, which is what Lean's `Real.log 0 = 0` gives.
-/
import SpinCodes.Structured.Induction

namespace Spin

open Real

/-! ## The two-term log-sum inequality -/

private lemma log_ge_one_sub_inv {t : ℝ} (ht : 0 < t) : 1 - 1 / t ≤ Real.log t := by
  have h := Real.log_le_sub_one_of_pos (x := 1 / t) (by positivity)
  rw [one_div, Real.log_inv] at h
  have : (1 : ℝ) / t = t⁻¹ := one_div t
  rw [this]
  linarith

/-- **Log-sum, two terms.**  Zero numerators are allowed, with `0 log 0 = 0`. -/
theorem logSum2 {a₁ a₂ b₁ b₂ : ℝ} (ha₁ : 0 ≤ a₁) (ha₂ : 0 ≤ a₂)
    (hb₁ : 0 < b₁) (hb₂ : 0 < b₂) :
    (a₁ + a₂) * Real.log ((a₁ + a₂) / (b₁ + b₂))
      ≤ a₁ * Real.log (a₁ / b₁) + a₂ * Real.log (a₂ / b₂) := by
  have hB : (0 : ℝ) < b₁ + b₂ := by linarith
  rcases eq_or_lt_of_le (add_nonneg ha₁ ha₂) with hA | hA
  · -- every numerator is zero
    have h1 : a₁ = 0 := by linarith
    have h2 : a₂ = 0 := by linarith
    simp [h1, h2]
  · -- the generic case
    have key : ∀ a b : ℝ, 0 ≤ a → 0 < b →
        a - b * ((a₁ + a₂) / (b₁ + b₂))
          ≤ a * Real.log (a / b) - a * Real.log ((a₁ + a₂) / (b₁ + b₂)) := by
      intro a b ha hb
      rcases eq_or_lt_of_le ha with ha0 | ha0
      · rw [← ha0]
        simp only [zero_mul, zero_sub, sub_zero, neg_nonpos]
        positivity
      · have hlog : Real.log (a / b) - Real.log ((a₁ + a₂) / (b₁ + b₂))
            = Real.log (a * (b₁ + b₂) / (b * (a₁ + a₂))) := by
          rw [Real.log_div (ne_of_gt ha0) (ne_of_gt hb),
            Real.log_div (ne_of_gt hA) (ne_of_gt hB),
            Real.log_div (by positivity) (by positivity),
            Real.log_mul (ne_of_gt ha0) (ne_of_gt hB),
            Real.log_mul (ne_of_gt hb) (ne_of_gt hA)]
          ring
        have ht : (0 : ℝ) < a * (b₁ + b₂) / (b * (a₁ + a₂)) := by positivity
        have hge := log_ge_one_sub_inv ht
        have hinv : 1 / (a * (b₁ + b₂) / (b * (a₁ + a₂)))
            = b * (a₁ + a₂) / (a * (b₁ + b₂)) := by
          field_simp
        rw [hinv] at hge
        have hmul := mul_le_mul_of_nonneg_left hge (le_of_lt ha0)
        have hsimp : a * (1 - b * (a₁ + a₂) / (a * (b₁ + b₂)))
            = a - b * ((a₁ + a₂) / (b₁ + b₂)) := by field_simp
        have hright : a * Real.log (a * (b₁ + b₂) / (b * (a₁ + a₂)))
            = a * Real.log (a / b) - a * Real.log ((a₁ + a₂) / (b₁ + b₂)) := by
          rw [← hlog]; ring
        rw [hsimp, hright] at hmul
        exact hmul
    have k1 := key a₁ b₁ ha₁ hb₁
    have k2 := key a₂ b₂ ha₂ hb₂
    have hcancel : b₁ * ((a₁ + a₂) / (b₁ + b₂)) + b₂ * ((a₁ + a₂) / (b₁ + b₂))
        = a₁ + a₂ := by
      field_simp
    linarith

/-! ## Binary relative entropy -/

/-- Binary relative entropy with natural logarithms. -/
noncomputable def binKL (u v : ℝ) : ℝ :=
  u * Real.log (u / v) + (1 - u) * Real.log ((1 - u) / (1 - v))

/-- **The chain-rule bound**: `D(αx ‖ py) ≤ D(α‖p) + α D(x‖y)`. -/
theorem binKL_chain {α x p y : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1)
    (hx0 : 0 < x) (hx1 : x ≤ 1) (hp0 : 0 < p) (hp1 : p < 1)
    (hy0 : 0 < y) (hy1 : y < 1) :
    binKL (α * x) (p * y) ≤ binKL α p + α * binKL x y := by
  have hpy : (0 : ℝ) < p * y := by positivity
  have hax : (0 : ℝ) < α * x := by positivity
  -- the positive term splits exactly
  have h1 : α * x * Real.log (α * x / (p * y))
      = α * x * Real.log (α / p) + α * x * Real.log (x / y) := by
    rw [Real.log_div (ne_of_gt hax) (ne_of_gt hpy),
      Real.log_div (ne_of_gt hα0) (ne_of_gt hp0),
      Real.log_div (ne_of_gt hx0) (ne_of_gt hy0),
      Real.log_mul (ne_of_gt hα0) (ne_of_gt hx0),
      Real.log_mul (ne_of_gt hp0) (ne_of_gt hy0)]
    ring
  -- the complementary term splits too, degenerating harmlessly at `x = 1`
  have h3 : α * (1 - x) * Real.log (α * (1 - x) / (p * (1 - y)))
      = α * (1 - x) * Real.log (α / p)
        + α * (1 - x) * Real.log ((1 - x) / (1 - y)) := by
    rcases eq_or_lt_of_le hx1 with hx | hx
    · simp [← hx]
    · have hx1' : (0 : ℝ) < 1 - x := by linarith
      have hy1' : (0 : ℝ) < 1 - y := by linarith
      rw [Real.log_div (by positivity) (by positivity),
        Real.log_div (ne_of_gt hα0) (ne_of_gt hp0),
        Real.log_div (ne_of_gt hx1') (ne_of_gt hy1'),
        Real.log_mul (ne_of_gt hα0) (ne_of_gt hx1'),
        Real.log_mul (ne_of_gt hp0) (ne_of_gt hy1')]
      ring
  -- log-sum on the complementary term
  have hsum : (1 : ℝ) - α * x = α * (1 - x) + (1 - α) := by ring
  have hsumb : (1 : ℝ) - p * y = p * (1 - y) + (1 - p) := by ring
  have h2 : (1 - α * x) * Real.log ((1 - α * x) / (1 - p * y))
      ≤ α * (1 - x) * Real.log (α * (1 - x) / (p * (1 - y)))
        + (1 - α) * Real.log ((1 - α) / (1 - p)) := by
    rw [hsum, hsumb]
    exact logSum2 (by nlinarith) (by linarith) (by nlinarith) (by linarith)
  unfold binKL
  rw [h1]
  have hxy : α * (x * Real.log (x / y) + (1 - x) * Real.log ((1 - x) / (1 - y)))
      = α * x * Real.log (x / y) + α * (1 - x) * Real.log ((1 - x) / (1 - y)) := by
    ring
  rw [hxy, h3] at *
  linarith

end Spin
