/-
Empty-epoch estimates for `app:imt-fixed`.

    "During empty epochs a live state has transition `P = (I+Π)/2`, where `Π`
     replaces it by a uniform nonzero state. ... For any initial live state,
     `E[f(q_j) | q_i] = μ + 2^{-(j-i)}(f(q_i) - μ)`."

The mixing step is a fair coin between "stay" and "resample uniformly", so
the conditional mean obeys the scalar recurrence `m_{g+1} = (m_g + μ)/2`.
Everything the appendix needs from the empty epochs follows from that
recurrence and two finite geometric sums, with no probability space and no
matrix exponentials:

*  `|E S_g − μg| ≤ 2t`, because `∑ 2^{-i} ≤ 2`;
*  `Var(S_g) ≤ 3t²g`, because each row of the covariance sum is damped by
   `∑_{d≥1} 2^{-d} ≤ 1`.

Also here: `G_Q ≤ G_3`, which is monotonicity of `(1 + 1/Q)` against a
logarithm that is nonnegative on the relevant range.
-/
import SpinCodes.Structured.Collatz

namespace Spin

open Finset

/-! ## The mixing recurrence -/

/-- **The empty-epoch conditional mean.**  A fair coin between "stay" and
"resample" gives `m_g = μ + 2^{-g}(m_0 - μ)`. -/
theorem geom_recurrence {m : ℕ → ℝ} {mu : ℝ}
    (h : ∀ g, m (g + 1) = (m g + mu) / 2) (g : ℕ) :
    m g = mu + (1 / 2 : ℝ) ^ g * (m 0 - mu) := by
  induction g with
  | zero => simp
  | succ n ih =>
      rw [h n, ih]
      ring

/-! ## Two finite geometric sums -/

private lemma geom_half_le_two (g : ℕ) : ∑ i ∈ range g, (1 / 2 : ℝ) ^ i ≤ 2 := by
  rw [geom_sum_eq (by norm_num : (1 / 2 : ℝ) ≠ 1)]
  have h : (0 : ℝ) ≤ (1 / 2 : ℝ) ^ g := by positivity
  have he : ((1 / 2 : ℝ) ^ g - 1) / ((1 / 2 : ℝ) - 1) = 2 - 2 * (1 / 2 : ℝ) ^ g := by
    ring_nf
  rw [he]
  linarith

private lemma geom_half_nonneg (g : ℕ) : (0 : ℝ) ≤ ∑ i ∈ range g, (1 / 2 : ℝ) ^ i :=
  Finset.sum_nonneg fun i _ => by positivity

private lemma geom_half_succ_le_one (g : ℕ) :
    ∑ d ∈ range g, (1 / 2 : ℝ) ^ (d + 1) ≤ 1 := by
  have h : ∑ d ∈ range g, (1 / 2 : ℝ) ^ (d + 1)
      = (1 / 2 : ℝ) * ∑ d ∈ range g, (1 / 2 : ℝ) ^ d := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun d _ => by ring
  rw [h]
  linarith [geom_half_le_two g]

/-! ## The two bounds the appendix uses -/

/-- **`|E S_g − μg| ≤ 2t`.**  The deviation is the initial deviation damped by
a geometric sum bounded by `2`. -/
theorem mean_sum_bound {m : ℕ → ℝ} {mu t : ℝ}
    (h : ∀ g, m (g + 1) = (m g + mu) / 2) (ht : |m 0 - mu| ≤ t) (g : ℕ) :
    |(∑ i ∈ range g, m i) - mu * g| ≤ 2 * t := by
  have hsum : (∑ i ∈ range g, m i)
      = mu * g + (m 0 - mu) * ∑ i ∈ range g, (1 / 2 : ℝ) ^ i := by
    rw [Finset.sum_congr rfl (fun i _ => geom_recurrence h i), Finset.sum_add_distrib,
      Finset.sum_const, Finset.card_range, nsmul_eq_mul, ← Finset.sum_mul]
    ring
  rw [hsum, add_sub_cancel_left, abs_mul,
    abs_of_nonneg (geom_half_nonneg g)]
  have ht0 : 0 ≤ t := le_trans (abs_nonneg _) ht
  nlinarith [geom_half_le_two g, geom_half_nonneg g, abs_nonneg (m 0 - mu)]

/-- **`Var(S_g) ≤ 3t²g`.**  Each covariance row is damped by `∑_{d≥1} 2^{-d} ≤ 1`,
so the off-diagonal contributes at most twice the diagonal. -/
theorem var_sum_bound {g : ℕ} {t : ℝ} (v : ℕ → ℝ)
    (hv0 : ∀ i, 0 ≤ v i) (hv : ∀ i, v i ≤ t ^ 2) :
    (∑ i ∈ range g, v i)
        + 2 * ∑ i ∈ range g, ∑ d ∈ range g, (1 / 2 : ℝ) ^ (d + 1) * v i
      ≤ 3 * t ^ 2 * g := by
  have hinner : ∀ i ∈ range g,
      ∑ d ∈ range g, (1 / 2 : ℝ) ^ (d + 1) * v i ≤ v i := by
    intro i _
    rw [← Finset.sum_mul]
    calc (∑ d ∈ range g, (1 / 2 : ℝ) ^ (d + 1)) * v i
        ≤ 1 * v i := mul_le_mul_of_nonneg_right (geom_half_succ_le_one g) (hv0 i)
      _ = v i := one_mul _
  have hoff : ∑ i ∈ range g, ∑ d ∈ range g, (1 / 2 : ℝ) ^ (d + 1) * v i
      ≤ ∑ i ∈ range g, v i := Finset.sum_le_sum hinner
  have hdiag : ∑ i ∈ range g, v i ≤ t ^ 2 * g := by
    calc ∑ i ∈ range g, v i ≤ ∑ _i ∈ range g, t ^ 2 :=
          Finset.sum_le_sum fun i _ => hv i
      _ = t ^ 2 * g := by
          rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]
          ring
  linarith

/-! ## `G_Q ≤ G_3`

`G_Q = sup_y {-τy + (1 + 1/Q) log(1 - y + y/σ)}`.  On `0 ≤ y ≤ 1` with
`0 < σ < 1` the argument of the logarithm is at least `1`, so the logarithm is
nonnegative and the whole expression decreases in `Q`. -/

theorem log_arg_ge_one {sigma y : ℝ} (hsig0 : 0 < sigma) (hsig1 : sigma < 1)
    (hy0 : 0 ≤ y) : (1 : ℝ) ≤ 1 - y + y / sigma := by
  have h : y ≤ y / sigma := by
    rw [le_div_iff₀ hsig0]
    nlinarith
  linarith

/-- **`G_Q` is decreasing in `Q`**, pointwise in `y`. -/
theorem G_antitone {sigma tau y : ℝ} (hsig0 : 0 < sigma) (hsig1 : sigma < 1)
    (hy0 : 0 ≤ y) {Q Q' : ℕ} (hQ' : 0 < Q') (hQQ : Q' ≤ Q) :
    -tau * y + (1 + 1 / (Q : ℝ)) * Real.log (1 - y + y / sigma)
      ≤ -tau * y + (1 + 1 / (Q' : ℝ)) * Real.log (1 - y + y / sigma) := by
  have hlog : 0 ≤ Real.log (1 - y + y / sigma) :=
    Real.log_nonneg (log_arg_ge_one hsig0 hsig1 hy0)
  have hQ'R : (0 : ℝ) < (Q' : ℝ) := by exact_mod_cast hQ'
  have hQR : (0 : ℝ) < (Q : ℝ) := by
    have : 0 < Q := lt_of_lt_of_le hQ' hQQ
    exact_mod_cast this
  have hinv : 1 / (Q : ℝ) ≤ 1 / (Q' : ℝ) := by
    apply one_div_le_one_div_of_le hQ'R
    exact_mod_cast hQQ
  nlinarith

end Spin
