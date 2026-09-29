/-
Rigorous rational enclosures of `Real.log`.

See `DECISIONS.md` D1 for why this is built here rather than imported.  The
hard part is already in Mathlib:

  `Real.abs_log_sub_add_sum_range_le : |x| < 1 →
     |(∑ i ∈ range n, x^(i+1)/(i+1)) + log (1 - x)| ≤ |x|^(n+1) / (1 - |x|)`

Everything below is the two-sided reading of that bound, plus the argument
reduction that makes it converge fast.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Numeric

open Finset Real

/-- Partial sum of the `log (1-x)` series, negated: the approximation to
`log (1-x)` is `-logSeries n x`. -/
noncomputable def logSeries (n : ℕ) (x : ℝ) : ℝ := ∑ i ∈ range n, x ^ (i + 1) / (i + 1)

/-- The proved remainder bound. -/
noncomputable def logRem (n : ℕ) (x : ℝ) : ℝ := |x| ^ (n + 1) / (1 - |x|)

/-- **Two-sided enclosure of `log (1-x)`.** -/
theorem log_one_sub_mem_Icc {x : ℝ} (h : |x| < 1) (n : ℕ) :
    -logSeries n x - logRem n x ≤ Real.log (1 - x) ∧
      Real.log (1 - x) ≤ -logSeries n x + logRem n x := by
  have hb := Real.abs_log_sub_add_sum_range_le h n
  rw [abs_le] at hb
  unfold logSeries logRem
  constructor <;> linarith [hb.1, hb.2]

/-- The enclosure in the form used for a rational argument `q = 1 - x`. -/
theorem log_mem_Icc {q : ℝ} (hq0 : 0 < q) (hq2 : q < 2) (n : ℕ) :
    -logSeries n (1 - q) - logRem n (1 - q) ≤ Real.log q ∧
      Real.log q ≤ -logSeries n (1 - q) + logRem n (1 - q) := by
  have h : |1 - q| < 1 := by
    rw [abs_lt]
    constructor <;> linarith
  have := log_one_sub_mem_Icc h n
  rwa [show (1 : ℝ) - (1 - q) = q by ring] at this

/-- The remainder shrinks geometrically: at `|x| ≤ 1/3` — which argument
reduction to `[2/3, 4/3]` always achieves — it is below `(3/2) · 3^{-(n+1)}`. -/
theorem logRem_le {x : ℝ} (hx : |x| ≤ 1 / 3) (n : ℕ) :
    logRem n x ≤ (3 / 2) * (1 / 3) ^ (n + 1) := by
  unfold logRem
  have h1 : |x| ^ (n + 1) ≤ (1 / 3 : ℝ) ^ (n + 1) :=
    pow_le_pow_left₀ (abs_nonneg x) hx _
  have h2 : (2 / 3 : ℝ) ≤ 1 - |x| := by linarith
  have h3 : (0 : ℝ) < 1 - |x| := by linarith
  rw [div_le_iff₀ h3]
  nlinarith [pow_nonneg (by norm_num : (0:ℝ) ≤ 1/3) (n + 1), abs_nonneg x]

/-! ## Worked bound

Evidence that the enclosure reaches the precision the certificates need.  The
majorant's largest residual is `-1.70e-9`, so roughly ten digits are required;
at `|x| = 1/10` twelve terms already give a remainder below `1e-12`. -/

example : logRem 11 (1 / 10 : ℝ) ≤ 1 / 10 ^ 11 := by
  unfold logRem
  rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/10)]
  norm_num

/-- `log (9/10)` is pinned to within `1e-11` by the degree-11 truncation. -/
example :
    -logSeries 11 (1 / 10 : ℝ) - 1 / 10 ^ 11 ≤ Real.log (9 / 10) ∧
      Real.log (9 / 10) ≤ -logSeries 11 (1 / 10 : ℝ) + 1 / 10 ^ 11 := by
  have hb := log_mem_Icc (q := (9 / 10 : ℝ)) (by norm_num) (by norm_num) 11
  have hx : (1 : ℝ) - 9 / 10 = 1 / 10 := by norm_num
  rw [hx] at hb
  have hr : logRem 11 (1 / 10 : ℝ) ≤ 1 / 10 ^ 11 := by
    unfold logRem
    rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1/10)]
    norm_num
  constructor <;> linarith [hb.1, hb.2]

end Spin.Numeric
