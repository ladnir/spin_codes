/-
Computable rational enclosures of `Real.log`.

`LogBounds` (T5) gives the *mathematical* enclosure, stated over `ℝ`.  That is
enough for a handful of constants (T6h, T6j used it directly), but the box
covers of `eq:structured-ba-dense-tail` and `eq:ba-spectrum-majorant` need
6,747 and 10,721 checks.  Driving those through `norm_num` one at a time would
cost hours of elaboration.

So this file adds the *computable* layer: `ℚ`-valued bounds that the kernel can
evaluate, with soundness proved against `Real.log` once.  A whole cover then
becomes a single `decide` over a list, the way the sparse Collatz certificate
was discharged in T3.

Argument reduction is by powers of two, with the shift supplied as a parameter
rather than searched for: the certificate generator knows it, and passing it in
keeps these definitions total and cheap.  All arguments here are probabilities
in `(0,1)`, so the shift is always upward.
-/
import SpinCodes.Numeric.LogBounds

set_option linter.unusedSectionVars false

namespace Spin.Numeric

open Finset

/-- A rational strictly below `Real.log 2`. -/
def log2Lo : ℚ := 6931471803 / 10000000000

/-- A rational strictly above `Real.log 2`. -/
def log2Hi : ℚ := 6931471808 / 10000000000

lemma log2Lo_lt : (log2Lo : ℝ) < Real.log 2 := by
  have h := Real.log_two_gt_d9
  unfold log2Lo
  push_cast
  norm_num at h ⊢
  linarith

lemma lt_log2Hi : Real.log 2 < (log2Hi : ℝ) := by
  have h := Real.log_two_lt_d9
  unfold log2Hi
  push_cast
  norm_num at h ⊢
  linarith

/-! ## The series, over `ℚ` -/

/-- Computable partial sum of the `log(1-x)` series. -/
def ratLogSeries (n : ℕ) (x : ℚ) : ℚ := ∑ i ∈ range n, x ^ (i + 1) / (i + 1)

/-- Computable remainder bound. -/
def ratLogRem (n : ℕ) (x : ℚ) : ℚ := |x| ^ (n + 1) / (1 - |x|)

lemma ratLogSeries_cast (n : ℕ) (x : ℚ) :
    ((ratLogSeries n x : ℚ) : ℝ) = logSeries n (x : ℝ) := by
  unfold ratLogSeries logSeries
  push_cast
  rfl

lemma ratLogRem_cast (n : ℕ) (x : ℚ) :
    ((ratLogRem n x : ℚ) : ℝ) = logRem n (x : ℝ) := by
  unfold ratLogRem logRem
  push_cast
  rfl

/-! ## Enclosures with a power-of-two shift -/

/-- Upper bound for `Real.log q`, reducing by `2^k`. -/
def logUpper (q : ℚ) (k n : ℕ) : ℚ :=
  (-ratLogSeries n (1 - q * 2 ^ k) + ratLogRem n (1 - q * 2 ^ k)) - k * log2Lo

/-- Lower bound for `Real.log q`, reducing by `2^k`. -/
def logLower (q : ℚ) (k n : ℕ) : ℚ :=
  (-ratLogSeries n (1 - q * 2 ^ k) - ratLogRem n (1 - q * 2 ^ k)) - k * log2Hi

/-- `log (q · 2^k) = log q + k log 2`. -/
private lemma log_mul_pow_two {q : ℚ} (hq : 0 < q) (k : ℕ) :
    Real.log ((q : ℝ) * 2 ^ k) = Real.log (q : ℝ) + k * Real.log 2 := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  rw [Real.log_mul (ne_of_gt hqR) (by positivity), Real.log_pow]

/-- **Soundness of the upper bound.** -/
theorem real_log_le_logUpper {q : ℚ} (hq : 0 < q) {k n : ℕ}
    (hlt : q * 2 ^ k < 2) : Real.log (q : ℝ) ≤ ((logUpper q k n : ℚ) : ℝ) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hshift : (0 : ℝ) < (q : ℝ) * 2 ^ k := by positivity
  have hlt2 : (q : ℝ) * 2 ^ k < 2 := by exact_mod_cast hlt
  have hser := (log_mem_Icc (q := (q : ℝ) * 2 ^ k) hshift hlt2 n).2
  have hsplit := log_mul_pow_two hq k
  have hl2 : (k : ℝ) * (log2Lo : ℝ) ≤ (k : ℝ) * Real.log 2 :=
    mul_le_mul_of_nonneg_left log2Lo_lt.le (by positivity)
  have hcast : ((1 : ℝ) - (q : ℝ) * 2 ^ k) = (((1 - q * 2 ^ k : ℚ) : ℝ)) := by
    push_cast
    ring
  rw [hcast] at hser
  unfold logUpper
  push_cast
  rw [ratLogSeries_cast, ratLogRem_cast]
  linarith [hser, hsplit, hl2]

/-- **Soundness of the lower bound.** -/
theorem logLower_le_real_log {q : ℚ} (hq : 0 < q) {k n : ℕ}
    (hlt : q * 2 ^ k < 2) : ((logLower q k n : ℚ) : ℝ) ≤ Real.log (q : ℝ) := by
  have hqR : (0 : ℝ) < (q : ℝ) := by exact_mod_cast hq
  have hshift : (0 : ℝ) < (q : ℝ) * 2 ^ k := by positivity
  have hlt2 : (q : ℝ) * 2 ^ k < 2 := by exact_mod_cast hlt
  have hser := (log_mem_Icc (q := (q : ℝ) * 2 ^ k) hshift hlt2 n).1
  have hsplit := log_mul_pow_two hq k
  have hl2 : (k : ℝ) * Real.log 2 ≤ (k : ℝ) * (log2Hi : ℝ) :=
    mul_le_mul_of_nonneg_left lt_log2Hi.le (by positivity)
  have hcast : ((1 : ℝ) - (q : ℝ) * 2 ^ k) = (((1 - q * 2 ^ k : ℚ) : ℝ)) := by
    push_cast
    ring
  rw [hcast] at hser
  unfold logLower
  push_cast
  rw [ratLogSeries_cast, ratLogRem_cast]
  linarith [hser, hsplit, hl2]

end Spin.Numeric
