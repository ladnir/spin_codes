import SpinCodes.Numeric.BAExponent
namespace Spin.Numeric
private lemma half_expand_strict {a c : ℝ} (ha : 0 < a) (hc : 0 < c) (hac : a / 2 < c) :
    c * hEnt (a / (2 * c))
      = -(a / 2) * Real.log a + (a / 2) * Real.log 2 + c * Real.log c
        - (c - a / 2) * Real.log (c - a / 2) := by
  have hca : 0 < c - a / 2 := by linarith
  have h2c : (0 : ℝ) < 2 * c := by linarith
  have hsplit : Real.log (a / (2 * c)) = Real.log a - Real.log 2 - Real.log c := by
    rw [Real.log_div (ne_of_gt ha) (ne_of_gt h2c), Real.log_mul (by norm_num) (ne_of_gt hc)]
    ring
  have hone : (1 : ℝ) - a / (2 * c) = (c - a / 2) / c := by
    field_simp
  have hsplit2 : Real.log ((c - a / 2) / c) = Real.log (c - a / 2) - Real.log c :=
    Real.log_div (ne_of_gt hca) (ne_of_gt hc)
  have hmul : c * (a / (2 * c)) = a / 2 := by field_simp
  have hmul2 : c * ((c - a / 2) / c) = c - a / 2 := by field_simp
  unfold hEnt
  rw [hone, hsplit, hsplit2]
  have : c * (-(a / (2 * c)) * (Real.log a - Real.log 2 - Real.log c)
      - (c - a / 2) / c * (Real.log (c - a / 2) - Real.log c))
      = -(c * (a / (2 * c))) * (Real.log a - Real.log 2 - Real.log c)
        - (c * ((c - a / 2) / c)) * (Real.log (c - a / 2) - Real.log c) := by
    ring
  rw [this, hmul, hmul2]
  ring

lemma half_expand_closed {a c : ℝ} (ha : 0 < a) (hc : 0 < c) (hac : a / 2 ≤ c) :
    c * hEnt (a / (2 * c))
      = -(a / 2) * Real.log a + (a / 2) * Real.log 2 + c * Real.log c
        - (c - a / 2) * Real.log (c - a / 2) := by
  rcases eq_or_lt_of_le hac with he | ht
  · have heq : a = 2 * c := by linarith
    rw [heq, div_self (by positivity : (2 : ℝ) * c ≠ 0), hEnt_one]
    rw [Real.log_mul (by norm_num : (2 : ℝ) ≠ 0) hc.ne']
    have he2 : (2 : ℝ) * c / 2 = c := by ring
    rw [he2]
    simp
    ring
  · exact half_expand_strict ha hc ht

theorem piBA_eq_piEval_closed {a c : ℝ} (ha : 0 ≤ a) 
    (hlo : a / 2 ≤ c) (hhi : c ≤ 1 - a / 2) :
    piBA a c = piEval a c := by
  rcases eq_or_lt_of_le ha with hzero | hpos
  · -- `a = 0`: both sides collapse to zero.
    subst_vars
    simp only [piBA, piEval, hEnt_zero]
    rw [show (0 : ℝ) / (2 * c) = 0 by ring, show (0 : ℝ) / (2 * (1 - c)) = 0 by ring]
    simp only [hEnt_zero]
    norm_num
  · -- the generic case
    have hc : 0 < c := by linarith
    have hd : 0 < 1 - c := by linarith
    have hdlo : a / 2 ≤ 1 - c := by linarith
    have h1 := half_expand_closed hpos hc hlo
    have h2 := half_expand_closed hpos hd hdlo
    unfold piBA piEval
    rw [h1, h2]
    unfold hEnt
    rw [show (1 : ℝ) - c - a / 2 = (1 - c) - a / 2 by ring]
    ring


end Spin.Numeric
