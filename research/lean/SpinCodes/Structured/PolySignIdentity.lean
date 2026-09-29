import SpinCodes.Structured.PolyIdentity
import SpinCodes.Structured.PolySignIdentityDefs

namespace Spin.Structured

theorem listEval_le_tailBound (p : List ℤ) {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    listEval p x ≤ (tailBound p : ℝ) := by
  cases p with
  | nil => simp [listEval, tailBound]
  | cons a p =>
    have hn : (0 : ℝ) ≤ (posSum p : ℝ) := by exact_mod_cast posSum_nonneg p
    have hh := mul_le_mul_of_nonneg_left (listEval_le_posSum p h0 h1) h0
    have hi := mul_le_mul_of_nonneg_right h1 hn
    simp only [listEval, tailBound, Int.cast_add]
    linarith

theorem listEval_nonpos_of_drop (p : List ℤ) (h : tailBound (dropInitialZeros p) ≤ 0)
    {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) : listEval p x ≤ 0 := by
  induction p with
  | nil => simp [listEval]
  | cons a p ih =>
    by_cases ha : a = 0
    · simp only [dropInitialZeros, ha, if_pos] at h
      have ht := ih h
      simpa only [listEval, ha, Int.cast_zero, zero_add] using
        mul_nonpos_of_nonneg_of_nonpos h0 ht
    · have hh := listEval_le_tailBound (a :: p) h0 h1
      simp only [dropInitialZeros, if_neg ha] at h
      exact hh.trans (by exact_mod_cast h)

namespace RatPoly

theorem checkNonpos_sound (p : RatPoly) (h : checkNonpos p = true)
    {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) : p.eval x ≤ 0 := by
  simp only [checkNonpos, Bool.and_eq_true, decide_eq_true_eq] at h
  exact div_nonpos_of_nonpos_of_nonneg (listEval_nonpos_of_drop p.num h.2 h0 h1)
    (Nat.cast_nonneg _)

theorem checkLE_sound (p q : RatPoly) (h : checkLE p q = true)
    {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) : p.eval x ≤ q.eval x := by
  simp only [checkLE, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨hp, hq⟩, hh⟩ := h
  have hv := checkNonpos_sound _ hh h0 h1
  rw [eval_add _ _ hp (mul_den_pos _ _ (by decide) hq), eval_mul, eval_constant] at hv
  norm_num at hv
  linarith

end RatPoly
end Spin.Structured
