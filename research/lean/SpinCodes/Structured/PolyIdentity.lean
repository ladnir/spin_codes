import SpinCodes.Structured.PolyCert
import SpinCodes.Structured.PolyIdentityDefs

/-! Soundness of polynomial identity certificates. A proposed coefficient
list is accepted only after exact integer arithmetic verifies its identity.
Common denominators need not be reduced or chosen in any particular way. -/

namespace Spin.Structured

theorem polyFactorial_eq (n : ℕ) : polyFactorial n = n.factorial := by
  induction n with
  | zero => rfl
  | succ n ih => simp [polyFactorial, ih, Nat.factorial_succ]

theorem polyChoose_eq (n k : ℕ) : polyChoose n k = n.choose k := by
  unfold polyChoose
  split
  · simp only [polyFactorial_eq]
    exact (Nat.choose_eq_factorial_div_factorial ‹k ≤ n›).symm
  · exact (Nat.choose_eq_zero_of_lt (by omega)).symm

theorem listEval_polyAdd (p q : List ℤ) (x : ℝ) :
    listEval (polyAdd p q) x = listEval p x + listEval q x := by
  induction p generalizing q with
  | nil => simp [polyAdd, listEval]
  | cons a p ih =>
    cases q with
    | nil => simp [polyAdd, listEval]
    | cons b q => simp [polyAdd, listEval, ih, Int.cast_add]; ring

theorem listEval_polyScale (a : ℤ) (p : List ℤ) (x : ℝ) :
    listEval (polyScale a p) x = a * listEval p x := by
  induction p with
  | nil => simp [polyScale, listEval]
  | cons b p ih => simp [polyScale, listEval, ih, Int.cast_mul]; ring

theorem listEval_polyMul (p q : List ℤ) (x : ℝ) :
    listEval (polyMul p q) x = listEval p x * listEval q x := by
  induction p with
  | nil => simp [polyMul, listEval]
  | cons a p ih =>
    simp [polyMul, listEval_polyAdd, listEval_polyScale, listEval, ih]
    ring

theorem listEval_zero_of_polyZero (p : List ℤ) (h : polyZero p = true) (x : ℝ) :
    listEval p x = 0 := by
  induction p with
  | nil => rfl
  | cons a p ih =>
    simp only [polyZero, Bool.and_eq_true, beq_iff_eq] at h
    simp [listEval, h.1, ih h.2]

theorem listEval_eq_of_polyEq (p q : List ℤ) (h : polyEq p q = true) (x : ℝ) :
    listEval p x = listEval q x := by
  induction p generalizing q with
  | nil => exact (listEval_zero_of_polyZero q h x).symm
  | cons a p ih =>
    cases q with
    | nil =>
      simp only [polyEq, Bool.and_eq_true, beq_iff_eq] at h
      simp [listEval, h.1, listEval_zero_of_polyZero p h.2 x]
    | cons b q =>
      simp only [polyEq, Bool.and_eq_true, beq_iff_eq] at h
      simp [listEval, h.1, ih q h.2]

namespace RatPoly

noncomputable def eval (p : RatPoly) (x : ℝ) : ℝ := listEval p.num x / p.den

theorem add_den_pos (p q : RatPoly) (hp : 0 < p.den) (hq : 0 < q.den) :
    0 < (add p q).den := Nat.lcm_pos hp hq

theorem mul_den_pos (p q : RatPoly) (hp : 0 < p.den) (hq : 0 < q.den) :
    0 < (mul p q).den := Nat.mul_pos hp hq

theorem eval_add (p q : RatPoly) (hp : 0 < p.den) (hq : 0 < q.den) (x : ℝ) :
    (add p q).eval x = p.eval x + q.eval x := by
  have hl := Nat.lcm_pos hp hq
  have hdp : ((Nat.lcm p.den q.den / p.den : ℕ) : ℝ) * p.den = Nat.lcm p.den q.den := by
    exact_mod_cast Nat.div_mul_cancel (Nat.dvd_lcm_left p.den q.den)
  have hdq : ((Nat.lcm p.den q.den / q.den : ℕ) : ℝ) * q.den = Nat.lcm p.den q.den := by
    exact_mod_cast Nat.div_mul_cancel (Nat.dvd_lcm_right p.den q.den)
  have hp' : (p.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hp)
  have hq' : (q.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hl' : (Nat.lcm p.den q.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hl)
  simp only [add, eval, listEval_polyAdd, listEval_polyScale, Int.cast_natCast]
  field_simp
  linear_combination (listEval p.num x * (q.den : ℝ)) * hdp +
    (listEval q.num x * (p.den : ℝ)) * hdq

theorem eval_mul (p q : RatPoly) (x : ℝ) :
    (mul p q).eval x = p.eval x * q.eval x := by
  simp only [mul, eval, listEval_polyMul, Nat.cast_mul]
  ring

@[simp] theorem eval_constant (a : ℤ) (d : ℕ) (x : ℝ) :
    (constant a d).eval x = (a : ℝ) / d := by
  simp [constant, eval, listEval]

@[simp] theorem eval_variable (x : ℝ) : X.eval x = x := by
  simp [X, eval, listEval]

theorem pow_den_pos (p : RatPoly) (hp : 0 < p.den) (n : ℕ) :
    0 < (p.pow n).den := by
  induction n with
  | zero => norm_num [pow, constant]
  | succ n ih => exact mul_den_pos _ _ ih hp

theorem eval_pow (p : RatPoly) (n : ℕ) (x : ℝ) : (p.pow n).eval x = p.eval x ^ n := by
  induction n with
  | zero => simp [pow]
  | succ n ih => simp [pow, eval_mul, ih, pow_succ]

theorem sum_den_pos (ps : List RatPoly) (h : ∀ p ∈ ps, 0 < p.den) :
    0 < (sum ps).den := by
  induction ps with
  | nil => norm_num [sum, constant]
  | cons p ps ih =>
    exact add_den_pos _ _ (h _ (List.mem_cons_self ..))
      (ih fun q hq => h q (List.mem_cons_of_mem _ hq))

theorem eval_sum (ps : List RatPoly) (h : ∀ p ∈ ps, 0 < p.den) (x : ℝ) :
    (sum ps).eval x = (ps.map (fun p => p.eval x)).sum := by
  induction ps with
  | nil => simp [sum]
  | cons p ps ih =>
    have ht : ∀ q ∈ ps, 0 < q.den := fun q hq => h q (List.mem_cons_of_mem _ hq)
    simp only [sum, List.map_cons, List.sum_cons]
    rw [eval_add _ _ (h _ (List.mem_cons_self ..)) (sum_den_pos ps ht), ih ht]

theorem checkEq_sound (p q : RatPoly) (h : checkEq p q = true) (x : ℝ) :
    p.eval x = q.eval x := by
  simp only [checkEq, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨hp, hq⟩, he⟩ := h
  have hp' : (p.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hp)
  have hq' : (q.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hh := listEval_eq_of_polyEq _ _ he x
  simp only [listEval_polyScale, Int.cast_natCast] at hh
  unfold eval
  field_simp
  nlinarith [hh]

theorem checkAdd_sound (p q r : RatPoly) (h : checkAdd p q r = true) (x : ℝ) :
    r.eval x = p.eval x + q.eval x := by
  simp only [checkAdd, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨hp, hq⟩, hr⟩, he⟩ := h
  have hp' : (p.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hp)
  have hq' : (q.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hr' : (r.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hr)
  have hh := listEval_eq_of_polyEq _ _ he x
  simp only [listEval_polyAdd, listEval_polyScale, Int.cast_mul, Int.cast_natCast] at hh
  unfold eval
  field_simp
  nlinarith [hh]

theorem checkMul_sound (p q r : RatPoly) (h : checkMul p q r = true) (x : ℝ) :
    r.eval x = p.eval x * q.eval x := by
  simp only [checkMul, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨hp, hq⟩, hr⟩, he⟩ := h
  have hp' : (p.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hp)
  have hq' : (q.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hq)
  have hr' : (r.den : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hr)
  have hh := listEval_eq_of_polyEq _ _ he x
  simp only [listEval_polyMul, listEval_polyScale, Int.cast_mul, Int.cast_natCast] at hh
  unfold eval
  field_simp
  nlinarith [hh]

end RatPoly
end Spin.Structured


