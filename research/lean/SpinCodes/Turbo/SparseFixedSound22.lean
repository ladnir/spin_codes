import SpinCodes.Turbo.SparseFixed22
import SpinCodes.Turbo.SparseShellCoefficient22
import Mathlib.Tactic.IntervalCases

/-! Soundness of the remaining core-only checks in `SparseFixedDefs22`.

Audit: `Numeric.Fixed`/`Structured.DenseFourierExactCheck` keep exact integer
comparisons in a Mathlib-free definitions module and prove soundness through
positive common denominators; this file follows the same split.

* Degree-≤2 polynomial matrix products: one step of binary powering of the
  regional occupancy polynomial (`orderedMatrixDP_add` supplies the exact
  convolution).
* Regional normalization by choose(16R,d), d ≤ 2.
* The q1 row DP. `orderedMatrixDP` prepends a step on the left; the integer
  checker appends on the right. The bridge is `orderedMatrixDP_add` with the
  last factor of length one, so no physical step is reordered.
* The absolute lower power trace for low-tilt Chernoff denominators. -/
noncomputable section
namespace Spin.Turbo
open Finset Matrix

namespace SparseFixedSound22

theorem all11 (f : Nat → Bool) (h : sparseAll11Bool22 f = true) (i : Fin 11) :
    f i.val = true := by
  simp only [sparseAll11Bool22, Bool.and_eq_true] at h
  fin_cases i <;> tauto

theorem all3 (f : Nat → Bool) (h : sparseAll3Bool22 f = true) (d : ℕ) (hd : d ≤ 2) :
    f d = true := by
  simp only [sparseAll3Bool22, Bool.and_eq_true] at h
  interval_cases d <;> tauto

theorem sum11 (f : Nat → Int) : sparseSum11Int22 f = ∑ i : Fin 11, f i.val := by
  simp [sparseSum11Int22, Fin.sum_univ_succ]
  ring

theorem sum11_cast (f : Nat → Int) :
    ((sparseSum11Int22 f : Int) : ℝ) = ∑ i : Fin 11, (f i.val : ℝ) := by
  rw [sum11, Int.cast_sum]

theorem dot_le (s : ℝ) (hs : 0 < s) (x y X Y : Fin 11 → ℝ)
    (hx0 : ∀ k, 0 ≤ x k) (hy0 : ∀ k, 0 ≤ y k)
    (hx : ∀ k, x k ≤ X k / s) (hy : ∀ k, y k ≤ Y k / s) :
    ∑ k, x k * y k ≤ (∑ k, X k * Y k) / (s * s) := by
  have _hs := hs
  calc ∑ k, x k * y k ≤ ∑ k, (X k / s) * (Y k / s) :=
        Finset.sum_le_sum fun k _ =>
          mul_le_mul (hx k) (hy k) (hy0 k) ((hx0 k).trans (hx k))
    _ = (∑ k, X k * Y k) / (s * s) := by
        rw [Finset.sum_div]
        exact Finset.sum_congr rfl fun k _ => div_mul_div_comm _ _ _ _

theorem scale_close (s U I : ℤ) (hs : 0 < s) (A : ℝ)
    (hA : A ≤ (I : ℝ) / ((s : ℝ) * s)) (hI : I ≤ s * U) : A ≤ (U : ℝ) / s := by
  have hsR : (0 : ℝ) < s := by exact_mod_cast hs
  refine hA.trans ?_
  rw [div_le_div_iff₀ (mul_pos hsR hsR) hsR]
  have h : (I : ℝ) ≤ s * U := by exact_mod_cast hI
  nlinarith [mul_le_mul_of_nonneg_right h hsR.le]

/-! ## Degree-≤2 polynomial products -/

theorem antidiag1 {M : Type*} [AddCommMonoid M] (f : ℕ × ℕ → M) :
    ∑ p ∈ antidiagonal 1, f p = f (0, 1) + f (1, 0) := by
  rw [Finset.Nat.sum_antidiagonal_succ]
  simp

theorem antidiag2 {M : Type*} [AddCommMonoid M] (f : ℕ × ℕ → M) :
    ∑ p ∈ antidiagonal 2, f p = f (0, 2) + f (1, 1) + f (2, 0) := by
  rw [Finset.Nat.sum_antidiagonal_succ, Finset.Nat.sum_antidiagonal_succ]
  simp [add_assoc]

theorem polyMul_sound (scale : ℤ) (hs : 0 < scale) (a b u : Array ℤ)
    (hcheck : sparseIntPolyMulUpperCheck22 scale a b u = true)
    (F G : ℕ → Matrix (Fin 11) (Fin 11) ℝ)
    (hF0 : ∀ d i j, 0 ≤ F d i j) (hG0 : ∀ d i j, 0 ≤ G d i j)
    (hF : ∀ d, d ≤ 2 → ∀ i j, F d i j ≤ (sparseIntPolyMatrix22 a d i j : ℝ) / scale)
    (hG : ∀ d, d ≤ 2 → ∀ i j, G d i j ≤ (sparseIntPolyMatrix22 b d i j : ℝ) / scale)
    (d : ℕ) (hd : d ≤ 2) (i j : Fin 11) :
    (∑ p ∈ antidiagonal d, F p.1 * G p.2) i j ≤
      (sparseIntPolyMatrix22 u d i j : ℝ) / scale := by
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hc : sparseIntPolyProduct22 a b d i.val j.val ≤
      scale * sparseIntPolyMatrixGet22 u d i.val j.val :=
    of_decide_eq_true (all11 _ (all11 _ (all3 _ hcheck d hd) i) j)
  have hpair : ∀ r t, r ≤ 2 → t ≤ 2 → (F r * G t) i j ≤
      (sparseIntPolyDot22 a b r t i.val j.val : ℝ) / ((scale : ℝ) * scale) := by
    intro r t hr ht
    rw [Matrix.mul_apply]
    unfold sparseIntPolyDot22
    rw [sum11_cast]
    simp only [Int.cast_mul]
    exact dot_le _ hsR _ _ _ _ (fun k => hF0 r i k) (fun k => hG0 t k j)
      (fun k => hF r hr i k) (fun k => hG t ht k j)
  apply scale_close scale _ _ hs _ _ hc
  interval_cases d
  · rw [Finset.Nat.antidiagonal_zero, Finset.sum_singleton]
    exact hpair 0 0 (by norm_num) (by norm_num)
  · rw [antidiag1, Matrix.add_apply]
    show _ ≤ ((sparseIntPolyDot22 a b 0 1 i.val j.val +
      sparseIntPolyDot22 a b 1 0 i.val j.val : ℤ) : ℝ) / _
    rw [Int.cast_add, add_div]
    exact add_le_add (hpair 0 1 (by norm_num) (by norm_num))
      (hpair 1 0 (by norm_num) (by norm_num))
  · rw [antidiag2, Matrix.add_apply, Matrix.add_apply]
    show _ ≤ ((sparseIntPolyDot22 a b 0 2 i.val j.val +
      sparseIntPolyDot22 a b 1 1 i.val j.val +
      sparseIntPolyDot22 a b 2 0 i.val j.val : ℤ) : ℝ) / _
    rw [Int.cast_add, Int.cast_add, add_div, add_div]
    exact add_le_add (add_le_add (hpair 0 2 (by norm_num) (by norm_num))
      (hpair 1 1 (by norm_num) (by norm_num))) (hpair 2 0 (by norm_num) (by norm_num))

/-! ## Regional normalization -/

theorem regionDenominator_eq (R d : ℕ) (hd : d ≤ 2) :
    sparseRegionDenominator22 R d = Nat.choose (R * 16) d := by
  interval_cases d
  · show 1 = Nat.choose (R * 16) 0
    rw [Nat.choose_zero_right]
  · show 16 * R = Nat.choose (R * 16) 1
    rw [Nat.choose_one_right, Nat.mul_comm]
  · show (16 * R) * (16 * R - 1) / 2 = Nat.choose (R * 16) 2
    rw [Nat.choose_two_right, Nat.mul_comm R 16]

theorem normalize_sound (R : ℕ) (hR : 0 < R) (scale : ℤ) (hs : 0 < scale)
    (u result : Array ℤ)
    (hcheck : sparseIntRegionNormalizeCheck22 R u result = true)
    (X : ℕ → Matrix (Fin 11) (Fin 11) ℝ)
    (hX : ∀ d, d ≤ 2 → ∀ i j, X d i j ≤ (sparseIntPolyMatrix22 u d i j : ℝ) / scale)
    (d : ℕ) (hd : d ≤ 2) (i j : Fin 11) :
    ((1 / (Nat.choose (R * 16) d : ℝ)) • X d) i j ≤
      (sparseIntPolyMatrix22 result d i j : ℝ) / scale := by
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hc := of_decide_eq_true (all11 _ (all11 _ (all3 _ hcheck d hd) i) j)
  rw [regionDenominator_eq R d hd] at hc
  have hpos : (0 : ℝ) < (Nat.choose (R * 16) d : ℝ) := by
    exact_mod_cast Nat.choose_pos (by omega)
  have hcR : (sparseIntPolyMatrix22 u d i j : ℝ) ≤
      (Nat.choose (R * 16) d : ℝ) * (sparseIntPolyMatrix22 result d i j : ℝ) := by
    exact_mod_cast hc
  rw [Matrix.smul_apply, smul_eq_mul]
  calc 1 / (Nat.choose (R * 16) d : ℝ) * X d i j
      ≤ 1 / (Nat.choose (R * 16) d : ℝ) * ((sparseIntPolyMatrix22 u d i j : ℝ) / scale) :=
        mul_le_mul_of_nonneg_left (hX d hd i j) (by positivity)
    _ ≤ _ := by
        rw [div_mul_div_comm, one_mul, div_le_div_iff₀ (mul_pos hpos hsR) hsR]
        nlinarith [mul_le_mul_of_nonneg_right hcR hsR.le]

/-! ## q1 row DP -/

variable {I : Type*} [Fintype I] [DecidableEq I]

theorem listNonneg (ms : List (Matrix I I ℝ))
    (hm : ∀ M ∈ ms, ∀ i j, 0 ≤ M i j) : ∀ i j, 0 ≤ ms.prod i j := by
  induction ms with
  | nil =>
    intro i j
    simp only [List.prod_nil, Matrix.one_apply]
    split_ifs
    · exact zero_le_one
    · exact le_refl 0
  | cons M ms ih =>
    intro i j
    rw [List.prod_cons, Matrix.mul_apply]
    exact sum_nonneg (fun k _ => mul_nonneg (hm M (List.mem_cons_self) i k)
      (ih (fun N hN => hm N (List.mem_cons_of_mem _ hN)) k j))

theorem dp_one_bool (W : Bool → Matrix I I ℝ) (b : ℕ) :
    orderedMatrixDP sparseBoolDegree22 W 1 b =
      if b = 0 then W false else if b = 1 then W true else 0 := by
  rcases b with _ | _ | b <;> simp [orderedMatrixDP, sparseBoolDegree22]

theorem dp_succ_zero (W : Bool → Matrix I I ℝ) (n : ℕ) :
    orderedMatrixDP sparseBoolDegree22 W (n + 1) 0 =
      orderedMatrixDP sparseBoolDegree22 W n 0 * W false := by
  rw [orderedMatrixDP_add, Finset.Nat.antidiagonal_zero, Finset.sum_singleton, dp_one_bool]
  simp

theorem dp_succ_succ (W : Bool → Matrix I I ℝ) (n v : ℕ) :
    orderedMatrixDP sparseBoolDegree22 W (n + 1) (v + 1) =
      orderedMatrixDP sparseBoolDegree22 W n (v + 1) * W false +
        orderedMatrixDP sparseBoolDegree22 W n v * W true := by
  rw [orderedMatrixDP_add, Finset.Nat.sum_antidiagonal_succ']
  rcases v with _ | v
  · simp [dp_one_bool]
  · rw [Finset.Nat.sum_antidiagonal_succ']
    simp [dp_one_bool]

theorem q1Rows_sound (scale : ℤ) (hs : 0 < scale) (r0 r1 : Array ℤ)
    (W : Bool → Matrix (Fin 11) (Fin 11) ℝ) (hW0 : ∀ b i j, 0 ≤ W b i j)
    (hr0 : ∀ i j, W false i j ≤ (sparseIntMatrix22 r0 i j : ℝ) / scale)
    (hr1 : ∀ i j, W true i j ≤ (sparseIntMatrix22 r1 i j : ℝ) / scale)
    (rows : ℕ → Array ℤ) (N : ℕ)
    (hinit : sparseIntQ1InitialCheck22 scale (rows 0) = true)
    (hstep : ∀ n, n < N →
      sparseIntQ1StepUpperCheck22 scale r0 r1 (rows n) (rows (n + 1)) = true) :
    ∀ n, n ≤ N → ∀ v, v < 65 → ∀ j : Fin 11,
      orderedMatrixDP sparseBoolDegree22 W n v 0 j ≤
        (sparseIntRows22 (rows n) v j : ℝ) / scale := by
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  have hdp0 : ∀ n q i j, 0 ≤ orderedMatrixDP sparseBoolDegree22 W n q i j := by
    intro n q i j
    rw [← orderedMatrixCoefficient_eq_dp]
    unfold orderedMatrixCoefficient
    rw [Matrix.sum_apply]
    apply Finset.sum_nonneg
    intro xs _
    split_ifs
    · exact listNonneg _ (fun M hM => by
        obtain ⟨step, rfl⟩ := List.mem_ofFn.mp hM
        exact hW0 _) i j
    · exact le_rfl
  intro n
  induction n with
  | zero =>
    intro _ v hv j
    have hc := of_decide_eq_true
      (all11 _ (List.all_eq_true.mp hinit v (List.mem_range.mpr hv)) j)
    have hcR : ((if v = 0 ∧ j.val = 0 then scale else 0 : ℤ) : ℝ) ≤
        (sparseIntRows22 (rows 0) v j : ℝ) := by exact_mod_cast hc
    show orderedMatrixDP sparseBoolDegree22 W 0 v 0 j ≤ _
    simp only [orderedMatrixDP]
    by_cases hv0 : v = 0
    · subst hv0
      rw [if_pos rfl, Matrix.one_apply]
      by_cases hj : (0 : Fin 11) = j
      · subst hj
        rw [if_pos rfl, le_div_iff₀ hsR, one_mul]
        simpa using hcR
      · have hjv : ¬ (j.val = 0) := fun h => hj (Fin.ext h.symm)
        rw [if_neg hj]
        simp only [hjv, and_false, if_false, Int.cast_zero] at hcR
        exact div_nonneg hcR hsR.le
    · rw [if_neg hv0, Matrix.zero_apply]
      simp only [hv0, false_and, if_false, Int.cast_zero] at hcR
      exact div_nonneg hcR hsR.le
  | succ n ih =>
    intro hn v hv j
    have ihn := ih (by omega)
    have hc := of_decide_eq_true
      (all11 _ (List.all_eq_true.mp (hstep n (by omega)) v (List.mem_range.mpr hv)) j)
    apply scale_close scale _ _ hs _ _ hc
    rcases v with _ | v
    · rw [dp_succ_zero, Matrix.mul_apply, if_pos rfl, add_zero, sum11_cast]
      simp only [Int.cast_mul]
      exact dot_le _ hsR _ _ _ _ (fun k => hdp0 n 0 0 k) (fun k => hW0 false k j)
        (fun k => ihn 0 (by norm_num) k) (fun k => hr0 k j)
    · rw [dp_succ_succ, Matrix.add_apply, Matrix.mul_apply, Matrix.mul_apply,
        if_neg (Nat.succ_ne_zero v), Nat.add_sub_cancel, Int.cast_add, add_div,
        sum11_cast, sum11_cast]
      simp only [Int.cast_mul]
      exact add_le_add
        (dot_le _ hsR _ _ _ _ (fun k => hdp0 n (v + 1) 0 k) (fun k => hW0 false k j)
          (fun k => ihn (v + 1) hv k) (fun k => hr0 k j))
        (dot_le _ hsR _ _ _ _ (fun k => hdp0 n v 0 k) (fun k => hW0 true k j)
          (fun k => ihn v (by omega) k) (fun k => hr1 k j))

theorem q1Rows_terminal (scale : ℤ) (hs : 0 < scale) (r0 r1 : Array ℤ)
    (W : Bool → Matrix (Fin 11) (Fin 11) ℝ) (hW0 : ∀ b i j, 0 ≤ W b i j)
    (hr0 : ∀ i j, W false i j ≤ (sparseIntMatrix22 r0 i j : ℝ) / scale)
    (hr1 : ∀ i j, W true i j ≤ (sparseIntMatrix22 r1 i j : ℝ) / scale)
    (rows : ℕ → Array ℤ) (N : ℕ)
    (hinit : sparseIntQ1InitialCheck22 scale (rows 0) = true)
    (hstep : ∀ n, n < N →
      sparseIntQ1StepUpperCheck22 scale r0 r1 (rows n) (rows (n + 1)) = true)
    (v : ℕ) (hv : v < 65) :
    ∑ j : Fin 11, orderedMatrixCoefficient sparseBoolDegree22 W N v 0 j ≤
      ∑ j : Fin 11, (sparseIntRows22 (rows N) v j : ℝ) / scale := by
  rw [orderedMatrixCoefficient_eq_dp]
  exact Finset.sum_le_sum fun j _ =>
    q1Rows_sound scale hs r0 r1 W hW0 hr0 hr1 rows N hinit hstep N le_rfl v hv j

/-! ## Lower power trace -/

private def traceStep (e : ℕ) (step : Bool × ℤ) : ℕ := 2 * e + if step.1 then 1 else 0

private theorem power_trace_aux (scale base : ℤ) (hs : 0 < scale) (hb : 0 ≤ base) :
    ∀ (trace : List (Bool × ℤ)) (prev : ℤ) (e : ℕ), 0 ≤ prev →
      (prev : ℝ) / scale ≤ ((base : ℝ) / scale) ^ e →
      sparseIntPowerTraceCheck22 scale base prev trace = true →
      (sparsePowerTraceLast22 prev trace : ℝ) / scale ≤
        ((base : ℝ) / scale) ^ (trace.foldl traceStep e)
  | [], prev, e, _, h, _ => by simpa [sparsePowerTraceLast22] using h
  | (bit, next) :: tail, prev, e, _, h, hc => by
    simp only [sparseIntPowerTraceCheck22, Bool.and_eq_true, decide_eq_true_eq] at hc
    obtain ⟨⟨hp0, hn0, hstep⟩, hrest⟩ := hc
    have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
    have hx : (0 : ℝ) ≤ (base : ℝ) / scale := div_nonneg (by exact_mod_cast hb) hsR.le
    have hpp : ((prev : ℝ) / scale) ^ 2 ≤ (((base : ℝ) / scale) ^ e) ^ 2 :=
      pow_le_pow_left₀ (div_nonneg (by exact_mod_cast hp0) hsR.le) h 2
    have hnext : (next : ℝ) / scale ≤ ((base : ℝ) / scale) ^ (traceStep e (bit, next)) := by
      cases bit
      · simp only [Bool.false_eq_true, if_false] at hstep
        have hstepR : (scale : ℝ) * next ≤ (prev : ℝ) * prev := by exact_mod_cast hstep
        have h1 : (next : ℝ) / scale ≤ ((prev : ℝ) / scale) ^ 2 := by
          rw [div_pow, div_le_div_iff₀ hsR (by positivity)]
          nlinarith [mul_le_mul_of_nonneg_left hstepR hsR.le]
        calc (next : ℝ) / scale ≤ ((prev : ℝ) / scale) ^ 2 := h1
          _ ≤ (((base : ℝ) / scale) ^ e) ^ 2 := hpp
          _ = _ := by simp only [traceStep, Bool.false_eq_true, if_false]; ring
      · simp only [if_true] at hstep
        have hstepR : (scale : ℝ) * scale * next ≤ (prev : ℝ) * prev * base := by
          exact_mod_cast hstep
        have h1 : (next : ℝ) / scale ≤ ((prev : ℝ) / scale) ^ 2 * ((base : ℝ) / scale) := by
          rw [div_pow, div_mul_div_comm, div_le_div_iff₀ hsR (by positivity)]
          nlinarith [mul_le_mul_of_nonneg_left hstepR hsR.le]
        calc (next : ℝ) / scale ≤ ((prev : ℝ) / scale) ^ 2 * ((base : ℝ) / scale) := h1
          _ ≤ (((base : ℝ) / scale) ^ e) ^ 2 * ((base : ℝ) / scale) :=
            mul_le_mul_of_nonneg_right hpp hx
          _ = _ := by simp only [traceStep, if_true]; ring
    exact power_trace_aux scale base hs hb tail next (traceStep e (bit, next)) hn0 hnext hrest

/-- A passing certificate gives a strictly positive lower bound for
(base/scale)^exponent. -/
theorem powerCertificate_sound (scale base : ℤ) (exponent : ℕ)
    (trace : List (Bool × ℤ))
    (hcheck : sparseIntPowerCertificateCheck22 scale base exponent trace = true) :
    0 < (sparsePowerTraceLast22 scale trace : ℝ) / scale ∧
      (sparsePowerTraceLast22 scale trace : ℝ) / scale ≤
        ((base : ℝ) / scale) ^ exponent := by
  simp only [sparseIntPowerCertificateCheck22, Bool.and_eq_true, decide_eq_true_eq] at hcheck
  obtain ⟨⟨hs, hb, hexp, hpos⟩, htr⟩ := hcheck
  have hsR : (0 : ℝ) < scale := by exact_mod_cast hs
  refine ⟨div_pos (by exact_mod_cast hpos) hsR, ?_⟩
  have h := power_trace_aux scale base hs hb trace scale 0 hs.le
    (by rw [div_self hsR.ne', pow_zero]) htr
  rw [← hexp]
  exact h

end SparseFixedSound22

open SparseFixedSound22 in
#print axioms polyMul_sound
open SparseFixedSound22 in
#print axioms normalize_sound
open SparseFixedSound22 in
#print axioms q1Rows_terminal
open SparseFixedSound22 in
#print axioms powerCertificate_sound
end Spin.Turbo
