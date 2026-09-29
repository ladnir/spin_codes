/-
Exact rational polynomial sign certificates (`eq:imt-sparse-collatz`).

The verifier's `sign_certificate` removes an exact zero at `x = 0`, then bounds
what remains by

    upper := c₀ + Σ_{i≥1} max(cᵢ, 0),

and succeeds when `upper < 0`.  The paper describes exactly this: "After
removing that factor, its constant coefficient plus the sum of its positive
remaining coefficients is strictly negative."  All seven sparse residuals are
certified this way (`method: positive_power_tail`, `zero_order: 1`, degree 257).

The engine here is denominator-agnostic: coefficients are integers, because
clearing a common positive denominator does not change a sign.  Evaluation is
Horner, which is also how the bound is proved.

The concrete data is checked by the kernel via `decide`.  Compiled evaluation
is deliberately avoided: it would add `Lean.ofReduceBool` to the axiom
closure, which `scripts/check.sh` rejects.
-/
import Mathlib
import SpinCodes.Structured.PolyCertDefs

set_option linter.unusedSectionVars false

namespace Spin.Structured

/-- Horner evaluation of an integer coefficient list, lowest degree first. -/
def listEval : List ℤ → ℝ → ℝ
  | [], _ => 0
  | a :: t, x => (a : ℝ) + x * listEval t x

lemma posPart_nonneg (a : ℤ) : 0 ≤ posPart a := by
  unfold posPart; split <;> simp_all

lemma le_posPart (a : ℤ) : a ≤ posPart a := by
  unfold posPart
  split
  · exact le_rfl
  · omega

lemma posSum_nonneg (c : List ℤ) : 0 ≤ posSum c := by
  induction c with
  | nil => simp [posSum]
  | cons a t ih =>
      have := posPart_nonneg a
      simp only [posSum]
      omega

/-- On `[0,1]` the polynomial is at most the sum of the positive parts of its
coefficients. -/
lemma listEval_le_posSum (c : List ℤ) {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) :
    listEval c x ≤ (posSum c : ℝ) := by
  induction c with
  | nil => simp [listEval, posSum]
  | cons a t ih =>
      have hts : (0:ℝ) ≤ (posSum t : ℝ) := by exact_mod_cast posSum_nonneg t
      have hstep : x * listEval t x ≤ (posSum t : ℝ) := by
        by_cases hpos : 0 ≤ listEval t x
        · calc x * listEval t x ≤ 1 * listEval t x :=
                mul_le_mul_of_nonneg_right hx1 hpos
            _ = listEval t x := one_mul _
            _ ≤ (posSum t : ℝ) := ih
        · have hneg : listEval t x ≤ 0 := (not_le.mp hpos).le
          have : x * listEval t x ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hx0 hneg
          linarith
      have hmax : (a : ℝ) ≤ ((posPart a : ℤ) : ℝ) := by exact_mod_cast le_posPart a
      have hps : posSum (a :: t) = posPart a + posSum t := rfl
      rw [listEval, hps, Int.cast_add]
      linarith

/-- **Positive-power-tail certificate.**  If `c₀ + Σ_{i≥1} max(cᵢ,0) < 0` then
the polynomial is strictly negative on all of `[0,1]`. -/
theorem listEval_neg_of_tailBound (c : List ℤ) (h : tailBound c < 0)
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) : listEval c x < 0 := by
  cases c with
  | nil => simp [tailBound] at h
  | cons a t =>
      have hstep : x * listEval t x ≤ (posSum t : ℝ) := by
        have hts : (0:ℝ) ≤ (posSum t : ℝ) := by exact_mod_cast posSum_nonneg t
        by_cases hpos : 0 ≤ listEval t x
        · calc x * listEval t x ≤ 1 * listEval t x :=
                mul_le_mul_of_nonneg_right hx1 hpos
            _ = listEval t x := one_mul _
            _ ≤ (posSum t : ℝ) := listEval_le_posSum t hx0 hx1
        · have hneg : listEval t x ≤ 0 := (not_le.mp hpos).le
          have : x * listEval t x ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hx0 hneg
          linarith
      have hlt : ((a + posSum t : ℤ) : ℝ) < 0 := by
        have : (a + posSum t) < 0 := h
        exact_mod_cast this
      rw [listEval]
      push_cast at hlt
      linarith

/-- With an exact zero of order one removed, the original residual is strictly
negative on `(0,1]`.  This is the form the certificate is used in: the stored
list is the residual divided by `x`. -/
theorem residual_neg_of_tailBound (c : List ℤ) (h : tailBound c < 0)
    {x : ℝ} (hx0 : 0 < x) (hx1 : x ≤ 1) : x * listEval c x < 0 :=
  mul_neg_of_pos_of_neg hx0 (listEval_neg_of_tailBound c h hx0.le hx1)

/-- Clearing a common positive denominator preserves the sign, so an integer
coefficient list certifies the rational polynomial it was scaled from. -/
theorem rat_poly_neg_of_int_scaled (c : List ℤ) (D : ℝ) (hD : 0 < D)
    (p : ℝ → ℝ) (hp : ∀ x, D * p x = listEval c x) (h : tailBound c < 0)
    {x : ℝ} (hx0 : 0 ≤ x) (hx1 : x ≤ 1) : p x < 0 := by
  have h1 : D * p x < 0 := by
    rw [hp]; exact listEval_neg_of_tailBound c h hx0 hx1
  nlinarith

end Spin.Structured
