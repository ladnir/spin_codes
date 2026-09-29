/-
One step of the IMT recurrence.

    Y_i = X_i + A Q_i,      Q_{i+1} = M_i Q_i + C X_i

with `M_i` the transvection and `X_i` the step input.  Combining the
transvection law with the syndrome shift `C X` gives the row of `T_j(z)`
leaving a state, and the combination is *exact*:

    2 ∑_x ∑_p wgt(x) f(M_p q + C x)
      = |u^⊥| · ( M ∑_x wgt(x) f(q + C x)
                  + ∑_x wgt(x) ( ∑_w f(w) − f(C x) ) ).

The second bracket is the refresh branch: flat over every state, with the
single state `C x` removed, because the refresh lands on `w + C x` for `w`
uniform over the *nonzero* states.

Two consequences are the paper's own entries.  Taking `f` to be the indicator
of zero gives the zero column, `½ L_j(q;z) + (1/2M) E[z^{wt} · 1(C x ≠ 0)]` —
the second factor is where "termination also requires a nonzero syndrome,
giving the minimum in the zero column" comes from.  And from state zero the
transvection is the identity, so that row carries no `½` at all, which is the
appendix's "the zero row is exact".

Nothing here needs `A` or `C` to be linear, or even to be anything but
functions: `A` enters only through the weight `wgt`, and `C` only as the shift.
Input weights are left as an arbitrary `wgt` over an arbitrary finite set of
inputs, so no normalisation and no probability measure is involved.
-/
import SpinCodes.Structured.Transfer

namespace Spin

open Finset
open scoped symmDiff

variable {s t : ℕ}

/-! ## Reindexing by a shift -/

lemma sum_symmDiff_shift (c : Finset (Fin s)) (f : Finset (Fin s) → ℝ) :
    ∑ w : Finset (Fin s), f (w ∆ c) = ∑ w : Finset (Fin s), f w := by
  classical
  refine Finset.sum_nbij' (fun w => w ∆ c) (fun w => w ∆ c) ?_ ?_ ?_ ?_ ?_
  · intro w _; exact Finset.mem_univ _
  · intro w _; exact Finset.mem_univ _
  · intro w _; exact symmDiff_cancel_tail w c
  · intro w _; exact symmDiff_cancel_tail w c
  · intro w _; rfl

/-- The refresh branch reaches every state except the shift itself. -/
lemma sum_nonzero_shift (c : Finset (Fin s)) (f : Finset (Fin s) → ℝ) :
    ∑ w ∈ nonzeroStates s, f (w ∆ c) = (∑ w : Finset (Fin s), f w) - f c := by
  classical
  have h1 : ∑ w : Finset (Fin s), f (w ∆ c)
      = f ((∅ : Finset (Fin s)) ∆ c)
        + ∑ w ∈ (univ : Finset (Finset (Fin s))).erase ∅, f (w ∆ c) :=
    (Finset.add_sum_erase _ (fun w => f (w ∆ c)) (Finset.mem_univ _)).symm
  have h3 : (∅ : Finset (Fin s)) ∆ c = c := by
    rw [show (∅ : Finset (Fin s)) = ⊥ from rfl, bot_symmDiff]
  rw [h3] at h1
  have h2 := sum_symmDiff_shift c f
  rw [nonzeroStates_eq_erase]
  linarith

/-! ## The step out of a nonzero state -/

/-- **One step from a nonzero state**, exactly. -/
theorem step_sum {q : Finset (Fin s)} (hq : q ≠ ∅)
    (Xs : Finset (Finset (Fin t))) (wgt : Finset (Fin t) → ℝ)
    (Cmap : Finset (Fin t) → Finset (Fin s)) (f : Finset (Fin s) → ℝ) :
    2 * ∑ x ∈ Xs, ∑ p ∈ transPairs s, wgt x * f (act p.1 p.2 q ∆ Cmap x)
      = ((perp q).card : ℝ)
        * (((nonzeroStates s).card : ℝ) * ∑ x ∈ Xs, wgt x * f (q ∆ Cmap x)
            + ∑ x ∈ Xs, wgt x * ((∑ w : Finset (Fin s), f w) - f (Cmap x))) := by
  classical
  have hper : ∀ x ∈ Xs,
      2 * ∑ p ∈ transPairs s, wgt x * f (act p.1 p.2 q ∆ Cmap x)
        = ((perp q).card : ℝ)
          * (wgt x * (((nonzeroStates s).card : ℝ) * f (q ∆ Cmap x)
              + ((∑ w : Finset (Fin s), f w) - f (Cmap x)))) := by
    intro x _
    have h := sum_act_eq hq (fun w => f (w ∆ Cmap x))
    rw [sum_nonzero_shift] at h
    calc 2 * ∑ p ∈ transPairs s, wgt x * f (act p.1 p.2 q ∆ Cmap x)
        = wgt x * (2 * ∑ p ∈ transPairs s, f (act p.1 p.2 q ∆ Cmap x)) := by
          rw [← Finset.mul_sum]; ring
      _ = wgt x * (((perp q).card : ℝ)
            * (((nonzeroStates s).card : ℝ) * f (q ∆ Cmap x)
                + ((∑ w : Finset (Fin s), f w) - f (Cmap x)))) := by rw [h]
      _ = _ := by ring
  rw [Finset.mul_sum, Finset.sum_congr rfl hper, ← Finset.mul_sum]
  congr 1
  rw [Finset.mul_sum, ← Finset.sum_add_distrib]
  exact Finset.sum_congr rfl fun x _ => by ring

/-- **The zero column.**  The lazy branch terminates exactly when `C x = q`,
and the refresh branch terminates only when the syndrome is nonzero. -/
theorem step_zero_column {q : Finset (Fin s)} (hq : q ≠ ∅)
    (Xs : Finset (Finset (Fin t))) (wgt : Finset (Fin t) → ℝ)
    (Cmap : Finset (Fin t) → Finset (Fin s)) :
    2 * ∑ x ∈ Xs, ∑ p ∈ transPairs s,
        wgt x * (if act p.1 p.2 q ∆ Cmap x = ∅ then 1 else 0)
      = ((perp q).card : ℝ)
        * (((nonzeroStates s).card : ℝ)
              * ∑ x ∈ Xs, wgt x * (if Cmap x = q then 1 else 0)
            + ∑ x ∈ Xs, wgt x * (if Cmap x = ∅ then 0 else 1)) := by
  classical
  have h := step_sum hq Xs wgt Cmap (fun w => if w = ∅ then (1 : ℝ) else 0)
  have htotal : ∑ w : Finset (Fin s), (if w = ∅ then (1 : ℝ) else 0) = 1 := by
    rw [Finset.sum_ite_eq' univ (∅ : Finset (Fin s)) (fun _ => (1 : ℝ))]
    simp
  have hlazy : ∀ x, (if q ∆ Cmap x = ∅ then (1 : ℝ) else 0)
      = (if Cmap x = q then 1 else 0) := by
    intro x
    by_cases hc : Cmap x = q
    · rw [if_pos hc, if_pos (by rw [hc]; simpa using symmDiff_self q)]
    · rw [if_neg hc, if_neg]
      intro hcon
      exact hc (symmDiff_eq_bot.mp (by simpa using hcon)).symm
  have hrefresh : ∀ x, (1 : ℝ) - (if Cmap x = ∅ then (1 : ℝ) else 0)
      = (if Cmap x = ∅ then 0 else 1) := by
    intro x
    by_cases hc : Cmap x = ∅ <;> simp [hc]
  rw [htotal] at h
  have hA : ∑ x ∈ Xs, wgt x * (if q ∆ Cmap x = ∅ then (1 : ℝ) else 0)
      = ∑ x ∈ Xs, wgt x * (if Cmap x = q then (1 : ℝ) else 0) :=
    Finset.sum_congr rfl fun x _ => by rw [hlazy x]
  have hB : ∑ x ∈ Xs, wgt x * ((1 : ℝ) - (if Cmap x = ∅ then (1 : ℝ) else 0))
      = ∑ x ∈ Xs, wgt x * (if Cmap x = ∅ then (0 : ℝ) else 1) :=
    Finset.sum_congr rfl fun x _ => by rw [hrefresh x]
  rw [h, hA, hB]

/-- **The `min` in the zero column.**  The refresh term is bounded both by the
total weight (giving `m_{d,j}`) and by the count of nonzero syndromes (giving
`ν_j`), so it is bounded by their minimum. -/
theorem sum_refresh_le_min (Xs : Finset (Finset (Fin t))) (wgt : Finset (Fin t) → ℝ)
    (Cmap : Finset (Fin t) → Finset (Fin s))
    (hw0 : ∀ x ∈ Xs, 0 ≤ wgt x) (hw1 : ∀ x ∈ Xs, wgt x ≤ 1) :
    ∑ x ∈ Xs, wgt x * (if Cmap x = ∅ then (0 : ℝ) else 1)
      ≤ min (∑ x ∈ Xs, wgt x)
          ((Xs.filter (fun x => Cmap x ≠ ∅)).card : ℝ) := by
  classical
  refine le_min (Finset.sum_le_sum fun x hx => ?_) ?_
  · by_cases hc : Cmap x = ∅ <;> simp [hc, hw0 x hx]
  · have hcount : ((Xs.filter (fun x => Cmap x ≠ ∅)).card : ℝ)
        = ∑ x ∈ Xs, (if Cmap x = ∅ then (0 : ℝ) else 1) := by
      rw [Finset.card_filter]
      push_cast
      exact Finset.sum_congr rfl fun x _ => by by_cases hc : Cmap x = ∅ <;> simp [hc]
    rw [hcount]
    refine Finset.sum_le_sum fun x hx => ?_
    by_cases hc : Cmap x = ∅
    · simp [hc]
    · simp only [hc, if_false, mul_one]
      exact hw1 x hx

/-! ## The step out of zero

`M_i` fixes the zero state, so this row carries no `½`: it is the input law
pushed through `C`, exactly. -/

@[simp] lemma act_empty (u v : Finset (Fin s)) : act u v (∅ : Finset (Fin s)) = ∅ := by
  unfold act
  simp

/-- **The zero row is exact.** -/
theorem step_from_zero (Xs : Finset (Finset (Fin t))) (wgt : Finset (Fin t) → ℝ)
    (Cmap : Finset (Fin t) → Finset (Fin s)) (f : Finset (Fin s) → ℝ) :
    ∑ x ∈ Xs, ∑ p ∈ transPairs s,
        wgt x * f (act p.1 p.2 (∅ : Finset (Fin s)) ∆ Cmap x)
      = ((transPairs s).card : ℝ) * ∑ x ∈ Xs, wgt x * f (Cmap x) := by
  classical
  have hshift : ∀ x, (∅ : Finset (Fin s)) ∆ Cmap x = Cmap x := by
    intro x
    rw [show (∅ : Finset (Fin s)) = ⊥ from rfl, bot_symmDiff]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [Finset.sum_congr rfl (fun p _ => by rw [act_empty, hshift x]),
    Finset.sum_const, nsmul_eq_mul]

end Spin
