/-
The negative-binomial bound behind `eq:structured-ba-sparse-moment`.

The paper uses

    Σ_{h ≥ ℓ} C(h-1, ℓ-1) ζ^h = (ζ/(1-ζ))^ℓ

to evaluate the `ζ`-moment of the accumulator output weight.  What the proof
actually needs is the *inequality* for partial sums — the moment is a finite
sum over `h ≤ b` — so no convergence argument is required, and the bound is
proved by a nested induction instead.

Writing the sum as `Σ_{j < H} C(j, k) z^{j+1}` (so `h = j+1`, `k = ℓ-1`)
removes every truncated subtraction; the `h = 0` term of the paper's form is
exactly where `ℕ` truncation would otherwise inject a spurious `C(0,0) = 1`.

The step is Pascal:

    f k (H+1) = z · (f (k-1) H + f k H),

and the bound `A_k = (z/(1-z))^{k+1}` satisfies `z(A_{k-1} + A_k) = A_k`
exactly, so the induction closes with no slack lost.
-/
import Mathlib

set_option linter.unusedSectionVars false

namespace Spin.Structured

open Finset

variable {z : ℝ}

/-- Partial sums of the negative-binomial series, reindexed to avoid `ℕ`
subtraction: `negBinomSum z k H = Σ_{j<H} C(j,k) z^{j+1}`. -/
noncomputable def negBinomSum (z : ℝ) (k H : ℕ) : ℝ :=
  ∑ j ∈ range H, (j.choose k : ℝ) * z ^ (j + 1)

lemma negBinomSum_zero_right (z : ℝ) (k : ℕ) : negBinomSum z k 0 = 0 := by
  simp [negBinomSum]

lemma negBinomSum_nonneg (hz : 0 ≤ z) (k H : ℕ) : 0 ≤ negBinomSum z k H :=
  Finset.sum_nonneg fun j _ => by positivity

/-- Pascal's rule, transported to the partial sums. -/
lemma negBinomSum_succ (z : ℝ) (k H : ℕ) :
    negBinomSum z (k + 1) (H + 1)
      = z * (negBinomSum z k H + negBinomSum z (k + 1) H) := by
  unfold negBinomSum
  rw [Finset.sum_range_succ']
  have hzero : ((0 : ℕ).choose (k + 1) : ℝ) * z ^ (0 + 1) = 0 := by simp
  rw [hzero, add_zero, mul_add, Finset.mul_sum, Finset.mul_sum,
    ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  have hp : ((j + 1).choose (k + 1) : ℝ) = (j.choose k : ℝ) + (j.choose (k + 1) : ℝ) := by
    rw [Nat.choose_succ_succ]
    push_cast
    ring
  rw [hp]
  ring

/-- The geometric base case, `k = 0`. -/
lemma negBinomSum_zero_le (hz0 : 0 ≤ z) (hz1 : z < 1) (H : ℕ) :
    negBinomSum z 0 H ≤ z / (1 - z) := by
  have h1z : 0 < 1 - z := by linarith
  have hgeom : negBinomSum z 0 H = z * ∑ j ∈ range H, z ^ j := by
    unfold negBinomSum
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun j _ => by simp [pow_succ]; ring
  rw [hgeom]
  have hsum : ∑ j ∈ range H, z ^ j ≤ 1 / (1 - z) := by
    have hmul : (∑ j ∈ range H, z ^ j) * (1 - z) = 1 - z ^ H := by
      have := geom_sum_mul z H
      linear_combination -this
    rw [le_div_iff₀ h1z, hmul]
    have : 0 ≤ z ^ H := pow_nonneg hz0 H
    linarith
  calc z * ∑ j ∈ range H, z ^ j ≤ z * (1 / (1 - z)) :=
        mul_le_mul_of_nonneg_left hsum hz0
    _ = z / (1 - z) := by ring

/-- **The negative-binomial bound.**

`Σ_{j<H} C(j,k) z^{j+1} ≤ (z/(1-z))^{k+1}` for `0 ≤ z < 1`, uniformly in `H`.
In the paper's indexing this is `Σ_{h ≥ ℓ} C(h-1,ℓ-1) z^h ≤ (z/(1-z))^ℓ`. -/
theorem negBinomSum_le (hz0 : 0 ≤ z) (hz1 : z < 1) (k H : ℕ) :
    negBinomSum z k H ≤ (z / (1 - z)) ^ (k + 1) := by
  have h1z : 0 < 1 - z := by linarith
  have hA : (0 : ℝ) ≤ z / (1 - z) := by positivity
  induction k generalizing H with
  | zero =>
      simpa using negBinomSum_zero_le hz0 hz1 H
  | succ k ih =>
      induction H with
      | zero =>
          rw [negBinomSum_zero_right]
          positivity
      | succ n ihH =>
          rw [negBinomSum_succ]
          have hA1 : z * (1 + z / (1 - z)) = z / (1 - z) := by
            field_simp
            ring
          have hstep : z * ((z / (1 - z)) ^ (k + 1) + (z / (1 - z)) ^ (k + 1 + 1))
              = (z / (1 - z)) ^ (k + 1 + 1) := by
            have hsplit : (z / (1 - z)) ^ (k + 1 + 1)
                = (z / (1 - z)) ^ (k + 1) * (z / (1 - z)) := by ring
            rw [hsplit]
            calc z * ((z / (1 - z)) ^ (k + 1) + (z / (1 - z)) ^ (k + 1) * (z / (1 - z)))
                = (z / (1 - z)) ^ (k + 1) * (z * (1 + z / (1 - z))) := by ring
              _ = (z / (1 - z)) ^ (k + 1) * (z / (1 - z)) := by rw [hA1]
          calc z * (negBinomSum z k n + negBinomSum z (k + 1) n)
              ≤ z * ((z / (1 - z)) ^ (k + 1) + (z / (1 - z)) ^ (k + 1 + 1)) := by
                refine mul_le_mul_of_nonneg_left (add_le_add (ih n) ihH) hz0
            _ = (z / (1 - z)) ^ (k + 1 + 1) := hstep

/-- The paper's form: with `ℓ ≥ 1` and `h = j+1`, the moment sum over any range
of output weights is at most `(ζ/(1-ζ))^ℓ`. -/
theorem moment_sum_le (hz0 : 0 ≤ z) (hz1 : z < 1) {ℓ : ℕ} (hℓ : 1 ≤ ℓ) (H : ℕ) :
    ∑ j ∈ range H, (j.choose (ℓ - 1) : ℝ) * z ^ (j + 1) ≤ (z / (1 - z)) ^ ℓ := by
  have h := negBinomSum_le hz0 hz1 (ℓ - 1) H
  rwa [show ℓ - 1 + 1 = ℓ from by omega] at h


/-! ## Sanity check

At `z = 1/2` the bound is exactly `1`, which is the statement that a negative
binomial distribution sums to one — so the bound is tight in the limit and the
constant is not an artefact. -/

example (k H : ℕ) : negBinomSum (1 / 2 : ℝ) k H ≤ 1 := by
  have h := negBinomSum_le (z := (1 / 2 : ℝ)) (by norm_num) (by norm_num) k H
  norm_num at h
  exact h

end Spin.Structured
