import Mathlib

/-!
# Kernel-cheap binomial coefficients for the tail certificates

The tail certificates check numerical inequalities involving `C(N, q)` with
`decide +kernel`, so the kernel must evaluate the binomial. The original
definition `N! / (q! (N-q)!)` makes the kernel build `N!` (about 2 Mbit at
`N = 2^17`) through `N` sequential products in every term check, which needs
~18 GB per check at K = 2^24. Here

  `C(N, q) = (∏_{i<k} (N - k + 1 + i)) / k!`,  `k = min(q, N - q)`,

and both products are evaluated by a balanced binary split (`rprod`), so the
kernel performs `k` multiplications on balanced operands.
-/

namespace Spin.Turbo

open Finset

/-- `rprod d a len = ∏_{i<len} (a+i)` whenever `len ≤ 2^d`, by balanced splitting. -/
def rprod : ℕ → ℕ → ℕ → ℕ
  | 0, a, len => if len = 0 then 1 else a
  | d + 1, a, len =>
    if len ≤ 1 then (if len = 0 then 1 else a)
    else rprod d a (len / 2) * rprod d (a + len / 2) (len - len / 2)

theorem rprod_eq : ∀ (d a len : ℕ), len ≤ 2 ^ d → rprod d a len = ∏ i ∈ range len, (a + i)
  | 0, a, len, h => by
    have : len ≤ 1 := by simpa using h
    interval_cases len <;> simp [rprod]
  | d + 1, a, len, h => by
    by_cases h1 : len ≤ 1
    · interval_cases len <;> simp [rprod]
    · have hm : len / 2 ≤ 2 ^ d := by
        rw [pow_succ] at h; omega
      have hn : len - len / 2 ≤ 2 ^ d := by
        rw [pow_succ] at h; omega
      rw [rprod, if_neg h1, rprod_eq d a (len / 2) hm, rprod_eq d (a + len / 2) (len - len / 2) hn]
      have hsplit : len = len / 2 + (len - len / 2) := by omega
      conv_rhs => rw [hsplit, prod_range_add]
      congr 1
      refine prod_congr rfl fun i _ => ?_
      ring

theorem descFactorial_eq_prod_shift (n k : ℕ) (hk : k ≤ n) :
    n.descFactorial k = ∏ i ∈ range k, (n - k + 1 + i) := by
  induction k with
  | zero => simp
  | succ k ih =>
    rw [Nat.descFactorial_succ, ih (by omega), prod_range_succ', mul_comm]
    congr 1
    · refine prod_congr rfl fun i _ => ?_
      omega
    · omega

/-- `min(q, N - q)`. -/
def binomK (N q : ℕ) : ℕ := if 2 * q ≤ N then q else N - q

/-- Kernel-cheap natural-number binomial `C(N, q)` (valid for `N ≤ 2^d`, `q ≤ N`). -/
def binomFast (N d q : ℕ) : ℕ :=
  rprod d (N - binomK N q + 1) (binomK N q) / rprod d 1 (binomK N q)

theorem binomFast_eq (N d q : ℕ) (hq : q ≤ N) (hd : N ≤ 2 ^ d) :
    binomFast N d q = Nat.choose N q := by
  have hk : binomK N q ≤ N := by unfold binomK; split_ifs <;> omega
  have hk2 : binomK N q ≤ 2 ^ d := le_trans hk hd
  have hfac : rprod d 1 (binomK N q) = (binomK N q).factorial := by
    rw [rprod_eq d 1 _ hk2, ← prod_range_add_one_eq_factorial]
    exact prod_congr rfl fun i _ => by ring
  unfold binomFast
  rw [rprod_eq d _ _ hk2, ← descFactorial_eq_prod_shift N _ hk, hfac,
    ← Nat.choose_eq_descFactorial_div_factorial]
  unfold binomK; split_ifs
  · rfl
  · exact Nat.choose_symm hq

end Spin.Turbo
