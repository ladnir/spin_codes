import Mathlib

/-! Kernel-cheap binomial for the tail certificates.

`binomQ` used `N! / (q! (N-q)!)`, which the kernel recomputes in every term check;
at N = 2^17 (K = 2^24) that is ~2 Mbit through 131072 sequential products (~18 GB
per check). Here `C(N,q) = (∏_{i<k} (N-k+1+i)) / k!` with `k = min(q, N-q)`, and
both products use a balanced binary split (`rprod`) so intermediate operands stay
balanced. Experiment file; not part of the build. -/

namespace BinomFast

open Finset

/-- `rprod d a len = ∏_{i<len} (a+i)` for `len ≤ 2^d`, by balanced splitting. -/
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

theorem descFactorial_eq_prod (n k : ℕ) (hk : k ≤ n) :
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
def kOf (q : ℕ) : ℕ := if 2 * q ≤ 131072 then q else 131072 - q

/-- Kernel-cheap `C(131072, q)`. -/
def binomQ (q : ℕ) : ℤ :=
  ((rprod 17 (131072 - kOf q + 1) (kOf q) / rprod 17 1 (kOf q) : ℕ) : ℤ)

theorem binomQ_nat (q : ℕ) (hq : q ≤ 131072) :
    rprod 17 (131072 - kOf q + 1) (kOf q) / rprod 17 1 (kOf q) = Nat.choose 131072 q := by
  have hk : kOf q ≤ 131072 := by unfold kOf; split_ifs <;> omega
  have hk2 : kOf q ≤ 2 ^ 17 := by norm_num; omega
  have hfac : rprod 17 1 (kOf q) = (kOf q).factorial := by
    rw [rprod_eq 17 1 _ hk2, ← prod_range_add_one_eq_factorial]
    exact prod_congr rfl fun i _ => by ring
  rw [rprod_eq 17 _ _ hk2, ← descFactorial_eq_prod 131072 _ hk, hfac,
    ← Nat.choose_eq_descFactorial_div_factorial]
  unfold kOf; split_ifs
  · rfl
  · exact Nat.choose_symm hq

theorem binomQ_eq (q : ℕ) (hq : q ≤ 131072) : (binomQ q : ℝ) = (Nat.choose 131072 q : ℝ) := by
  unfold binomQ
  rw [binomQ_nat q hq]
  simp

example : 0 < binomQ 3 := by decide +kernel
example : 0 < binomQ 77227 := by decide +kernel
example : 0 < binomQ 131072 := by decide +kernel

end BinomFast
