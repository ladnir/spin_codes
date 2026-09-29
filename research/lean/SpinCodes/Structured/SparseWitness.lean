import SpinCodes.Structured.Occupation

/-! The exact affine column used in the sparse Collatz certificate. -/
noncomputable section
namespace Spin.Imt.Occupation.Sparse

open Finset

def count : Fin 5 → ℝ := ![5166, 110288, 293455, 110128, 5250]
def level : Fin 5 → ℕ := ![48, 56, 64, 72, 80]
def correction : Fin 5 → ℝ := ![
  116184452864 / 2257055535,
  92653723520 / 3613910291,
  -(26853143552 / 769273207925),
  -(461216096896 / 18043337105),
  -(351546952448 / 6881266875)]

def witness (α : ℝ) : Coords 5 where
  Z := 1
  D := (1 + 900 * α) / 1024
  S := fun i => (1 + correction i * α) / 1024

theorem count_pos (i : Fin 5) : 0 < count i := by
  fin_cases i <;> norm_num [count]

theorem count_sum : ∑ i, count i = 524287 := by
  norm_num [count, Fin.sum_univ_succ]

theorem correction_bounds (i : Fin 5) : -900 ≤ correction i ∧ correction i ≤ 900 := by
  fin_cases i <;> norm_num [correction]

theorem correction_mean : ∑ i, count i * correction i = 0 := by
  norm_num [count, correction, Fin.sum_univ_succ]

theorem liveAverage_witness (α : ℝ) :
    liveAverage count 524287 (witness α) = 1 / 1024 := by
  have hs : (∑ i, count i * (witness α).S i) =
      ((∑ i, count i) + (∑ i, count i * correction i) * α) / 1024 := by
    simp only [witness, Finset.sum_mul]
    rw [← Finset.sum_add_distrib, Finset.sum_div]
    apply Finset.sum_congr rfl
    intro i _
    ring
  rw [liveAverage, hs, count_sum, correction_mean]
  norm_num

theorem witness_lower {α : ℝ} (h0 : 0 ≤ α) (h1 : α ≤ 1 / 10000) :
    1 / 2048 ≤ (witness α).Z ∧ 1 / 2048 ≤ (witness α).D ∧
      ∀ i, 1 / 2048 ≤ (witness α).S i := by
  refine ⟨by norm_num [witness], ?_, fun i => ?_⟩
  · dsimp [witness]
    linarith
  · have hl := (correction_bounds i).1
    have hm := mul_nonneg (by linarith : 0 ≤ correction i + 900) h0
    dsimp [witness]
    nlinarith

theorem witness_nonneg {α : ℝ} (h0 : 0 ≤ α) (h1 : α ≤ 1 / 10000) :
    (witness α).Nonneg := by
  obtain ⟨hz, hd, hs⟩ := witness_lower h0 h1
  exact ⟨by linarith, by linarith, fun i => by linarith [hs i]⟩

theorem live_pair_nonneg {α : ℝ} (h0 : 0 ≤ α) :
    0 ≤ (witness α).D + liveAverage count 524287 (witness α) := by
  rw [liveAverage_witness]
  dsimp [witness]
  positivity

end Spin.Imt.Occupation.Sparse


