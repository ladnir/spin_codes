import SpinCodes.UniformMixer.ScalableTransferPaperOrientation

/-!
# Sparse-occupation Collatz witness

The paper's rational sparse witness is oriented toward the row transfer.
Only the three componentwise numerical inequalities remain as assumptions
of the final theorem in this file. No numerical inequality is asserted here.
-/

noncomputable section

namespace Spin.UniformMixer

def sparseBeta (α : ℝ) : ℝ := (4 / 5) * α
def sparseTilt (α : ℝ) : ℝ := 1 - (8 / 5) * α
def sparseContraction (α : ℝ) : ℝ := 1 - 48 * α
def sparseWitness (α : ℝ) : Fin 3 → ℝ :=
  ![1, (1 + 24 * α) / 1024, 1 / 1024]

theorem sparseBeta_nonneg {α : ℝ} (hα : 0 ≤ α) :
    0 ≤ sparseBeta α := by unfold sparseBeta; positivity

theorem sparseBeta_le_one {α : ℝ} (hα : α ≤ 1 / 10000) :
    sparseBeta α ≤ 1 := by unfold sparseBeta; linarith

theorem sparseTilt_nonneg {α : ℝ} (hα : α ≤ 1 / 10000) :
    0 ≤ sparseTilt α := by unfold sparseTilt; linarith

theorem sparseTilt_le_one {α : ℝ} (hα : 0 ≤ α) :
    sparseTilt α ≤ 1 := by unfold sparseTilt; linarith

theorem sparseContraction_nonneg {α : ℝ} (hα : α ≤ 1 / 10000) :
    0 ≤ sparseContraction α := by unfold sparseContraction; linarith

theorem sparseWitness_dominates_ones {α : ℝ} (hα : 0 ≤ α)
    (i : Fin 3) :
    (1 : ℝ) ≤ 1024 * sparseWitness α i := by
  fin_cases i <;> simp [sparseWitness]; linarith

theorem scalable_sparse_collatz_from_numeric_gate
    (α : ℝ) (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000)
    (hgate : ∀ i : Fin 3,
      (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
        (sparseWitness α) i ≤
      sparseContraction α * sparseWitness α i)
    (n : ℕ) :
    ∑ r : State22,
      (bernoulliStepKernel (sparseBeta α) (sparseTilt α)
        (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1) ^ n) 0 r ≤
      1024 * sparseContraction α ^ n := by
  have h := scalable_bernoulli_paper_collatz_bound
    (sparseBeta α) (sparseTilt α)
    (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1)
    (sparseTilt_nonneg hα1) (sparseTilt_le_one hα0)
    (sparseWitness α) 1024 (sparseContraction α)
    (by norm_num) (sparseContraction_nonneg hα1)
    (sparseWitness_dominates_ones hα0) hgate n
  simpa [sparseWitness] using h

end Spin.UniformMixer
