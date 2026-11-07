import SpinCodes.UniformMixer.ScalableSparseMatrixRows

/-! The sparse Collatz inequalities on the rational interval. -/

noncomputable section

namespace Spin.UniformMixer

private theorem sparse_small_divisor (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    min (spectrumF (sparseBeta α) (sparseTilt α))
        (spectrumu (sparseBeta α)) / (liveCount State22 : ℝ) ≤
      α / 4096 := by
  have hu := sparse_spectrumu_le α hα0 hα1
  have hm : min (spectrumF (sparseBeta α) (sparseTilt α))
      (spectrumu (sparseBeta α)) ≤ 52 * α :=
    (min_le_right _ _).trans hu
  have hD : (liveCount State22 : ℝ) = 4194303 := by
    norm_num [state22_liveCount]
  rw [hD]
  apply (div_le_iff₀ (by norm_num : (0 : ℝ) < 4194303)).2
  nlinarith only [hm, hα0]

private theorem sparse_small_divisor_uniform (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    min (spectrumU (sparseBeta α) (sparseTilt α))
        (spectrumu (sparseBeta α)) / (liveCount State22 : ℝ) ≤
      α / 4096 := by
  have hu := sparse_spectrumu_le α hα0 hα1
  have hm : min (spectrumU (sparseBeta α) (sparseTilt α))
      (spectrumu (sparseBeta α)) ≤ 52 * α :=
    (min_le_right _ _).trans hu
  have hD : (liveCount State22 : ℝ) = 4194303 := by
    norm_num [state22_liveCount]
  rw [hD]
  apply (div_le_iff₀ (by norm_num : (0 : ℝ) < 4194303)).2
  nlinarith only [hm, hα0]

theorem sparse_paper_row_zero_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
      (sparseWitness α) 0 ≤
      sparseContraction α * sparseWitness α 0 := by
  rw [sparse_paper_row_zero]
  have hH := sparse_spectrumH_le α hα0 hα1
  have hG := sparse_spectrumG_le α hα0 hα1
  have hv : 0 ≤ (1 + 24 * α) / 1024 := by linarith
  have hGv := mul_le_mul_of_nonneg_right hG hv
  have hsq := sparse_square_le hα0 hα1
  simp only [sparseContraction, sparseWitness, Matrix.cons_val_zero,
    mul_one]
  nlinarith only [hH, hGv, hsq]

theorem sparse_paper_row_fresh_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
      (sparseWitness α) 1 ≤
      sparseContraction α * sparseWitness α 1 := by
  rw [sparse_paper_row_fresh]
  have hmin := sparse_small_divisor α hα0 hα1
  have hF := sparse_spectrumF_le α hα0 hα1
  have hFdiv : spectrumF (sparseBeta α) (sparseTilt α) / 1024 ≤
      (1 - 25 * α) / 1024 :=
    div_le_div_of_nonneg_right hF (by norm_num)
  have hsq := sparse_square_le hα0 hα1
  change _ ≤ (1 - 48 * α) * ((1 + 24 * α) / 1024)
  nlinarith only [hmin, hFdiv, hsq]

theorem sparse_paper_row_uniform_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
      (sparseWitness α) 2 ≤
      sparseContraction α * sparseWitness α 2 := by
  rw [sparse_paper_row_uniform]
  have hmin := sparse_small_divisor_uniform α hα0 hα1
  have hU := sparse_spectrumU_le α hα0 hα1
  have hUdiv : spectrumU (sparseBeta α) (sparseTilt α) / 1024 ≤
      (1 - 50 * α) / 1024 :=
    div_le_div_of_nonneg_right hU (by norm_num)
  change _ ≤ (1 - 48 * α) * (1 / 1024)
  linarith only [hmin, hUdiv, hα0]

theorem scalable_sparse_numeric_gate (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    ∀ i : Fin 3,
      (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
        (sparseWitness α) i ≤
        sparseContraction α * sparseWitness α i := by
  intro i
  fin_cases i
  · exact sparse_paper_row_zero_le α hα0 hα1
  · exact sparse_paper_row_fresh_le α hα0 hα1
  · exact sparse_paper_row_uniform_le α hα0 hα1

theorem scalable_sparse_collatz_bound
    (α : ℝ) (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000)
    (n : ℕ) :
    ∑ r : State22,
      (bernoulliStepKernel (sparseBeta α) (sparseTilt α)
        (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1) ^ n) 0 r ≤
      1024 * sparseContraction α ^ n :=
  scalable_sparse_collatz_from_numeric_gate α hα0 hα1
    (scalable_sparse_numeric_gate α hα0 hα1) n

end Spin.UniformMixer
