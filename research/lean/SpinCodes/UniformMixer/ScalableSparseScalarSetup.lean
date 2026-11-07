import SpinCodes.UniformMixer.ScalableSparseSpectrumMeanUpper
import SpinCodes.UniformMixer.ScalableSparseSpectrumHoldUpper
import SpinCodes.UniformMixer.ScalableSparseFeedbackUpper

/-! Rational inequalities for the sparse interval. -/

noncomputable section

namespace Spin.UniformMixer

def sparseLoss (α : ℝ) : ℝ :=
  (1 - sparseBeta α) * (1 - sparseTilt α)

theorem sparseLoss_formula (α : ℝ) :
    sparseLoss α = (8 / 5) * α - (32 / 25) * α ^ 2 := by
  unfold sparseLoss sparseBeta sparseTilt
  ring

theorem sparse_square_le {α : ℝ}
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    α ^ 2 ≤ α / 10000 := by
  nlinarith [mul_le_mul_of_nonneg_left hα1 hα0]

theorem sparseLoss_nonneg {α : ℝ}
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    0 ≤ sparseLoss α := by
  unfold sparseLoss
  exact mul_nonneg (by unfold sparseBeta; linarith)
    (by unfold sparseTilt; linarith)

theorem sparseLoss_le_linear {α : ℝ} (_hα0 : 0 ≤ α) :
    sparseLoss α ≤ (8 / 5) * α := by
  rw [sparseLoss_formula]
  nlinarith [sq_nonneg α]

theorem sparseLoss_square_le {α : ℝ}
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    (sparseLoss α) ^ 2 ≤ (64 / 25) * α ^ 2 := by
  have h := pow_le_pow_left₀ (sparseLoss_nonneg hα0 hα1)
    (sparseLoss_le_linear hα0) 2
  nlinarith only [h]

theorem sparse_tail_le {α : ℝ}
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    (2 ^ 64 : ℝ) * (sparseBeta α) ^ 8 ≤ α / 1000 := by
  have hβ : sparseBeta α ≤ α := by unfold sparseBeta; linarith
  have hpow : (sparseBeta α) ^ 8 ≤ α ^ 8 :=
    pow_le_pow_left₀ (sparseBeta_nonneg hα0) hβ _
  have hα7 : α ^ 7 ≤ (1 / 10000 : ℝ) ^ 7 :=
    pow_le_pow_left₀ hα0 hα1 _
  have hpow8 : α ^ 8 ≤ α * (1 / 10000 : ℝ) ^ 7 := by
    calc
      α ^ 8 = α * α ^ 7 := by ring
      _ ≤ α * (1 / 10000 : ℝ) ^ 7 :=
        mul_le_mul_of_nonneg_left hα7 hα0
  have hconst : (2 ^ 64 : ℝ) * (1 / 10000 : ℝ) ^ 7 ≤ 1 / 1000 := by
    norm_num
  calc
    (2 ^ 64 : ℝ) * (sparseBeta α) ^ 8 ≤
        (2 ^ 64 : ℝ) * α ^ 8 :=
      mul_le_mul_of_nonneg_left hpow (by positivity)
    _ ≤ (2 ^ 64 : ℝ) * (α * (1 / 10000 : ℝ) ^ 7) :=
      mul_le_mul_of_nonneg_left hpow8 (by positivity)
    _ = α * ((2 ^ 64 : ℝ) * (1 / 10000 : ℝ) ^ 7) := by ring
    _ ≤ α * (1 / 1000) := mul_le_mul_of_nonneg_left hconst hα0
    _ = α / 1000 := by ring

end Spin.UniformMixer
