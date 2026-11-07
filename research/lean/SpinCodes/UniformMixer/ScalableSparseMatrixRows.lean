import SpinCodes.UniformMixer.ScalableSparseScalarBounds

/-! The three rows of the paper-oriented sparse transfer. -/

noncomputable section

namespace Spin.UniformMixer

theorem sparse_paper_row_zero (α : ℝ) :
    (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
      (sparseWitness α) 0 =
      spectrumH (sparseBeta α) (sparseTilt α) +
        spectrumG (sparseBeta α) (sparseTilt α) *
          ((1 + 24 * α) / 1024) := by
  simp only [paperSpectrumTransfer, spectrumTransfer, Matrix.mulVec,
    dotProduct, Fin.sum_univ_three, Matrix.transpose_apply]
  dsimp only [weightedTransfer, sparseWitness]
  norm_num

theorem sparse_paper_row_fresh (α : ℝ) :
    (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
      (sparseWitness α) 1 =
      min (spectrumF (sparseBeta α) (sparseTilt α))
          (spectrumu (sparseBeta α)) /
          (liveCount State22 : ℝ) +
        spectrumF (sparseBeta α) (sparseTilt α) / 1024 := by
  simp only [paperSpectrumTransfer, spectrumTransfer, Matrix.mulVec,
    dotProduct, Fin.sum_univ_three, Matrix.transpose_apply]
  dsimp only [weightedTransfer, sparseWitness]
  change
    (min (spectrumF (sparseBeta α) (sparseTilt α))
        (spectrumu (sparseBeta α)) / (liveCount State22 : ℝ)) * 1 +
      (0 : ℝ) * ((1 + 24 * α) / 1024) +
      spectrumF (sparseBeta α) (sparseTilt α) * (1 / 1024) = _
  ring

theorem sparse_paper_row_uniform (α : ℝ) :
    (paperSpectrumTransfer (sparseBeta α) (sparseTilt α)).mulVec
      (sparseWitness α) 2 =
      min (spectrumU (sparseBeta α) (sparseTilt α))
          (spectrumu (sparseBeta α)) /
          (liveCount State22 : ℝ) +
        spectrumU (sparseBeta α) (sparseTilt α) / 1024 := by
  simp only [paperSpectrumTransfer, spectrumTransfer, Matrix.mulVec,
    dotProduct, Fin.sum_univ_three, Matrix.transpose_apply]
  dsimp only [weightedTransfer, sparseWitness]
  change
    (min (spectrumU (sparseBeta α) (sparseTilt α))
        (spectrumu (sparseBeta α)) / (liveCount State22 : ℝ)) * 1 +
      (0 : ℝ) * ((1 + 24 * α) / 1024) +
      spectrumU (sparseBeta α) (sparseTilt α) * (1 / 1024) = _
  ring

end Spin.UniformMixer
