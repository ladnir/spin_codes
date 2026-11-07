import SpinCodes.UniformMixer.ScalableSparseScalarSetup

/-! Coarse first-order bounds with ample slack for the sparse witness. -/

noncomputable section

namespace Spin.UniformMixer

theorem sparse_spectrumH_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    spectrumH (sparseBeta α) (sparseTilt α) ≤ 1 - 50 * α := by
  have hH := spectrumH_sparse_quadratic
    (sparseBeta α) (sparseTilt α)
    (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1)
    (sparseTilt_nonneg hα1) (sparseTilt_le_one hα0)
  have hsq := sparse_square_le hα0 hα1
  have htail := sparse_tail_le hα0 hα1
  have hb2 : (sparseBeta α) ^ 2 = (16 / 25) * α ^ 2 := by
    unfold sparseBeta
    ring
  have hpoly :
      1 - 64 * sparseBeta α + 2016 * (sparseBeta α) ^ 2 +
          (2 ^ 64 : ℝ) * (sparseBeta α) ^ 8 ≤ 1 - 50 * α := by
    rw [hb2]
    unfold sparseBeta
    change (2 ^ 64 : ℝ) * ((4 / 5) * α) ^ 8 ≤ α / 1000 at htail
    nlinarith only [hsq, htail]
  exact hH.trans hpoly

theorem sparse_spectrumG_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    spectrumG (sparseBeta α) (sparseTilt α) ≤ 52 * α := by
  have h := spectrumG_sparse_upper
    (sparseBeta α) (sparseTilt α)
    (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1)
    (sparseTilt_nonneg hα1) (sparseTilt_le_one hα0)
  have hnum : 64 * sparseBeta α ≤ 52 * α := by
    unfold sparseBeta
    linarith
  exact h.trans hnum

theorem sparse_spectrumu_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    spectrumu (sparseBeta α) ≤ 52 * α := by
  have h := spectrumu_sparse_upper
    (sparseBeta α) (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1)
  have hnum : 64 * sparseBeta α ≤ 52 * α := by
    unfold sparseBeta
    linarith
  exact h.trans hnum

theorem sparse_spectrumF_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    spectrumF (sparseBeta α) (sparseTilt α) ≤ 1 - 25 * α := by
  have hF := spectrumF_sparse_upper
    (sparseBeta α) (sparseTilt α)
    (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1)
    (sparseTilt_nonneg hα1) (sparseTilt_le_one hα0)
  change spectrumF (sparseBeta α) (sparseTilt α) ≤
    1 - 16 * sparseLoss α + 120 * (sparseLoss α) ^ 2 at hF
  have hsq := sparse_square_le hα0 hα1
  have hLoss : (8 / 5) * α - (32 / 25) * α ^ 2 ≤ sparseLoss α :=
    le_of_eq (sparseLoss_formula α).symm
  have hLossSq := sparseLoss_square_le hα0 hα1
  have hpoly :
      1 - 16 * sparseLoss α + 120 * (sparseLoss α) ^ 2 ≤
        1 - 25 * α := by
    nlinarith only [hsq, hLoss, hLossSq]
  exact hF.trans hpoly

theorem sparse_spectrumU_le (α : ℝ)
    (hα0 : 0 ≤ α) (hα1 : α ≤ 1 / 10000) :
    spectrumU (sparseBeta α) (sparseTilt α) ≤ 1 - 50 * α := by
  have hU := spectrumU_sparse_upper
    (sparseBeta α) (sparseTilt α)
    (sparseBeta_nonneg hα0) (sparseBeta_le_one hα1)
    (sparseTilt_nonneg hα1) (sparseTilt_le_one hα0)
  change spectrumU (sparseBeta α) (sparseTilt α) ≤
    1 - 32 * sparseLoss α + 2016 * (sparseLoss α) ^ 2 at hU
  have hsq := sparse_square_le hα0 hα1
  have hLoss : (8 / 5) * α - (32 / 25) * α ^ 2 ≤ sparseLoss α :=
    le_of_eq (sparseLoss_formula α).symm
  have hLossSq := sparseLoss_square_le hα0 hα1
  have hpoly :
      1 - 32 * sparseLoss α + 2016 * (sparseLoss α) ^ 2 ≤
        1 - 50 * α := by
    nlinarith only [hsq, hLoss, hLossSq]
  exact hU.trans hpoly

end Spin.UniformMixer
