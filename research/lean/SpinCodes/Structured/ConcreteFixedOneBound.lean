import SpinCodes.Structured.ConcreteFixedSmallKernelFormula

/-! A conservative rational weighted bound for the actual Q=1 continuum kernel. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder
attribute [local instance] Classical.propDecidable

theorem fugacityContinuum_one (θ : ℝ) (u : Fin 1 → ℝ) :
    fugacityContinuum θ u = continuumRegionKernel θ 0 + u 0 • continuumRegionKernel θ 1 := by
  have hs : (univ : Finset (Finset (Fin 1))) = {∅, {0}} := by decide
  simp only [fugacityContinuum, hs]
  simp [markCoefficient]

theorem exp_neg_le_half_of_one_le {γ : ℝ} (hγ : 1 ≤ γ) : Real.exp (-γ) ≤ 1/2 := by
  rw [Real.exp_neg]
  rw [one_div]
  apply (inv_le_inv₀ (Real.exp_pos γ) (by norm_num : (0:ℝ) < 2)).mpr
  have h := Real.add_one_le_exp γ
  linarith

theorem continuum_one_weighted_bound :
    Spin.RowNormLe (1/1000) (101/100) (fugacityContinuum 2 (fun _ : Fin 1 => 1)) := by
  let γ : ℝ := 2 * epochMean / 128
  have hγ : 1 ≤ γ := by norm_num [γ, epochMean]
  have hγpos : 0 < γ := by linarith
  have he0 : 0 ≤ Real.exp (-γ) := (Real.exp_pos _).le
  have he : Real.exp (-γ) ≤ 1/2 := exp_neg_le_half_of_one_le hγ
  have hf0 : 0 ≤ (1 - Real.exp (-γ)) / γ := div_nonneg (by linarith) hγpos.le
  have hf : (1 - Real.exp (-γ)) / γ ≤ 1 := by
    apply (div_le_iff₀ hγpos).mpr
    linarith
  rw [fugacityContinuum_one, continuumRegionKernel_zero, continuumRegionKernel_one hγpos.ne']
  change Spin.RowNormLe (1/1000) (101/100)
    (timeEmpty γ 1 + (1:ℝ) • !![0, (1-Real.exp (-γ))/γ;
      (1/524287)*((1-Real.exp (-γ))/γ), (1-1/524287)*Real.exp (-γ)])
  unfold Spin.RowNormLe
  simp only [one_smul, Matrix.add_apply, timeEmpty, Matrix.of_apply, Matrix.cons_val_zero, Matrix.cons_val_one, mul_one, add_zero, zero_add]
  rw [abs_of_nonneg (by norm_num : (0:ℝ) ≤ 1), abs_of_nonneg hf0,
    abs_of_nonneg (mul_nonneg (by norm_num : (0:ℝ) ≤ 1/524287) hf0),
    abs_of_nonneg (add_nonneg he0 (mul_nonneg (by norm_num : (0:ℝ) ≤ 1-1/524287) he0))]
  constructor <;> nlinarith

end Spin.Structured.Placement



