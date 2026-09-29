import SpinCodes.Structured.ConcreteFixedTwoKernel
import SpinCodes.Structured.ConcreteFixedKernelNonnegative

/-! Conservative rational Q=2 continuum contraction. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder
attribute [local instance] Classical.propDecidable

theorem fugacityContinuum_two (θ : ℝ) :
    fugacityContinuum θ (fun _ : Fin 2 => 1) =
      continuumRegionKernel θ 0 + (2:ℝ) • continuumRegionKernel θ 1 + continuumRegionKernel θ 2 := by
  have hs : (univ : Finset (Finset (Fin 2))) = {∅, {0}, {1}, {0,1}} := by decide
  simp only [fugacityContinuum, hs]
  rw [Finset.sum_insert (by decide : (∅ : Finset (Fin 2)) ∉ ({{0}, {1}, {0,1}} : Finset (Finset (Fin 2)))),
    Finset.sum_insert (by decide : ({0} : Finset (Fin 2)) ∉ ({{1}, {0,1}} : Finset (Finset (Fin 2)))),
    Finset.sum_insert (by decide : ({1} : Finset (Fin 2)) ∉ ({{0,1}} : Finset (Finset (Fin 2)))), Finset.sum_singleton]
  simp [markCoefficient, two_smul]
  abel

theorem exp_neg_le_quarter_of_two_le {γ : ℝ} (hγ : 2 ≤ γ) : Real.exp (-γ) ≤ 1/4 := by
  have he := exp_neg_le_half_of_one_le (show 1 ≤ γ/2 by linarith)
  have he0 := (Real.exp_pos (-(γ/2))).le
  have heq : Real.exp (-γ) = Real.exp (-(γ/2)) ^ 2 := by
    rw [pow_two, ← Real.exp_add]
    congr 1
    ring
  rw [heq]
  nlinarith

theorem continuum_two_weighted_bound :
    Spin.RowNormLe (1/1000) (101/100) (fugacityContinuum 4 (fun _ : Fin 2 => 1)) := by
  let γ : ℝ := 4 * epochMean / 128
  have hγ : 2 ≤ γ := by norm_num [γ,epochMean]
  have hγpos : 0 < γ := by linarith
  have he0 : 0 ≤ Real.exp (-γ) := (Real.exp_pos _).le
  have he : Real.exp (-γ) ≤ 1/4 := exp_neg_le_quarter_of_two_le hγ
  have hf : (1 - Real.exp (-γ)) / γ ≤ 1/2 := by
    apply (div_le_iff₀ hγpos).mpr
    linarith
  have hg : 2*(γ-1+Real.exp (-γ))/γ^2 ≤ 1 := by
    apply (div_le_iff₀ (sq_pos_of_pos hγpos)).mpr
    nlinarith [sq_nonneg (γ-1)]
  have hd : 2*(1-(1+γ)*Real.exp (-γ))/γ^2 ≤ 1/2 := by
    apply (div_le_iff₀ (sq_pos_of_pos hγpos)).mpr
    have hh : 0 ≤ (1+γ)*Real.exp (-γ) := mul_nonneg (by linarith) he0
    nlinarith [sq_nonneg (γ-2)]
  have h00 := continuumRegionKernel_two_le (show (4:ℝ)*epochMean/128 ≠ 0 from hγpos.ne') 0 0
  have h01 := continuumRegionKernel_two_le (show (4:ℝ)*epochMean/128 ≠ 0 from hγpos.ne') 0 1
  have h10 := continuumRegionKernel_two_le (show (4:ℝ)*epochMean/128 ≠ 0 from hγpos.ne') 1 0
  have h11 := continuumRegionKernel_two_le (show (4:ℝ)*epochMean/128 ≠ 0 from hγpos.ne') 1 1
  change continuumRegionKernel 4 2 0 0 ≤ (1/524287) * (2*(γ-1+Real.exp (-γ))/γ^2) at h00
  change continuumRegionKernel 4 2 0 1 ≤ 2*(1-(1+γ)*Real.exp (-γ))/γ^2 at h01
  change continuumRegionKernel 4 2 1 0 ≤ (1/524287)*(2*(1-(1+γ)*Real.exp (-γ))/γ^2) at h10
  change continuumRegionKernel 4 2 1 1 ≤ (1/524287)*(2*(1-(1+γ)*Real.exp (-γ))/γ^2)+Real.exp (-γ) at h11
  unfold Spin.RowNormLe
  simp only [abs_of_nonneg (fugacityContinuum_nonneg 4 (fun _ : Fin 2 => 1) (by intro i; norm_num) _ _)]
  rw [fugacityContinuum_two, continuumRegionKernel_zero, continuumRegionKernel_one hγpos.ne']
  change (timeEmpty γ 1 + (2:ℝ) • !![0,(1-Real.exp (-γ))/γ;
    (1/524287)*((1-Real.exp (-γ))/γ),(1-1/524287)*Real.exp (-γ)] + continuumRegionKernel 4 2) 0 0 +
      (timeEmpty γ 1 + (2:ℝ) • !![0,(1-Real.exp (-γ))/γ;
    (1/524287)*((1-Real.exp (-γ))/γ),(1-1/524287)*Real.exp (-γ)] + continuumRegionKernel 4 2) 0 1 * (1/1000) ≤ 101/100 ∧ _
  simp only [Matrix.add_apply, Matrix.smul_apply, smul_eq_mul, timeEmpty, Matrix.of_apply,
    Matrix.cons_val_zero, Matrix.cons_val_one, mul_one, add_zero, zero_add, mul_zero]
  constructor
  · nlinarith
  · change 2*((1/524287)*((1-Real.exp (-γ))/γ)) + continuumRegionKernel 4 2 1 0 +
      (Real.exp (-γ)+2*((1-1/524287)*Real.exp (-γ))+continuumRegionKernel 4 2 1 1)*(1/1000) ≤ (101/100)*(1/1000)
    nlinarith

end Spin.Structured.Placement


