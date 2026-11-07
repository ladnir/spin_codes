import SpinCodes.UniformMixer.ScalableFixedOneNumeric
import SpinCodes.Structured.ConcreteFixedTwoKernel
import Mathlib.LinearAlgebra.Matrix.Reindex

/-! Compare the actual two-site Scalable kernel with the integrated
two-site envelope through the zero/live index equivalence. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement MeasureTheory

def toFinTwo : Matrix Bool Bool ℝ → Matrix (Fin 2) (Fin 2) ℝ :=
  Matrix.reindexAlgEquiv ℝ ℝ finTwoEquiv.symm

theorem toFinTwo_mul (A B : Matrix Bool Bool ℝ) :
    toFinTwo (A * B) = toFinTwo A * toFinTwo B :=
  map_mul (Matrix.reindexAlgEquiv ℝ ℝ finTwoEquiv.symm) A B

theorem toFinTwo_apply (A : Matrix Bool Bool ℝ)
    (i j : Fin 2) :
    toFinTwo A i j = A (finTwoEquiv i) (finTwoEquiv j) := by
  rfl

theorem toFinTwo_gap (γ u : ℝ) :
    toFinTwo (timeRegionGap γ u) = timeEmpty γ u := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    norm_num [toFinTwo_apply, finTwoEquiv,
      timeRegionGap, coarseD, timeEmpty]

theorem toFinTwo_impulse :
    toFinTwo coarseP =
      !![0, 1; 1 / (liveCount State22 : ℝ),
        1 - 1 / (liveCount State22 : ℝ)] := by
  ext i j
  fin_cases i <;> fin_cases j <;> rfl

theorem finMatrix_mul_nonneg
    (A B : Matrix (Fin 2) (Fin 2) ℝ)
    (hA : ∀ i j, 0 ≤ A i j) (hB : ∀ i j, 0 ≤ B i j) :
    ∀ i j, 0 ≤ (A * B) i j := by
  intro i j
  rw [Matrix.mul_apply]
  exact Finset.sum_nonneg (fun k _ => mul_nonneg (hA i k) (hB k j))

theorem finMatrix_mul_mono
    (A A' B B' : Matrix (Fin 2) (Fin 2) ℝ)
    (hA : ∀ i j, A i j ≤ A' i j)
    (hB : ∀ i j, B i j ≤ B' i j)
    (hA' : ∀ i j, 0 ≤ A' i j)
    (hB0 : ∀ i j, 0 ≤ B i j) :
    ∀ i j, (A * B) i j ≤ (A' * B') i j := by
  intro i j
  simp only [Matrix.mul_apply]
  apply Finset.sum_le_sum
  intro k _
  exact mul_le_mul (hA i k) (hB k j) (hB0 k j) (hA' i k)

theorem scalable_impulse_nonneg :
    ∀ i j, 0 ≤ toFinTwo coarseP i j := by
  rw [toFinTwo_impulse]
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [state22_liveCount]

theorem plus_impulse_nonneg :
    ∀ i j, 0 ≤ Spin.Pplus (1 / 524287) i j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [Spin.Pplus]

theorem scalable_impulse_le_plus :
    ∀ i j, toFinTwo coarseP i j ≤
      Spin.Pplus (1 / 524287) i j := by
  rw [toFinTwo_impulse]
  intro i j
  fin_cases i <;> fin_cases j <;>
    norm_num [state22_liveCount, Spin.Pplus]

private theorem timeEmpty_nonneg (γ u : ℝ) :
    ∀ i j, 0 ≤ timeEmpty γ u i j := by
  intro i j
  fin_cases i <;> fin_cases j <;>
    simp [timeEmpty]; positivity

theorem continuumSiteProduct_two (γ : ℝ) (x : Fin 2 → ℝ) :
    continuumSiteProduct γ x =
      timeRegionGap γ (x 0) * coarseP *
        timeRegionGap γ (x 1 - x 0) * coarseP *
        timeRegionGap γ (1 - x 1) := by
  have h0 : positionGaps x (0 : Fin 3) = x 0 := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
  have h1 : positionGaps x (1 : Fin 3) = x 1 - x 0 := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
    rfl
  have h2 : positionGaps x (Fin.last 2) = 1 - x 1 := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
    rfl
  change timeRegionGap γ (positionGaps x (0 : Fin 3)) * coarseP *
    (timeRegionGap γ (positionGaps x (1 : Fin 3)) * coarseP *
      timeRegionGap γ (positionGaps x (Fin.last 2))) = _
  rw [h0, h1, h2]
  simp only [Matrix.mul_assoc]

private theorem scalable_site_two_le_plus (γ : ℝ) (x : Fin 2 → ℝ)
    (i j : Fin 2) :
    toFinTwo (continuumSiteProduct γ x) i j ≤
      plusSiteTwo γ x i j := by
  let D0 := timeEmpty γ (x 0)
  let D1 := timeEmpty γ (x 1 - x 0)
  let D2 := timeEmpty γ (1 - x 1)
  let P := toFinTwo coarseP
  let L := Spin.Pplus (1 / 524287)
  have hD0 : ∀ i j, 0 ≤ D0 i j := timeEmpty_nonneg _ _
  have hD1 : ∀ i j, 0 ≤ D1 i j := timeEmpty_nonneg _ _
  have hD2 : ∀ i j, 0 ≤ D2 i j := timeEmpty_nonneg _ _
  have hP : ∀ i j, 0 ≤ P i j := scalable_impulse_nonneg
  have hL : ∀ i j, 0 ≤ L i j := plus_impulse_nonneg
  have hPL : ∀ i j, P i j ≤ L i j := scalable_impulse_le_plus
  have h1 : ∀ i j, 0 ≤ (D0 * L) i j :=
    finMatrix_mul_nonneg D0 L hD0 hL
  have h2 : ∀ i j, 0 ≤ (D0 * L * D1) i j :=
    finMatrix_mul_nonneg (D0 * L) D1 h1 hD1
  have h3 : ∀ i j, 0 ≤ (D0 * L * D1 * L) i j :=
    finMatrix_mul_nonneg (D0 * L * D1) L h2 hL
  have h01 : ∀ i j, (D0 * P) i j ≤ (D0 * L) i j :=
    finMatrix_mul_mono D0 D0 P L (fun _ _ => le_refl _) hPL hD0 hP
  have h012 : ∀ i j, (D0 * P * D1) i j ≤
      (D0 * L * D1) i j :=
    finMatrix_mul_mono (D0 * P) (D0 * L) D1 D1 h01
      (fun _ _ => le_refl _) h1 hD1
  have h0123 : ∀ i j, (D0 * P * D1 * P) i j ≤
      (D0 * L * D1 * L) i j :=
    finMatrix_mul_mono (D0 * P * D1) (D0 * L * D1) P L
      h012 hPL h2 hP
  have h : (D0 * P * D1 * P * D2) i j ≤
      (D0 * L * D1 * L * D2) i j :=
    finMatrix_mul_mono (D0 * P * D1 * P)
      (D0 * L * D1 * L) D2 D2 h0123
      (fun _ _ => le_refl _) h3 hD2 i j
  simpa only [D0, D1, D2, P, L, plusSiteTwo,
    continuumSiteProduct_two, toFinTwo_mul,
    toFinTwo_gap] using h

theorem scalableContinuumTwo_two_le {θ : ℝ}
    (hγ : regionGamma θ ≠ 0) (i j : Fin 2) :
    scalableContinuumTwo θ 2 (finTwoEquiv i) (finTwoEquiv j) ≤
      enlargedKernelTwo (regionGamma θ) i j := by
  rw [← plusSiteTwo_integral hγ i j]
  unfold scalableContinuumTwo
  norm_num only [Nat.factorial_succ, Nat.factorial_zero,
    Nat.reduceAdd, Nat.reduceMul, Nat.cast_ofNat]
  apply mul_le_mul_of_nonneg_left _ (by norm_num)
  have hi (F : (Fin 2 → ℝ) → ℝ) (hF : Continuous F) :
      IntegrableOn F (orderedSiteDomain 2) := by
    apply hF.integrableOn_Icc.mono_set
    intro x hx
    exact ⟨fun k => (hx.1 k).1, fun k => (hx.1 k).2.le⟩
  apply setIntegral_mono_on
    (hi _ ((continuous_apply (finTwoEquiv j)).comp
      ((continuous_apply (finTwoEquiv i)).comp
        (continuous_continuumSiteProduct _ 2))))
    (hi _ ((continuous_apply j).comp
      ((continuous_apply i).comp (continuous_plusSiteTwo _))))
    (orderedSiteDomain_measurable 2)
  intro x _
  change continuumSiteProduct (regionGamma θ) x
      (finTwoEquiv i) (finTwoEquiv j) ≤ plusSiteTwo (regionGamma θ) x i j
  simpa only [toFinTwo_apply] using
    scalable_site_two_le_plus (regionGamma θ) x i j

end Spin.UniformMixer
