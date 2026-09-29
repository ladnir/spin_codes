import SpinCodes.Structured.ConcreteFixedProfile
import Mathlib.MeasureTheory.Integral.Pi

/-! Evaluating the zero- and one-impulse continuum kernels of the actual encoder. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset ConcreteEncoder MeasureTheory
attribute [local instance] Classical.propDecidable

theorem orderedSiteDomain_one (x : Fin 1 → ℝ) :
    x ∈ orderedSiteDomain 1 ↔ x 0 ∈ Set.Ico (0:ℝ) 1 := by
  constructor
  · intro hx; exact hx.1 0
  · intro hx
    constructor
    · intro i; simpa only [Subsingleton.elim i 0, Set.mem_Ico] using hx
    · intro i j hij; exact False.elim (by have := Subsingleton.elim i j; omega)

theorem orderedSite_integral_one (F : (Fin 1 → ℝ) → ℝ) :
    (∫ x in orderedSiteDomain 1, F x) = ∫ t in (0:ℝ)..1, F (fun _ => t) := by
  rw [intervalIntegral.integral_of_le (by norm_num : (0:ℝ) ≤ 1), ← integral_Ico_eq_integral_Ioc]
  rw [← integral_indicator (orderedSiteDomain_measurable 1), ← integral_indicator measurableSet_Ico]
  have he := (volume_preserving_funUnique (Fin 1) ℝ).symm.integral_comp' ((orderedSiteDomain 1).indicator F)
  convert! he.symm.trans ?_ using 1
  apply integral_congr_ae
  filter_upwards [] with t
  have heq : (MeasurableEquiv.funUnique (Fin 1) ℝ).symm t = (fun _ : Fin 1 => t) := by
    funext i
    exact @uniqueElim_const (Fin 1) ℝ inferInstance t i
  rw [heq]
  simp only [Set.indicator, orderedSiteDomain_one]
theorem siteProduct_one (γ t : ℝ) :
    siteProduct γ (fun _ : Fin 1 => t) = timeEmpty γ t * impulseMatrix * timeEmpty γ (1-t) := by
  have h0 : positionGaps (fun _ : Fin 1 => t) (0 : Fin 2) = t := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
  have h1 : positionGaps (fun _ : Fin 1 => t) (Fin.last 1) = 1-t := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
    rfl
  change timeEmpty γ (positionGaps (fun _ : Fin 1 => t) (0 : Fin 2)) * impulseMatrix *
    timeEmpty γ (positionGaps (fun _ : Fin 1 => t) (Fin.last 1)) = _
  rw [h0, h1]
theorem continuumRegionKernel_one_integral (θ : ℝ) (i j : Fin 2) :
    continuumRegionKernel θ 1 i j =
      ∫ t in (0:ℝ)..1, (timeEmpty (θ * epochMean / 128) t * impulseMatrix *
        timeEmpty (θ * epochMean / 128) (1-t)) i j := by
  simp only [continuumRegionKernel, Nat.factorial_one, Nat.cast_one, one_mul]
  rw [orderedSite_integral_one]
  simp only [siteProduct_one]

end Spin.Structured.Placement





