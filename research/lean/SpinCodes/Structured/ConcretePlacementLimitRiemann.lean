import SpinCodes.Structured.ConcretePlacementSimplexLattice
import Mathlib.Analysis.SpecialFunctions.Choose

/-! The exact sorted-block finite sum converges to its ordered-simplex integral. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder Filter MeasureTheory Asymptotics
open scoped Topology

theorem pow_div_choose_tendsto (a : Nat) :
    Tendsto (fun n : Nat => (n : ℝ)^a / (n.choose a : ℝ)) atTop (𝓝 (a.factorial : ℝ)) := by
  have hc : ∀ᶠ n : Nat in atTop, (n.choose a : ℝ) ≠ 0 := by
    filter_upwards [eventually_ge_atTop a] with n hn
    exact_mod_cast Nat.choose_ne_zero hn
  have hh := (isEquivalent_iff_tendsto_one hc).mp (isEquivalent_choose a).symm
  have hfac : (a.factorial : ℝ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero a
  have h := hh.const_mul (a.factorial : ℝ)
  simpa only [mul_one] using h.congr (fun n => by
    change (a.factorial : ℝ) * (((n : ℝ)^a / a.factorial) / (n.choose a : ℝ)) = _
    field_simp)

theorem placement_normalizer_tendsto (a : Nat) :
    Tendsto (fun R : Nat => (128 : ℝ)^a * (R : ℝ)^a / ((128 * R).choose a : ℝ))
      atTop (𝓝 (a.factorial : ℝ)) := by
  have hn : Tendsto (fun R : Nat => 128 * R) atTop atTop :=
    tendsto_atTop_mono (fun R => by change R ≤ 128 * R; omega) tendsto_id
  simpa only [Function.comp_def, Nat.cast_mul, Nat.cast_ofNat, mul_pow] using (pow_div_choose_tendsto a).comp hn

theorem weighted_orderedSite_sum_tendsto {a : Nat} (F : (Fin a → ℝ) → ℝ) (hF : Continuous F) :
    Tendsto (fun R : Nat => ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
        ∑ B : BlockSubset R a, F (normalizedSites B))
      atTop (𝓝 ((a.factorial : ℝ) * ∫ x in orderedSiteDomain a, F x)) := by
  have h := (placement_normalizer_tendsto a).mul (orderedSite_sum_tendsto F hF)
  apply h.congr'
  filter_upwards [eventually_gt_atTop 0] with R hR
  have hRp : (R : ℝ) ≠ 0 := by exact_mod_cast Nat.ne_of_gt hR
  field_simp

def positionGaps {a : Nat} (x : Fin a → ℝ) (i : Fin (a + 1)) : ℝ :=
  (Fin.snoc x (1 : ℝ) : Fin (a + 1) → ℝ) i - (Fin.cons (0 : ℝ) x : Fin (a + 1) → ℝ) i

def siteProduct (γ : ℝ) {a : Nat} (x : Fin a → ℝ) : Matrix (Fin 2) (Fin 2) ℝ :=
  simplexProduct γ (positionGaps x)

theorem continuous_timeEmpty (γ : ℝ) : Continuous (timeEmpty γ) := by
  apply continuous_pi
  intro i
  apply continuous_pi
  intro j
  fin_cases i <;> fin_cases j <;> simp only [timeEmpty, Matrix.cons_val_zero, Matrix.cons_val_one] <;> fun_prop

theorem continuous_timeProduct {X : Type*} [TopologicalSpace X] {a : Nat} (γ : ℝ)
    (gaps : X → Fin a → ℝ) (last : X → ℝ)
    (hg : ∀ i, Continuous (fun x => gaps x i)) (hl : Continuous last) :
    Continuous (fun x => timeProduct γ (gaps x) (last x)) := by
  induction a with
  | zero => exact (continuous_timeEmpty γ).comp hl
  | succ a ih =>
    exact (((continuous_timeEmpty γ).comp (hg 0)).mul continuous_const).mul
      (ih (fun x => Fin.tail (gaps x)) (fun i => hg i.succ))

theorem continuous_positionGaps {a : Nat} (i : Fin (a + 1)) :
    Continuous (fun x : Fin a → ℝ => positionGaps x i) := by
  have hr : Continuous (fun x : Fin a → ℝ => (Fin.snoc x (1 : ℝ) : Fin (a + 1) → ℝ) i) := by
    refine Fin.lastCases ?_ (fun j => ?_) i
    · simpa only [Fin.snoc_last] using (continuous_const : Continuous (fun _ : Fin a → ℝ => (1 : ℝ)))
    · simpa only [Fin.snoc_castSucc] using (continuous_apply j : Continuous (fun x : Fin a → ℝ => x j))
  have hl : Continuous (fun x : Fin a → ℝ => (Fin.cons (0 : ℝ) x : Fin (a + 1) → ℝ) i) := by
    refine Fin.cases ?_ (fun j => ?_) i
    · simpa only [Fin.cons_zero] using (continuous_const : Continuous (fun _ : Fin a → ℝ => (0 : ℝ)))
    · simpa only [Fin.cons_succ] using (continuous_apply j : Continuous (fun x : Fin a → ℝ => x j))
  exact hr.sub hl

theorem continuous_siteProduct (γ : ℝ) (a : Nat) : Continuous (siteProduct γ (a := a)) := by
  exact continuous_timeProduct γ (fun x i => positionGaps x i.castSucc)
    (fun x => positionGaps x (Fin.last a))
    (fun i => continuous_positionGaps i.castSucc) (continuous_positionGaps (Fin.last a))

def continuumRegionKernel (θ : ℝ) (a : Nat) : Matrix (Fin 2) (Fin 2) ℝ := fun i j =>
  (a.factorial : ℝ) * ∫ x in orderedSiteDomain a, siteProduct (θ * epochMean / 128) x i j

theorem siteProduct_sum_tendsto (θ : ℝ) (a : Nat) (i j : Fin 2) :
    Tendsto (fun R : Nat => ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
        ∑ B : BlockSubset R a, siteProduct (θ * epochMean / 128) (normalizedSites B) i j)
      atTop (𝓝 (continuumRegionKernel θ a i j)) := by
  apply weighted_orderedSite_sum_tendsto
  exact (continuous_apply j).comp ((continuous_apply i).comp (continuous_siteProduct _ a))

end Spin.Structured.Placement


