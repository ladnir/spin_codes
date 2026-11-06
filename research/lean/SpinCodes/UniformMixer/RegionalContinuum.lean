import SpinCodes.UniformMixer.RegionalGridLimit

/-!
# Continuum two-state regional product

The gap factors are diagonal, and each marked step applies the same
zero/live impulse matrix. This is the integrand of the paper's
ordered-simplex regional kernel. We retain the full-state boundary maps.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement Filter MeasureTheory
open scoped Topology

def timeRegionGap (γ u : ℝ) : Matrix Bool Bool ℝ :=
  coarseD (Real.exp (-γ * u))

def timeRegionProduct (γ : ℝ) :
    {a : ℕ} → (Fin a → ℝ) → ℝ → Matrix Bool Bool ℝ
  | 0, _, last => timeRegionGap γ last
  | _ + 1, gaps, last =>
      timeRegionGap γ (gaps 0) * coarseP *
        timeRegionProduct γ (Fin.tail gaps) last

theorem continuous_timeRegionGap (γ : ℝ) :
    Continuous (timeRegionGap γ) := by
  apply continuous_pi
  intro b
  apply continuous_pi
  intro d
  cases b <;> cases d <;>
    simp [timeRegionGap, coarseD] <;> fun_prop

theorem continuous_timeRegionProduct {X : Type*} [TopologicalSpace X]
    {a : ℕ} (γ : ℝ) (gaps : X → Fin a → ℝ) (last : X → ℝ)
    (hg : ∀ i, Continuous (fun x => gaps x i))
    (hl : Continuous last) :
    Continuous (fun x => timeRegionProduct γ (gaps x) (last x)) := by
  induction a with
  | zero => exact (continuous_timeRegionGap γ).comp hl
  | succ a ih =>
      exact (((continuous_timeRegionGap γ).comp (hg 0)).mul
        continuous_const).mul
        (ih (fun x => Fin.tail (gaps x))
          (fun i => hg i.succ))

def continuumSiteProduct (γ : ℝ) {a : ℕ}
    (x : Fin a → ℝ) : Matrix Bool Bool ℝ :=
  timeRegionProduct γ (fun i => positionGaps x i.castSucc)
    (positionGaps x (Fin.last a))

theorem continuous_continuumSiteProduct (γ : ℝ) (a : ℕ) :
    Continuous (continuumSiteProduct γ (a := a)) := by
  exact continuous_timeRegionProduct γ
    (fun x i => positionGaps x i.castSucc)
    (fun x => positionGaps x (Fin.last a))
    (fun i => continuous_positionGaps i.castSucc)
    (continuous_positionGaps (Fin.last a))

def regionGamma (θ : ℝ) : ℝ :=
  θ * liveWeightAverage coordMap / 64

def continuumEndpoint (θ : ℝ) (a : ℕ) (q r : State22) : ℝ :=
  (a.factorial : ℝ) * ∫ x in orderedSiteDomain a,
    (liftM * continuumSiteProduct (regionGamma θ) x * projectM) q r

theorem continuous_continuumEndpointIntegrand (θ : ℝ)
    (a : ℕ) (q r : State22) :
    Continuous (fun x : Fin a → ℝ =>
      (liftM * continuumSiteProduct (regionGamma θ) x * projectM) q r) := by
  have hM : Continuous (fun x : Fin a → ℝ =>
      liftM * continuumSiteProduct (regionGamma θ) x * projectM) :=
    ((continuous_const : Continuous (fun _ : Fin a → ℝ => liftM)).matrix_mul
      (continuous_continuumSiteProduct (regionGamma θ) a)).matrix_mul
      (continuous_const : Continuous (fun _ : Fin a → ℝ => projectM))
  exact (continuous_apply r).comp ((continuous_apply q).comp hM)

theorem idealRegionalGrid_tendsto (θ : ℝ) (a : ℕ)
    (q r : State22) :
    Tendsto (fun R : ℕ =>
      (∑ B : BlockSubset (64 * R) a,
        (liftM * continuumSiteProduct (regionGamma θ)
          (normalizedSites B) * projectM) q r) /
        (((64 * R).choose a : ℕ) : ℝ)) atTop
      (𝓝 (continuumEndpoint θ a q r)) := by
  exact regional_grid_limit _
    (continuous_continuumEndpointIntegrand θ a q r)

end Spin.UniformMixer
