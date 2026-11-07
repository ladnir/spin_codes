import SpinCodes.UniformMixer.ScalableFixedProfile
import SpinCodes.Structured.ConcreteFixedInsertionPathBound
import SpinCodes.Structured.ConcreteFixedSmallKernelFormula

/-! The continuum regional kernel factors through the zero/live classes.
This exposes the two-state matrix used by fixed-occupation certificates. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement MeasureTheory

attribute [local instance] Classical.propDecidable

def scalableContinuumTwo (θ : ℝ) (a : ℕ) : Matrix Bool Bool ℝ :=
  fun b d => (a.factorial : ℝ) *
    ∫ x in orderedSiteDomain a,
      continuumSiteProduct (regionGamma θ) x b d

theorem scalableContinuumTwo_zero (θ : ℝ) :
    scalableContinuumTwo θ 0 = timeRegionGap (regionGamma θ) 1 := by
  ext b d
  simp only [scalableContinuumTwo, Nat.factorial_zero, Nat.cast_one,
    one_mul, orderedSiteDomain_zero, Measure.restrict_univ]
  rw [Measure.volume_pi_eq_dirac (fun i : Fin 0 => Fin.elim0 i), integral_dirac]
  simp [continuumSiteProduct, timeRegionProduct, positionGaps, Fin.snoc]

theorem continuumSiteProduct_one (γ t : ℝ) :
    continuumSiteProduct γ (fun _ : Fin 1 => t) =
      timeRegionGap γ t * coarseP * timeRegionGap γ (1 - t) := by
  have h0 : positionGaps (fun _ : Fin 1 => t) (0 : Fin 2) = t := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
  have h1 : positionGaps (fun _ : Fin 1 => t) (Fin.last 1) = 1 - t := by
    norm_num [positionGaps, Fin.snoc, Fin.cons]
    rfl
  change timeRegionGap γ (positionGaps (fun _ : Fin 1 => t) (0 : Fin 2)) *
    coarseP * timeRegionGap γ
      (positionGaps (fun _ : Fin 1 => t) (Fin.last 1)) = _
  rw [h0, h1]

theorem scalableContinuumTwo_one_integral (θ : ℝ) (b d : Bool) :
    scalableContinuumTwo θ 1 b d =
      ∫ t in (0 : ℝ)..1,
        (timeRegionGap (regionGamma θ) t * coarseP *
          timeRegionGap (regionGamma θ) (1 - t)) b d := by
  simp only [scalableContinuumTwo, Nat.factorial_one, Nat.cast_one, one_mul]
  rw [orderedSite_integral_one]
  simp only [continuumSiteProduct_one]

def scalableOneIntegrand (γ t : ℝ) : Matrix Bool Bool ℝ :=
  fun b d =>
    if b then
      if d then (1 - 1 / (liveCount State22 : ℝ)) * Real.exp (-γ)
      else (1 / (liveCount State22 : ℝ)) * Real.exp (-γ * t)
    else if d then Real.exp (-γ * (1 - t)) else 0

theorem scalableOneIntegrand_eq (γ t : ℝ) :
    timeRegionGap γ t * coarseP * timeRegionGap γ (1 - t) =
      scalableOneIntegrand γ t := by
  have he := Spin.exp_prod_one γ t
  ext b d
  cases b <;> cases d <;>
    simp only [Matrix.mul_apply, Fintype.sum_bool] <;>
    dsimp only [timeRegionGap, coarseD, coarseP,
      scalableOneIntegrand] <;>
    simp only [Bool.false_eq_true, Bool.true_eq_false,
      ite_true, ite_false]
  all_goals
    try rw [← he]
    ring

def scalableOneFormula (γ : ℝ) : Matrix Bool Bool ℝ :=
  fun b d =>
    if b then
      if d then (1 - 1 / (liveCount State22 : ℝ)) * Real.exp (-γ)
      else (1 / (liveCount State22 : ℝ)) *
        ((1 - Real.exp (-γ)) / γ)
    else if d then (1 - Real.exp (-γ)) / γ else 0

theorem scalableContinuumTwo_one {θ : ℝ}
    (hγ : regionGamma θ ≠ 0) :
    scalableContinuumTwo θ 1 = scalableOneFormula (regionGamma θ) := by
  ext b d
  rw [scalableContinuumTwo_one_integral]
  simp only [scalableOneIntegrand_eq]
  cases b <;> cases d
  · change (∫ _t in (0 : ℝ)..1, (0 : ℝ)) = 0
    simp
  · change (∫ t in (0 : ℝ)..1,
        Real.exp (-(regionGamma θ) * (1 - t))) =
        (1 - Real.exp (-(regionGamma θ))) / regionGamma θ
    exact Spin.integral_expNeg_rev hγ
  · change (∫ t in (0 : ℝ)..1,
        (1 / (liveCount State22 : ℝ)) *
          Real.exp (-(regionGamma θ) * t)) =
        (1 / (liveCount State22 : ℝ)) *
          ((1 - Real.exp (-(regionGamma θ))) / regionGamma θ)
    rw [intervalIntegral.integral_const_mul, Spin.integral_expNeg hγ]
  · change (∫ _t in (0 : ℝ)..1,
        (1 - 1 / (liveCount State22 : ℝ)) *
          Real.exp (-(regionGamma θ))) =
        (1 - 1 / (liveCount State22 : ℝ)) *
          Real.exp (-(regionGamma θ))
    simp only [intervalIntegral.integral_const, sub_zero, one_smul]

theorem continuumEndpoint_eq_lift_two (θ : ℝ) (a : ℕ)
    (q r : State22) :
    continuumEndpoint θ a q r =
      (liftM * scalableContinuumTwo θ a * projectM) q r := by
  let c : Bool := if q = 0 then false else true
  have hpoint (x : Fin a → ℝ) :
      (liftM * continuumSiteProduct (regionGamma θ) x * projectM) q r =
        ∑ d : Bool,
          continuumSiteProduct (regionGamma θ) x c d * projectM d r := by
    rw [Matrix.mul_apply]
    simp only [liftM_mul_apply]
    rfl
  have hcontinuous (d : Bool) :
      Continuous (fun x : Fin a → ℝ =>
        continuumSiteProduct (regionGamma θ) x c d) :=
    (continuous_apply d).comp
      ((continuous_apply c).comp
        (continuous_continuumSiteProduct (regionGamma θ) a))
  have hintegrable (d : Bool) :
      IntegrableOn (fun x : Fin a → ℝ =>
        continuumSiteProduct (regionGamma θ) x c d * projectM d r)
        (orderedSiteDomain a) :=
    continuous_integrable_ordered _
      ((hcontinuous d).mul continuous_const)
  calc
    continuumEndpoint θ a q r =
        (a.factorial : ℝ) *
          ∫ x in orderedSiteDomain a,
            ∑ d : Bool,
              continuumSiteProduct (regionGamma θ) x c d * projectM d r := by
      unfold continuumEndpoint
      congr 1
      apply setIntegral_congr_fun (orderedSiteDomain_measurable a)
      intro x _
      exact hpoint x
    _ = (a.factorial : ℝ) *
        ∑ d : Bool,
          ∫ x in orderedSiteDomain a,
            continuumSiteProduct (regionGamma θ) x c d * projectM d r := by
      rw [integral_finsetSum Finset.univ (fun d _ => hintegrable d)]
    _ = ∑ d : Bool,
        scalableContinuumTwo θ a c d * projectM d r := by
      rw [Finset.mul_sum]
      apply Finset.sum_congr rfl
      intro d _
      rw [integral_mul_const]
      unfold scalableContinuumTwo
      ring
    _ = (liftM * scalableContinuumTwo θ a * projectM) q r := by
      rw [Matrix.mul_apply]
      simp only [liftM_mul_apply]
      rfl

def scalableFugacityTwo {Q : ℕ} (θ : ℝ)
    (u : Fin Q → ℝ) : Matrix Bool Bool ℝ :=
  ∑ A : Finset (Fin Q),
    markCoefficient u A • scalableContinuumTwo θ A.card

theorem scalableFugacityTwo_one (θ : ℝ) (u : Fin 1 → ℝ) :
    scalableFugacityTwo θ u =
      scalableContinuumTwo θ 0 + u 0 • scalableContinuumTwo θ 1 := by
  have hs : (Finset.univ : Finset (Finset (Fin 1))) =
      {∅, {0}} := by decide
  simp only [scalableFugacityTwo, hs]
  simp [markCoefficient]

theorem scalableFugacityTwo_two (θ : ℝ) :
    scalableFugacityTwo θ (fun _ : Fin 2 => 1) =
      scalableContinuumTwo θ 0 +
        (2 : ℝ) • scalableContinuumTwo θ 1 +
        scalableContinuumTwo θ 2 := by
  have hs : (Finset.univ : Finset (Finset (Fin 2))) =
      {∅, {0}, {1}, {0, 1}} := by decide
  simp only [scalableFugacityTwo, hs]
  rw [Finset.sum_insert (by decide :
      (∅ : Finset (Fin 2)) ∉
        ({{0}, {1}, {0, 1}} : Finset (Finset (Fin 2)))),
    Finset.sum_insert (by decide :
      ({0} : Finset (Fin 2)) ∉
        ({{1}, {0, 1}} : Finset (Finset (Fin 2)))),
    Finset.sum_insert (by decide :
      ({1} : Finset (Fin 2)) ∉
        ({{0, 1}} : Finset (Finset (Fin 2)))),
    Finset.sum_singleton]
  simp [markCoefficient, two_smul]
  abel

theorem scalableFugacityContinuum_eq_lift_two {Q : ℕ}
    (θ : ℝ) (u : Fin Q → ℝ) :
    scalableFugacityContinuum θ u =
      liftM * scalableFugacityTwo θ u * projectM := by
  ext q r
  have hstep (A : Finset (Fin Q)) :
      markCoefficient u A * continuumEndpoint θ A.card q r =
        (markCoefficient u A •
          (liftM * scalableContinuumTwo θ A.card * projectM)) q r := by
    rw [continuumEndpoint_eq_lift_two]
    rfl
  simp only [scalableFugacityContinuum, hstep,
    scalableFugacityTwo]
  rw [Matrix.mul_sum, Matrix.sum_mul]
  simp only [Matrix.mul_smul, Matrix.smul_mul,
    Matrix.sum_apply, Matrix.smul_apply, smul_eq_mul]

end Spin.UniformMixer
