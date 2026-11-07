import SpinCodes.UniformMixer.ScalableSparseNumericGate
import SpinCodes.UniformMixer.ScalableFixedFairSmall
import SpinCodes.Structured.ConcreteMarkedMoment

/-! Conditioned fair-row sparse bound for the actual sampled Scalable inner. -/

noncomputable section
set_option maxRecDepth 8192
set_option maxHeartbeats 1000000
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteMarked
  Spin.Structured.ConcreteRoute Spin.Structured.ConcreteNativeFamily

private theorem bernoulliKernel_congr {β β' z : ℝ}
    (h : β = β') (h0 : 0 ≤ β) (h1 : β ≤ 1)
    (h0' : 0 ≤ β') (h1' : β' ≤ 1) :
    bernoulliStepKernel β z h0 h1 =
      bernoulliStepKernel β' z h0' h1' := by
  subst β'
  rfl

theorem scalable_fair_sparse_moment (m : ℕ)
    (S : Finset (Fin (Lsched m))) {α : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) :
    (regionsLaw S (bsched m)).expect
      (scalableRoutedMoment m (sparseTilt α)) ≤
      (1 / markMass (Lsched m) S.card ((8 / 5) * α)) ^
        (bsched m) *
        (1024 * sparseContraction α ^
          (bsched m * scalableStepsPerRegion m)) := by
  let p : ℝ := (8 / 5) * α
  have hp0 : 0 < p := by dsimp [p]; positivity
  have hp1 : p < 1 := by dsimp [p]; linarith
  have hz : 0 ≤ sparseTilt α := sparseTilt_nonneg hα1
  have hm := markMass_pos S hp0 hp1
  have hd := regionsLaw_dominates S (bsched m) p hp0.le hp1.le hm
  have href :
      (iidLaw (Lsched m) (bsched m) (p / 2)
        (by positivity) (by linarith)).expect
        (scalableRoutedMoment m (sparseTilt α)) =
      ∑ r : State22,
        ((bernoulliStepKernel (sparseBeta α) (sparseTilt α)
          (sparseBeta_nonneg hα0.le)
          (sparseBeta_le_one hα1)) ^
            (bsched m * scalableStepsPerRegion m)) 0 r := by
    have hb : p / 2 = sparseBeta α := by
      dsimp [p, sparseBeta]
      ring
    have hpower := scalable_iid_moment_power m
      (p / 2) (sparseTilt α) (by positivity) (by linarith)
    have hk : bernoulliStepKernel (p / 2) (sparseTilt α)
        (by positivity) (by linarith) =
        bernoulliStepKernel (sparseBeta α) (sparseTilt α)
          (sparseBeta_nonneg hα0.le)
          (sparseBeta_le_one hα1) :=
      bernoulliKernel_congr hb _ _ _ _
    rw [hk] at hpower
    exact hpower
  calc
    _ ≤ (1 / markMass (Lsched m) S.card p) ^ (bsched m) *
        (iidLaw (Lsched m) (bsched m) (p / 2)
          (by positivity) (by linarith)).expect
          (scalableRoutedMoment m (sparseTilt α)) := by
      exact hd.expect_le (fun regions =>
        scalableRoutedMoment_nonneg m hz regions)
    _ = (1 / markMass (Lsched m) S.card p) ^ (bsched m) *
        (∑ r : State22,
          ((bernoulliStepKernel (sparseBeta α) (sparseTilt α)
            (sparseBeta_nonneg hα0.le)
            (sparseBeta_le_one hα1)) ^
              (bsched m * scalableStepsPerRegion m)) 0 r) := by
      rw [href]
    _ ≤ _ := by
      apply mul_le_mul_of_nonneg_left
      · exact scalable_sparse_collatz_bound α hα0.le hα1 _
      · positivity

theorem scalable_fair_sparse_probability (m : ℕ)
    (S : Finset (Fin (Lsched m))) {α : ℝ}
    (hα0 : 0 < α) (hα1 : α ≤ 1 / 10000) (d : ℕ) :
    ((regionsLaw S (bsched m)).prod
      (scalableStateMapLaw m)).prob
      (fun ω => sampledOutputWeight (scalableFlatInput m ω.1)
        (scalableFlatMaps m ω.2) 0 ≤ d) ≤
      (((1 / markMass (Lsched m) S.card ((8 / 5) * α)) ^
          (bsched m) *
          (1024 * sparseContraction α ^
            (bsched m * scalableStepsPerRegion m))) /
        (sparseTilt α) ^ d) := by
  apply FinPMF.prob_weight_le_of_moment _ _ d
    (by unfold sparseTilt; linarith)
    (by unfold sparseTilt; linarith)
  rw [FinPMF.expect_prod]
  simp_rw [scalable_fixed_route_moment]
  exact scalable_fair_sparse_moment m S hα0 hα1

#print axioms scalable_fair_sparse_moment
#print axioms scalable_fair_sparse_probability
end Spin.UniformMixer
