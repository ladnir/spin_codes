import SpinCodes.Turbo.RoutingPairedPhysical22

/-!
# Actual fixed-active-set fractional endpoint

Historical audit: `ConcreteFixedCoefficient` preserves ordered matrix factors
under finite averaging; `ScalableEnvelopeIterate` preserves the boundary state
until its terminal sum. Here the original routing law is disintegrated before
capping and taking fractional powers. Numerical row soundness is a premise.
-/
noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Turbo
open Finset Spin.UniformMixer

theorem turboActiveBadCount22_prob_le_paired_dp
    (n : ℕ) (points : Fin 8 → Field16) (hpoints : Function.Injective points)
    (S : Finset (Fin ((n * 2) * 16))) (d : ℕ)
    (z : ℝ) (hz0 : 0 < z) (hz1 : z ≤ 1)
    (B : Fin 17 → Matrix (Fin 11) (Fin 11) ℝ)
    (hs : ∀ j, physicalRowSound22 (Matrix.of (turboSelectedKernel22 j z)) (B j))
    (hB : ∀ j a b, 0 ≤ B j a b)
    {α : ℝ} (hα0 : 0 < α) (hα1 : α ≤ 1) :
    (turboSetup22 (n * 2) points).joint.prob
      (fun ω => 1 ≤ turboActiveBadCount22 (n * 2) points d S ω) ≤
      (turboOldUniformFactor22 ^ S.card / z ^ d) ^ α *
        ∑ k : Fin 11, ((pairedRegionDP22 B α n S.card) ^ 64) 0 k := by
  let c := turboOldUniformFactor22 ^ S.card / z ^ d
  let P := turboPotentialCountLaw22 (n * 2) S
  let K := fun counts : TurboPotentialCounts22 (n * 2) =>
    ∑ q : State22, turboCountKernel22 counts z 0 q
  let F := fun counts : TurboPotentialCounts22 (n * 2) =>
    ∑ k : Fin 11, (List.ofFn (fun region : Fin 64 =>
      regionalPairedMatrix22 B α n (counts region))).prod 0 k
  have hc : 0 ≤ c := div_nonneg (pow_nonneg turboOldUniformFactor22_pos.le _) (pow_nonneg hz0.le _)
  have hk (counts) : 0 ≤ K counts := by
    change 0 ≤ ∑ q : State22, turboCountKernel22 counts z 0 q
    rw [← turboCountStepArrayMoment22 counts z]
    exact FinPMF.expect_nonneg _ (fun _ => pow_nonneg hz0.le _)
  have hstart : (turboSetup22 (n * 2) points).joint.prob
      (fun ω => 1 ≤ turboActiveBadCount22 (n * 2) points d S ω) ≤
      P.expect (fun counts => (c * K counts) ^ α) := by
    rw [turboSetup22_count_disintegration (n * 2) points S]
    apply bind_prob_count_le_fractional_mean_on_support _ _ _ _ hα0 hα1
    intro counts hp
    have h := turboActiveBadCount22_fiber_le_kernel (n * 2) points hpoints S counts hp d z hz0 hz1
    convert h using 1 <;> dsimp [c, K] <;> ring
  have hstep : P.expect (fun counts => (c * K counts) ^ α) ≤
      P.expect (fun counts => c ^ α * F counts) := by
    apply FinPMF.expect_mono
    intro counts
    rw [Real.mul_rpow hc (hk counts)]
    exact mul_le_mul_of_nonneg_left
      (turboCountKernel22_fractional_le_pairs counts z B hs hB hα0 hα1)
      (Real.rpow_nonneg hc _)
  apply hstart.trans (hstep.trans ?_)
  have he : P.expect F = ∑ k : Fin 11, ((pairedRegionDP22 B α n S.card) ^ 64) 0 k := by
    rw [show P.expect F = ∑ k : Fin 11, P.expect (fun counts =>
      (List.ofFn (fun region : Fin 64 => regionalPairedMatrix22 B α n (counts region))).prod 0 k) by
        dsimp [F]
        rw [FinPMF.expect_sum]]
    apply sum_congr rfl
    intro k _
    exact congrFun (congrFun (turboPairedProfiles22_expect_eq_dp_power B α n S) 0) k
  have hfactor : P.expect (fun counts => c ^ α * F counts) = c ^ α * P.expect F := by
    simp only [FinPMF.expect, Finset.mul_sum]
    apply sum_congr rfl
    intro counts _
    ring
  rw [hfactor, he]

#print axioms turboActiveBadCount22_prob_le_paired_dp
end Spin.Turbo

