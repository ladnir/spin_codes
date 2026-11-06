import SpinCodes.UniformMixer.RegionalEnvelopeComparison
import SpinCodes.Structured.ConcretePlacementLimitSiteError

/-!
# Averaging the regional continuum comparison

The fixed-cardinality shuffle is uniform. We first extend the good-site
pointwise bound to every support by using nonnegativity of the ideal
integrand, and then average over all supports.
-/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured.Placement

attribute [local instance] Classical.propDecidable

theorem continuumSiteProduct_endpoint_nonneg {R a : ℕ} (hR : 0 < R)
    (B : BlockSubset (64 * R) a) {γ : ℝ} (hγ : 0 ≤ γ)
    (q r : State22) :
    0 ≤ (liftM * continuumSiteProduct γ (normalizedSites B) * projectM) q r := by
  have h64R : 0 < 64 * R := by omega
  have hK : Spin.Structured.FiniteKernel.Substochastic
      (continuumSiteProduct γ (normalizedSites B)) := by
    exact timeRegionProduct_substochastic hγ _ _
      (fun i => positionGaps_normalizedSites_nonneg h64R B i.castSucc)
      (positionGaps_normalizedSites_nonneg h64R B (Fin.last a))
  have hJ (b : Bool) : 0 ≤ liftM q b := by
    change 0 ≤ classLift q b
    cases b
    · change 0 ≤ (if q = 0 then (1 : ℝ) else 0)
      split_ifs <;> norm_num
    · change 0 ≤ (if q = 0 then (0 : ℝ) else 1)
      split_ifs <;> norm_num
  have hJK (b : Bool) :
      0 ≤ (liftM * continuumSiteProduct γ (normalizedSites B)) q b := by
    rw [Matrix.mul_apply]
    exact Finset.sum_nonneg (fun c _ => mul_nonneg (hJ c) (hK.nonneg c b))
  rw [Matrix.mul_apply]
  exact Finset.sum_nonneg (fun b _ =>
    mul_nonneg (hJK b) (projectM_unit_interval b r).1)

theorem regionalEnvelope_all_sites_bound {R a : ℕ} (hR : 0 < R)
    (B : BlockSubset (64 * R) a) (θ : ℝ) (hθ : 0 ≤ θ)
    (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    (if ¬ badRegionalSites B.val then
      coarseRegionalEnvelope θ B.val q r else 0) ≤
      (1 + C) ^ (a + 1) *
        ((liftM * continuumSiteProduct (regionGamma θ)
          (normalizedSites B) * projectM) q r +
          regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))) := by
  by_cases hg : ¬ badRegionalSites B.val
  · simpa only [if_pos hg] using
      good_coarseRegionalEnvelope_continuum hR B hg θ hθ q r
  · simp only [if_neg hg]
    have hI := continuumSiteProduct_endpoint_nonneg hR B
      (regionGamma_nonneg hθ) q r
    have hRp : (0 : ℝ) < R := by exact_mod_cast hR
    have hδ : 0 ≤ regionGamma θ *
        (3 * (a + 1 : ℕ) / (R : ℝ)) :=
      mul_nonneg (regionGamma_nonneg hθ) (by positivity)
    have hC : 0 ≤ ((64 * θ + 64 * θ ^ 2) /
        (64 * (R : ℝ))) / Real.exp (-θ) := by
      positivity
    have hpow : 0 ≤ (1 + ((64 * θ + 64 * θ ^ 2) /
        (64 * (R : ℝ))) / Real.exp (-θ)) ^ (a + 1) := by
      positivity
    exact mul_nonneg hpow (add_nonneg hI hδ)

theorem coarseRegionalEnvelope_average_bound {R a : ℕ}
    (hR : 0 < R) (ha : a ≤ 64 * R)
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    (∑ B : BlockSubset (64 * R) a,
      if ¬ badRegionalSites B.val then
        coarseRegionalEnvelope θ B.val q r else 0) /
        (((64 * R).choose a : ℕ) : ℝ) ≤
      (1 + C) ^ (a + 1) *
        ((∑ B : BlockSubset (64 * R) a,
          (liftM * continuumSiteProduct (regionGamma θ)
            (normalizedSites B) * projectM) q r) /
            (((64 * R).choose a : ℕ) : ℝ) +
          regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))) := by
  let D : ℝ := (((64 * R).choose a : ℕ) : ℝ)
  let c : ℝ :=
    (1 + ((64 * θ + 64 * θ ^ 2) /
      (64 * (R : ℝ))) / Real.exp (-θ)) ^ (a + 1)
  let ε : ℝ := regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))
  let F : BlockSubset (64 * R) a → ℝ :=
    fun B => (liftM * continuumSiteProduct (regionGamma θ)
      (normalizedSites B) * projectM) q r
  have hD : 0 < D := by
    dsimp [D]
    exact_mod_cast Nat.choose_pos ha
  have hsum :
      (∑ B : BlockSubset (64 * R) a,
        if ¬ badRegionalSites B.val then
          coarseRegionalEnvelope θ B.val q r else 0) ≤
      ∑ B : BlockSubset (64 * R) a, c * (F B + ε) := by
    apply Finset.sum_le_sum
    intro B _
    exact regionalEnvelope_all_sites_bound hR B θ hθ q r
  have hsumdiv := div_le_div_of_nonneg_right hsum hD.le
  have hnorm :
      (∑ B : BlockSubset (64 * R) a, c * (F B + ε)) / D =
        c * ((∑ B : BlockSubset (64 * R) a, F B) / D + ε) := by
    have hcard : Fintype.card (BlockSubset (64 * R) a) =
        (64 * R).choose a := card_blockSubset _ _
    simp only [mul_add, Finset.sum_add_distrib, ← Finset.mul_sum,
      Finset.sum_const, Finset.card_univ, hcard, nsmul_eq_mul]
    have hDne : (((64 * R).choose a : ℕ) : ℝ) ≠ 0 := by
      exact_mod_cast Nat.choose_ne_zero ha
    have hDne' : (((R * 64).choose a : ℕ) : ℝ) ≠ 0 := by
      simpa only [Nat.mul_comm] using hDne
    field_simp [hDne, hDne']
    ring
  dsimp [D, c, ε, F] at hsumdiv hnorm
  simpa only [hnorm] using hsumdiv

theorem averagedRegionKernel_continuum_finite {R a : ℕ}
    (hR : 0 < R) (ha : a ≤ 64 * R)
    (S : Finset (Fin (64 * R))) (hS : S.card = a)
    (θ : ℝ) (hθ : 0 ≤ θ) (q r : State22) :
    let L : ℝ := 64 * R
    let C := ((64 * θ + 64 * θ ^ 2) / L) / Real.exp (-θ)
    averagedRegionKernel S (Real.exp (-(θ / L))) q r ≤
      (1 + C) ^ (a + 1) *
        ((∑ B : BlockSubset (64 * R) a,
          (liftM * continuumSiteProduct (regionGamma θ)
            (normalizedSites B) * projectM) q r) /
            (((64 * R).choose a : ℕ) : ℝ) +
          regionGamma θ * (3 * (a + 1 : ℕ) / (R : ℝ))) +
      128 * (a : ℝ) / L +
      257 * (a : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1) := by
  have hbase := averagedRegionKernel_coarse_layer hR S hS θ hθ q r
  have hcompare := coarseRegionalEnvelope_average_bound hR ha θ hθ q r
  have hsum := add_le_add_left hcompare
    (128 * (a : ℝ) / (64 * (R : ℝ)) +
      257 * (a : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1))
  have hbase' :
      averagedRegionKernel S (Real.exp (-(θ / (64 * (R : ℝ))))) q r ≤
      (∑ T : BlockSubset (64 * R) a,
        if ¬ badRegionalSites T.val then
          coarseRegionalEnvelope θ T.val q r else 0) /
          (((64 * R).choose a : ℕ) : ℝ) +
        (128 * (a : ℝ) / (64 * (R : ℝ)) +
          257 * (a : ℝ) ^ 2 / (((64 * R : ℕ) : ℝ) - 1)) := by
    simpa only [add_assoc] using hbase
  have hfinal := hbase'.trans hsum
  simpa only [add_assoc] using hfinal

end Spin.UniformMixer
