import SpinCodes.Structured.ConcreteFixedLargeGapLaw
import SpinCodes.Structured.ConcreteFixedLargeSpacings
import SpinCodes.Structured.ConcretePlacementLimitSiteError
import Mathlib.Data.Nat.Choose.Bounds

/-! Exchangeability passes from the exact finite gaps to the actual ordered-site integral. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Filter MeasureTheory
open scoped Topology

def emptyLiveDuration {R Q : ℕ} (S : Finset (Fin (Q+1))) (B : BlockSubset R Q) : ℝ :=
  (∑ i ∈ S, (emptyGaps B i : ℝ))/(R:ℝ)

theorem liveDuration_empty_distance {R Q : ℕ} (hR : 0 < R)
    (S : Finset (Fin (Q+1))) (B : BlockSubset R Q) :
    |liveDuration S (normalizedSites B) - emptyLiveDuration S B| ≤ (Q:ℝ)/R := by
  unfold liveDuration emptyLiveDuration
  rw [sum_div, ← sum_sub_distrib]
  calc _ ≤ ∑ i ∈ S, |positionGaps (normalizedSites B) i - (emptyGaps B i:ℝ)/R| := abs_sum_le_sum_abs _ _
    _ ≤ ∑ i, |positionGaps (normalizedSites B) i - (emptyGaps B i:ℝ)/R| :=
      sum_le_sum_of_subset_of_nonneg (subset_univ S) (fun _ _ _ => abs_nonneg _)
    _ = (Q:ℝ)/R := by simpa only [abs_sub_comm] using positionGaps_normalizedSites_distance hR B

theorem ordered_block_lattice_mass_le_one {R Q : ℕ} (hR : 0 < R) :
    (1/(R:ℝ)^Q) * (∑ _B : BlockSubset R Q, (1:ℝ)) ≤ 1 := by
  simp only [sum_const, card_univ, card_blockSubset, nsmul_eq_mul, mul_one, one_div, inv_mul_eq_div]
  apply (div_le_one (by positivity : 0 < (R:ℝ)^Q)).mpr
  exact_mod_cast Nat.choose_le_pow R Q

theorem liveDuration_lattice_error {R Q : ℕ} (hR : 0 < R)
    (S : Finset (Fin (Q+1))) {F : ℝ → ℝ} {K : NNReal} (hF : LipschitzWith K F) :
    |(∑ B : BlockSubset R Q, F (liveDuration S (normalizedSites B)))/(R:ℝ)^Q -
      (∑ B : BlockSubset R Q, F (emptyLiveDuration S B))/(R:ℝ)^Q| ≤ (K:ℝ)*((Q:ℝ)/R) := by
  have h := weighted_restricted_sum_error (fun _ : BlockSubset R Q => True)
    (fun B => F (liveDuration S (normalizedSites B))) (fun B => F (emptyLiveDuration S B))
    (w := 1/(R:ℝ)^Q) (ε := (K:ℝ)*((Q:ℝ)/R)) (by positivity) (by positivity)
    (by simpa only [ite_true] using ordered_block_lattice_mass_le_one (Q := Q) hR)
    (fun B _ => ?_)
  · simpa only [ite_true, one_div, inv_mul_eq_div] using h
  · simpa only [Real.dist_eq] using hF.dist_le_mul_of_le (show dist (liveDuration S (normalizedSites B))
        (emptyLiveDuration S B) ≤ (Q:ℝ)/R by
        simpa only [Real.dist_eq] using liveDuration_empty_distance hR S B)

theorem sum_emptyLiveDuration_permutation {R Q : ℕ} (π : Equiv.Perm (Fin (Q+1)))
    (S : Finset (Fin (Q+1))) (F : ℝ → ℝ) :
    (∑ B : BlockSubset R Q, F (emptyLiveDuration (S.map π.toEmbedding) B)) =
      ∑ B : BlockSubset R Q, F (emptyLiveDuration S B) := by
  have h := sum_emptyGaps_permutation (R := R) π (fun g => F ((∑ i ∈ S, (g i:ℝ))/(R:ℝ)))
  simpa only [emptyLiveDuration, sum_map, Equiv.toEmbedding_apply] using h

/-- Actual continuum gap sums are exchangeable, proved from the discrete placement law. -/
theorem liveDuration_integral_permutation {Q : ℕ} (π : Equiv.Perm (Fin (Q+1)))
    (S : Finset (Fin (Q+1))) {F : ℝ → ℝ} {K : NNReal} (hF : LipschitzWith K F) :
    (∫ x in orderedSiteDomain Q, F (liveDuration (S.map π.toEmbedding) x)) =
      ∫ x in orderedSiteDomain Q, F (liveDuration S x) := by
  have ht1 := orderedSite_sum_tendsto (fun x => F (liveDuration (S.map π.toEmbedding) x))
    (hF.continuous.comp (continuous_liveDuration _))
  have ht2 := orderedSite_sum_tendsto (fun x => F (liveDuration S x))
    (hF.continuous.comp (continuous_liveDuration _))
  have hz : Tendsto (fun R : ℕ => 2*(K:ℝ)*((Q:ℝ)/R)) atTop (𝓝 0) := by
    have hh : Tendsto (fun R : ℕ => (Q:ℝ)/(R:ℝ)) atTop (𝓝 0) :=
      tendsto_const_nhds.div_atTop tendsto_natCast_atTop_atTop
    simpa only [mul_zero] using hh.const_mul (2*(K:ℝ))
  have he : ∀ᶠ R : ℕ in atTop,
      |(∑ B : BlockSubset R Q, F (liveDuration (S.map π.toEmbedding) (normalizedSites B)))/(R:ℝ)^Q -
      (∑ B : BlockSubset R Q, F (liveDuration S (normalizedSites B)))/(R:ℝ)^Q| ≤
        2*(K:ℝ)*((Q:ℝ)/R) := by
    filter_upwards [eventually_gt_atTop 0] with R hR
    have h1 := liveDuration_lattice_error hR (S.map π.toEmbedding) hF
    have h2 := liveDuration_lattice_error hR S hF
    rw [sum_emptyLiveDuration_permutation] at h1
    have ht := abs_sub_le
      ((∑ B : BlockSubset R Q, F (liveDuration (S.map π.toEmbedding) (normalizedSites B)))/(R:ℝ)^Q)
      ((∑ B : BlockSubset R Q, F (emptyLiveDuration S B))/(R:ℝ)^Q)
      ((∑ B : BlockSubset R Q, F (liveDuration S (normalizedSites B)))/(R:ℝ)^Q)
    rw [abs_sub_comm ((∑ B : BlockSubset R Q, F (emptyLiveDuration S B))/(R:ℝ)^Q)] at ht
    nlinarith
  have hh := le_of_tendsto_of_tendsto (ht1.sub ht2).abs hz he
  exact sub_eq_zero.mp (abs_eq_zero.mp (le_antisymm hh (abs_nonneg _)))

#print axioms liveDuration_integral_permutation
end Spin.Structured.Placement
