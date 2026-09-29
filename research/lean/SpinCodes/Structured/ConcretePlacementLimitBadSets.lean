import SpinCodes.Structured.ConcretePlacementLimitSiteError

/-! Removing bad block sets from the coarse Riemann sum costs their actual shuffle probability. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Routing ConcreteEncoder FiniteKernel
attribute [local instance] Classical.propDecidable

theorem weighted_bad_placement_mass {R a : Nat} (S : Finset (Fin (128 * R))) (hS : S.card = a) (H : Nat) :
    ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
      (∑ B : BlockSubset R a, if badBlockIndices (orderedBlocks B) H then (1 : ℝ) else 0) ≤
        (shuffleLaw S).prob (badBlockPlacement H) := by
  have he := shuffle_distinct_blocks_sum S hS (fun T => if badBlockPlacement H T then (1 : ℝ) else 0)
  simp only [badBlockPlacement_singletonSupport, Spin.FinPMF.expect_const] at he
  rw [← he, Spin.FinPMF.prob_eq_expect_indicator]
  apply Spin.FinPMF.expect_mono
  intro T
  split_ifs <;> norm_num

theorem weighted_omit_sum_le {Ω : Type*} [Fintype Ω] (E : Ω → Prop) [DecidablePred E]
    (f : Ω → ℝ) (hf : ∀ x, 0 ≤ f x ∧ f x ≤ 1) {w : ℝ} (hw : 0 ≤ w) :
    |w * (∑ x, if ¬ E x then f x else 0) - w * ∑ x, f x| ≤
      w * (∑ x, if E x then (1 : ℝ) else 0) := by
  rw [abs_sub_comm, ← mul_sub, abs_mul, abs_of_nonneg hw, ← sum_sub_distrib]
  apply mul_le_mul_of_nonneg_left _ hw
  calc _ ≤ ∑ x, |f x - if ¬ E x then f x else 0| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ _ := by
      apply sum_le_sum
      intro x _
      by_cases hx : E x
      · simpa only [hx, not_true_eq_false, ite_false, ite_true, sub_zero, abs_of_nonneg (hf x).1] using (hf x).2
      · simp [hx]

theorem substochastic_entry_unit {α : Type*} [Fintype α] [DecidableEq α]
    {K : Matrix α α ℝ} (hK : Substochastic K) (q r : α) : 0 ≤ K q r ∧ K q r ≤ 1 :=
  ⟨hK.nonneg q r, (Finset.single_le_sum (fun s _ => hK.nonneg q s) (mem_univ r)).trans (hK.row_le q)⟩

theorem lifted_entry_unit {K : Matrix (Fin 2) (Fin 2) ℝ} (hK : Substochastic K) (q r : State) :
    0 ≤ (liveLift * K * liveProjection) q r ∧ (liveLift * K * liveProjection) q r ≤ 1 := by
  rw [lifted_entry]
  split_ifs
  · exact substochastic_entry_unit hK _ _
  · have h := substochastic_entry_unit hK (stateClass q) 1
    exact ⟨div_nonneg h.1 (by norm_num), (div_le_self h.1 (by norm_num)).trans h.2⟩

theorem siteProduct_normalized_substochastic {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ) (B : BlockSubset R a) :
    Substochastic (siteProduct (θ * epochMean / 128) (normalizedSites B)) := by
  apply timeProduct_substochastic (by unfold epochMean; positivity)
  · exact fun i => positionGaps_normalizedSites_nonneg hR B i.castSucc
  · exact positionGaps_normalizedSites_nonneg hR B (Fin.last a)

def finiteGoodKernel (θ : ℝ) (R a H : Nat) (q r : State) : ℝ :=
  ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
    ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then
      (liveLift * placementProduct θ B * liveProjection) q r else 0

def finiteSiteKernel (θ : ℝ) (R a : Nat) (q r : State) : ℝ :=
  ((128 : ℝ)^a / ((128 * R).choose a : ℝ)) *
    ∑ B : BlockSubset R a,
      (liveLift * siteProduct (θ * epochMean / 128) (normalizedSites B) * liveProjection) q r

theorem finiteGoodKernel_site_error {R a : Nat} (hR : 0 < R) {θ : ℝ} (hθ : 0 ≤ θ)
    (S : Finset (Fin (128 * R))) (hS : S.card = a) (H : Nat) (q r : State) :
    |finiteGoodKernel θ R a H q r - finiteSiteKernel θ R a q r| ≤
      (θ * epochMean / 128) * ((a : ℝ) / R) + placementOmissionError a R H := by
  let w : ℝ := (128 : ℝ)^a / ((128 * R).choose a : ℝ)
  let G : BlockSubset R a → ℝ := fun B =>
    (liveLift * siteProduct (θ * epochMean / 128) (normalizedSites B) * liveProjection) q r
  let M : ℝ := w * ∑ B : BlockSubset R a, if ¬ badBlockIndices (orderedBlocks B) H then G B else 0
  have heps : 0 ≤ (θ * epochMean / 128) * ((a : ℝ) / R) := by unfold epochMean; positivity
  have hfirst : |finiteGoodKernel θ R a H q r - M| ≤ (θ * epochMean / 128) * ((a : ℝ) / R) := by
    exact weighted_restricted_sum_error (fun B : BlockSubset R a => ¬ badBlockIndices (orderedBlocks B) H)
      (fun B => (liveLift * placementProduct θ B * liveProjection) q r) G (w := w) (by positivity) heps
      (good_placement_mass_le_one S hS H)
      (fun B _ => lifted_rowError_entry heps (placementProduct_site_error hR hθ B) q r)
  have hsecond : |M - finiteSiteKernel θ R a q r| ≤ placementOmissionError a R H := by
    have h := weighted_omit_sum_le (fun B : BlockSubset R a => badBlockIndices (orderedBlocks B) H)
      G (fun B => lifted_entry_unit (siteProduct_normalized_substochastic hR hθ B) q r) (w := w) (by positivity)
    have hb := (weighted_bad_placement_mass S hS H).trans (shuffle_prob_badBlockPlacement hR S H)
    rw [hS] at hb
    exact h.trans hb
  exact (abs_sub_le _ M _).trans (add_le_add hfirst hsecond)

end Spin.Structured.Placement
