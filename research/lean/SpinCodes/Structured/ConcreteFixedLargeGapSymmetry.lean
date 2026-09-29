import SpinCodes.Structured.ConcreteFixedLargeGapLimit
import SpinCodes.Structured.ConcreteShufflePermutation
import Mathlib.Analysis.Calculus.ContDiff.RCLike

/-! Equal-cardinality live-spacing sets have the same actual continuum law. -/
noncomputable section
namespace Spin.Structured.Placement
open Finset Set MeasureTheory
attribute [local instance] Classical.propDecidable

theorem liveDuration_integral_permutation_on {Q : ℕ} (π : Equiv.Perm (Fin (Q+1)))
    (S : Finset (Fin (Q+1))) {F : ℝ → ℝ} {K : NNReal} (hF : LipschitzOnWith K F (Icc 0 1)) :
    (∫ x in orderedSiteDomain Q, F (liveDuration (S.map π.toEmbedding) x)) =
      ∫ x in orderedSiteDomain Q, F (liveDuration S x) := by
  obtain ⟨G,hG,hFG⟩ := hF.extend_real
  have he (T : Finset (Fin (Q+1))) :
      (∫ x in orderedSiteDomain Q, F (liveDuration T x)) =
        ∫ x in orderedSiteDomain Q, G (liveDuration T x) := by
    apply setIntegral_congr_fun (orderedSiteDomain_measurable Q)
    intro x hx
    exact hFG (liveDuration_mem T hx)
  rw [he,he]
  exact liveDuration_integral_permutation π S hG

theorem liveDuration_integral_card_eq {Q : ℕ} (S T : Finset (Fin (Q+1)))
    (hST : S.card = T.card) {F : ℝ → ℝ} (hF : ContDiffOn ℝ 1 F (Icc 0 1)) :
    (∫ x in orderedSiteDomain Q, F (liveDuration S x)) =
      ∫ x in orderedSiteDomain Q, F (liveDuration T x) := by
  obtain ⟨K,hK⟩ := hF.exists_lipschitzOnWith (by norm_num) (convex_Icc _ _) isCompact_Icc
  obtain ⟨π,hπ⟩ := exists_perm_comp (v := Routing.supportBool S) (v' := Routing.supportBool T)
    (by simpa only [Routing.wtF_supportBool] using hST)
  have hm : S.map π.symm.toEmbedding = T := (Routing.supportBool_shuffle S T π).mp hπ
  have he := liveDuration_integral_permutation_on π.symm S hK
  rw [hm] at he
  exact he.symm

private theorem fin_Iic_step {Q : ℕ} (i : Fin Q) :
    (Finset.Iic i.succ) = insert i.succ (Finset.Iic i.castSucc) := by
  ext j
  simp only [Finset.mem_Iic, Finset.mem_insert]
  constructor
  · intro h; by_cases he : j = i.succ
    · exact Or.inl he
    · right; have := j.isLt; change j.val ≤ i.val; change j.val ≤ i.val+1 at h
      have hn : j.val ≠ i.val+1 := fun hh => he (Fin.ext hh)
      omega
  · rintro (rfl|h)
    · exact le_rfl
    · exact h.trans (Fin.castSucc_lt_succ.le)

theorem liveDuration_prefix {Q : ℕ} (x : Fin Q → ℝ) (i : Fin (Q+1)) :
    liveDuration (Finset.Iic i) x = (Fin.snoc x (1:ℝ) : Fin (Q+1) → ℝ) i := by
  refine Fin.induction ?_ (fun j ih => ?_) i
  · have h : Finset.Iic (0:Fin (Q+1)) = {0} := by
      ext j
      simp only [Finset.mem_Iic, Finset.mem_singleton]
      constructor
      · intro hj; apply Fin.ext; exact Nat.eq_zero_of_le_zero hj
      · rintro rfl; exact le_rfl
    simp only [liveDuration,h,sum_singleton,positionGaps,Fin.cons_zero,sub_zero]
  · unfold liveDuration at *
    rw [fin_Iic_step, sum_insert (by simp), ih]
    simp only [positionGaps,Fin.cons_succ,Fin.snoc_castSucc]
    ring

theorem liveDuration_prefix_site {Q : ℕ} (x : Fin Q → ℝ) (i : Fin Q) :
    liveDuration (Finset.Iic i.castSucc) x = x i := by
  simpa only [Fin.snoc_castSucc] using liveDuration_prefix x i.castSucc

#print axioms liveDuration_integral_card_eq
#print axioms liveDuration_prefix_site
end Spin.Structured.Placement
