import SpinCodes.UniformMixer.ScalableFixedLargeRoute
import SpinCodes.UniformMixer.ScalableOuterCounting
import SpinCodes.Structured.DenseOccupationProfilePartition

/-! Exact selected-outer profile count for a Scalable routed functional. -/

noncomputable section

namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteRoute
  Spin.Structured.DenseOccupationFixed Finset

theorem scalable_route_expect_eq_of_card_eq {L b : ℕ}
    (rows rows' : Fin L → Finset (Fin b))
    (h : ∀ i, (rows i).card = (rows' i).card)
    (F : (Fin b → Finset (Fin L)) → ℝ) :
    (law rows).expect F = (law rows').expect F := by
  have hr : rowLaw rows = rowLaw rows' := by
    apply FinPMF.ext
    funext out
    simp only [rowLaw, piPMF_apply]
    apply Finset.prod_congr rfl
    intro i _
    simp only [Routing.shuffleLaw_apply, h i]
  unfold law
  rw [hr]

theorem scalable_profile_failure_eq (m : ℕ)
    (seed : NativeSeed m)
    (S : Finset (Fin (Lsched m)))
    (rows : Fin (Lsched m) →
      Finset (Fin (nativeBlocks m * 24)))
    (hz : ∀ i, i ∉ S → rows i = ∅) (d : ℕ) :
    (∑ x ∈ profileMessages seed S
        (fun i => (rows i).card),
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d)) =
      (∏ i ∈ S,
        (spectrum seed (rows i).card : ℝ)) *
        (law (fun i => (rows i).map
          (finCongr (native_width m)).toEmbedding)).expect
            (scalableRegionFailure m d) := by
  have he : ∀ x ∈ profileMessages seed S
      (fun i => (rows i).card),
      (law (fun i => (rowSupports seed x i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d) =
      (law (fun i => (rows i).map
        (finCongr (native_width m)).toEmbedding)).expect
          (scalableRegionFailure m d) := by
    intro x hx
    apply scalable_route_expect_eq_of_card_eq
    intro i
    simp only [Finset.card_map]
    have hi := (Finset.mem_filter.mp hx).2 i
    by_cases hs : i ∈ S
    · simp only [hs, ite_true] at hi
      simpa only [rowSupports, support_card] using hi.2
    · simp only [hs, ite_false] at hi
      simp only [rowSupports, hi, encode_zero,
        hz i hs, support_card, Finset.card_empty]
      exact (wtF_eq_zero _).mpr rfl
  rw [Finset.sum_congr rfl he,
    Finset.sum_const, nsmul_eq_mul,
    profileMessages_card]

theorem scalable_route_probability_eq (m d : ℕ)
    (rows : Fin (Lsched m) → Finset (Fin (bsched m))) :
    (law rows).expect (scalableRegionFailure m d) =
      ((law rows).prod (scalableStateMapLaw m)).prob
        (fun ω => sampledOutputWeight
          (scalableFlatInput m ω.1)
          (scalableFlatMaps m ω.2) 0 ≤ d) := by
  exact (FinPMF.prob_prod_eq_expect _ _ _).symm

end Spin.UniformMixer
