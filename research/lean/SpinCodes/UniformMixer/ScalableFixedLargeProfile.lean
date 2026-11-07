import SpinCodes.UniformMixer.ScalableFixedLargeCount
import SpinCodes.Structured.DenseOccupationFamilyRateProfiles

/-! The fixed-row profile bound for the actual Scalable route and sampled maps. -/

noncomputable section
namespace Spin.UniformMixer

open Spin.Structured Spin.Structured.ConcreteOuter
  Spin.Structured.ConcreteNativeFamily Spin.Structured.ConcreteRoute
  Spin.Structured.DenseOccupationFixed Spin.Structured.Placement
  Finset Filter

theorem scalable_fixed_large_profile_bound {Q : ℕ} (hQ : 3 ≤ Q) :
    ∀ᶠ m in atTop, ∀ seed, nativeGood m seed →
      ∀ S : Finset (Fin (Lsched m)), S.card = Q →
        ∀ w : WeightProfile S (nativeBlocks m * 24),
          (∑ x ∈ profileMessages seed S (profileWeights S w),
            (law (fun i => (rowSupports seed x i).map
              (finCongr (native_width m)).toEmbedding)).expect
                (scalableRegionFailure m (threshold m))) ≤
            fixedProfileCountBound Q (nativeBlocks m * 24)
              (Lsched m) (threshold m) * scalableTiltCorrection Q m := by
  classical
  filter_upwards [scalable_fixed_large_routed hQ] with m hroute
  intro seed hg S hS w
  have hk : 0 < nativeBlocks m := by
    have h := bsched_pos m
    rw [← native_width m] at h
    omega
  by_cases hempty : (profileMessages seed S (profileWeights S w)).Nonempty
  · obtain ⟨message, hm⟩ := hempty
    have hact : message ∈ activeMessages S := by
      have hf : message ∈ (activeMessages S).filter
          (fun x => messageProfile seed S x = w) := by
        rw [profile_fiber_eq]
        exact hm
      exact (mem_filter.mp hf).1
    have hs := selected_active_rows hk seed hg S hact
    let e : Fin Q ↪ Fin (Lsched m) := (S.orderEmbOfFin hS).toEmbedding
    have he : univ.map e = S := by
      rw [Finset.map_eq_image]
      exact Finset.image_orderEmbOfFin_univ S hS
    have hei (i : Fin Q) : e i ∈ S := by
      rw [← he]
      exact mem_map.mpr ⟨i, mem_univ _, rfl⟩
    have hoff : ∀ j, j ∉ Set.range e → rowSupports seed message j = ∅ := by
      intro j hj
      apply hs.1 j
      intro hjS
      rw [← he] at hjS
      obtain ⟨i, _, hi⟩ := mem_map.mp hjS
      exact hj ⟨i, hi⟩
    let rows := fun i => rowSupports seed message (e i)
    let rowsB := fun i => (rows i).map (finCongr (native_width m)).toEmbedding
    have hembB : embedRows e rowsB =
        fun j => (rowSupports seed message j).map
          (finCongr (native_width m)).toEmbedding := by
      funext j
      by_cases hj : j ∈ Set.range e
      · obtain ⟨i, rfl⟩ := hj
        simp only [embedRows_apply, rowsB, rows]
      · rw [embedRows_off e rowsB j hj, hoff j hj]
        simp
    obtain ⟨u, hu⟩ := choose_fixed_fugacities
      (fun i => ((rows i).card : ℝ) / (nativeBlocks m * 24 : ℕ))
      (fun i => hs.2.2 (e i) (hei i))
    have hp := hroute u e rowsB (threshold m)
    rw [hembB, ← scalable_route_probability_eq] at hp
    have hp' :
        (law (fun j => (rowSupports seed message j).map
          (finCongr (native_width m)).toEmbedding)).expect
            (scalableRegionFailure m (threshold m)) ≤
        (((fixedNormCost u + fixedNormSlack u) ^
              (nativeBlocks m * 24) / min 1 (3 / 1600 : ℝ)) /
            (profileChoices rows *
              profileCoefficient (fixedFugacity u)
                (fun i => (rows i).card))) /
          Real.exp (-(scalableLargeTheta Q / (Lsched m : ℝ))) ^
            (threshold m) := by
      simpa only [rowsB, rows, profileChoices, profileCoefficient,
        Finset.card_map, native_width] using hp
    have hc := scalable_tilted_counted_profile_bound hk seed
      (weightWindow (nativeBlocks m * 24)) hg rows u
      (fun i => hs.2.1 (e i) (hei i))
      (fun i => hs.2.2 (e i) (hei i)) hu
      (scalableLargeTheta Q) hp'
    have hw : ∀ i ∈ S,
        (rowSupports seed message i).card = profileWeights S w i := by
      intro i hi
      have h := (mem_filter.mp hm).2 i
      simp only [hi, ite_true] at h
      simpa only [rowSupports, support_card] using h.2
    rw [← profileMessages_congr seed S hw,
      scalable_profile_failure_eq m seed S _ hs.1 (threshold m)]
    have hprod :
        (∏ i ∈ S,
          (spectrum seed (rowSupports seed message i).card : ℝ)) =
        ∏ i, (spectrum seed (rows i).card : ℝ) := by
      rw [← he, Finset.prod_map]
    rw [hprod]
    simpa only [scalableTiltCorrection] using hc
  · rw [Finset.not_nonempty_iff_eq_empty.mp hempty, Finset.sum_empty]
    unfold fixedProfileCountBound fixedShellCost scalableTiltCorrection
    positivity

#print axioms scalable_fixed_large_profile_bound
end Spin.UniformMixer
