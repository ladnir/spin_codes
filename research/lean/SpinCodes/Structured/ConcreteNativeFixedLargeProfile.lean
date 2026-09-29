import SpinCodes.Structured.ConcreteNativeFixedLargeCount
import SpinCodes.Structured.ConcreteNativeFixedLargeNativeRoute
import SpinCodes.Structured.DenseOccupationFamilyRateProfiles

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open Finset Filter ConcreteOuter ConcreteRoute ConcreteFixedNumeric Placement DenseOccupationFixed

theorem native_fixed_large_profile_bound {Q : ℕ} (hQ : 0 < Q) (hK : FixedLargeNorm Q) :
    ∀ᶠ m in atTop, ∀ seed, nativeGood m seed → ∀ S : Finset (Fin (Lsched m)), S.card=Q →
      ∀ w : WeightProfile S (nativeBlocks m*24),
      (∑ x∈profileMessages seed S (profileWeights S w),
        failureProbability (tupleWiring m) (rowSupports seed x) (threshold m)) ≤
          fixedProfileCountBound Q (nativeBlocks m*24) (Lsched m) (threshold m) := by
  classical
  filter_upwards [native_fixed_profile_routed hQ hK] with m hroute
  intro seed hg S hS w
  have hk : 0 < nativeBlocks m := by
    have h := bsched_pos m
    rw [←native_width m] at h
    omega
  by_cases hempty : (profileMessages seed S (profileWeights S w)).Nonempty
  · obtain ⟨message,hm⟩ := hempty
    have hact : message∈activeMessages S := by
      have hf : message∈(activeMessages S).filter (fun x => messageProfile seed S x=w) := by
        rw [profile_fiber_eq]
        exact hm
      exact (mem_filter.mp hf).1
    have hs := selected_active_rows hk seed hg S hact
    let e : Fin Q ↪ Fin (Lsched m) := (S.orderEmbOfFin hS).toEmbedding
    have he : univ.map e=S := by
      rw [Finset.map_eq_image]
      exact Finset.image_orderEmbOfFin_univ S hS
    have hei (i : Fin Q) : e i∈S := by
      rw [←he]
      exact mem_map.mpr ⟨i,mem_univ _,rfl⟩
    have hoff : ∀ j, j∉Set.range e → rowSupports seed message j=∅ := by
      intro j hj
      apply hs.1 j
      intro hjS
      rw [←he] at hjS
      obtain ⟨i,_,hi⟩ := mem_map.mp hjS
      exact hj ⟨i,hi⟩
    let rows := fun i => rowSupports seed message (e i)
    have hemb : embedRows e rows=rowSupports seed message := embedRows_extract e _ hoff
    obtain ⟨u,hu⟩ := choose_fixed_fugacities (fun i => ((rows i).card:ℝ)/(nativeBlocks m*24:ℕ))
      (fun i => hs.2.2 (e i) (hei i))
    have hp := hroute u e rows
    rw [hemb] at hp
    have hc := fixed_counted_profile_bound hk seed (weightWindow (nativeBlocks m*24)) hg rows u
      (fun i => hs.2.1 (e i) (hei i)) (fun i => hs.2.2 (e i) (hei i)) hu hp
    have hw : ∀ i∈S, (rowSupports seed message i).card=profileWeights S w i := by
      intro i hi
      have h := (mem_filter.mp hm).2 i
      simp only [hi,ite_true] at h
      simpa only [rowSupports,support_card] using h.2
    rw [←profileMessages_congr seed S hw,profile_failure_eq seed S _ hs.1 (tupleWiring m) (threshold m)]
    have hprod : (∏ i∈S, (spectrum seed (rowSupports seed message i).card:ℝ))=
        ∏ i, (spectrum seed (rows i).card:ℝ) := by
      rw [←he,Finset.prod_map]
    rw [hprod]
    exact hc
  · rw [Finset.not_nonempty_iff_eq_empty.mp hempty,Finset.sum_empty]
    unfold fixedProfileCountBound fixedShellCost
    positivity

#print axioms native_fixed_large_profile_bound
end Spin.Structured.ConcreteNativeFamily
