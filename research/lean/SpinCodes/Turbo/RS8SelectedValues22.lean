import SpinCodes.Turbo.RS8ReferencePlacement22

/-!
# Uniform four-bit values on a selected regional set

This is the variable-length version of `PhysicalThinning22.selectedGiven22`.
The selected set is kept as latent data: a selected coordinate may carry
zero, so it cannot be recovered from the output's nonzero support.
-/

noncomputable section

namespace Spin.Turbo

open Finset

def regionalSelectedLabelWord22 {n : ℕ}
    (T : Finset (Fin n)) (a : T → Packet4) : Fin n → Packet4 :=
  fun i => if h : i ∈ T then a ⟨i, h⟩ else 0

def regionalSelectedGiven22 {n : ℕ}
    (T : Finset (Fin n)) : FinPMF (Fin n → Packet4) :=
  (FinPMF.uniform (T → Packet4)).map
    (regionalSelectedLabelWord22 T)

theorem regionalSelectedLabelWord22_support_subset {n : ℕ}
    (T : Finset (Fin n)) (a : T → Packet4) :
    regionalPacketSupport22 (regionalSelectedLabelWord22 T a) ⊆ T := by
  intro i hi
  by_contra hnot
  have hnonzero := (Finset.mem_filter.mp hi).2
  apply hnonzero
  simp [regionalSelectedLabelWord22, hnot]

theorem regionalSelectedLabelWord22_injective {n : ℕ}
    (T : Finset (Fin n)) :
    Function.Injective (regionalSelectedLabelWord22 T) := by
  intro a b hab
  funext i
  have h := congrFun hab i.val
  simpa [regionalSelectedLabelWord22, i.property] using h

theorem regionalSelectedLabelWord22_exists_iff {n : ℕ}
    (T : Finset (Fin n)) (x : Fin n → Packet4) :
    (∃ a : T → Packet4, regionalSelectedLabelWord22 T a = x) ↔
      regionalPacketSupport22 x ⊆ T := by
  constructor
  · rintro ⟨a, rfl⟩
    exact regionalSelectedLabelWord22_support_subset T a
  · intro hx
    refine ⟨fun i => x i.val, ?_⟩
    funext i
    by_cases hi : i ∈ T
    · simp [regionalSelectedLabelWord22, hi]
    · have hz : x i = 0 := by
        by_contra hn
        exact hi (hx (Finset.mem_filter.mpr
          ⟨Finset.mem_univ _, hn⟩))
      simp [regionalSelectedLabelWord22, hi, hz]

theorem regionalSelectedGiven22_apply {n : ℕ}
    (T : Finset (Fin n)) (x : Fin n → Packet4) :
    (regionalSelectedGiven22 T).p x =
      if regionalPacketSupport22 x ⊆ T then
        1 / (16 : ℝ) ^ T.card else 0 := by
  classical
  rw [regionalSelectedGiven22, FinPMF.map_p]
  by_cases hx : regionalPacketSupport22 x ⊆ T
  · rw [if_pos hx]
    obtain ⟨a, ha⟩ :=
      (regionalSelectedLabelWord22_exists_iff T x).mpr hx
    have he (b : T → Packet4) :
        regionalSelectedLabelWord22 T b = x ↔ b = a := by
      constructor
      · intro hb
        exact regionalSelectedLabelWord22_injective T
          (hb.trans ha.symm)
      · rintro rfl
        exact ha
    simp only [he, Finset.sum_ite_eq', Finset.mem_univ, if_pos]
    change 1 / (Fintype.card (T → Packet4) : ℝ) = _
    simp [Packet4]
  · rw [if_neg hx]
    apply Finset.sum_eq_zero
    intro a _
    have ha : regionalSelectedLabelWord22 T a ≠ x := by
      intro he
      exact hx ((regionalSelectedLabelWord22_exists_iff T x).mp
        ⟨a, he⟩)
    simp [ha]

theorem turboReferenceRegionPre_eq_selected_given
    {R : ℕ} (S : Finset (Fin (R * 16))) :
    Spin.piPMF (fun i : Fin (R * 16) =>
      turboReferenceCellLaw22 S i) =
        regionalSelectedGiven22 S := by
  apply FinPMF.ext
  funext x
  rw [turboReferenceRegionPre_apply, regionalSelectedGiven22_apply]

theorem turboReferenceGiven_shuffle
    {R : ℕ} (S : Finset (Fin (R * 16)))
    (σ : Equiv.Perm (Fin (R * 16))) :
    (Spin.piPMF (fun i : Fin (R * 16) =>
      turboReferenceCellLaw22 S i)).map
        (regionalValueShuffleEquiv22 σ) =
      regionalSelectedGiven22
        (Spin.Structured.Routing.shuffleSupport σ S) := by
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply, turboReferenceRegionPre_apply,
    regionalSelectedGiven22_apply]
  simp only [regionalValueShuffle_subset_iff,
    Spin.Structured.Routing.shuffleSupport_card]

theorem turboReferenceOneRegionLaw22_eq_selected_mixture
    (R : ℕ) (S : Finset (Fin (R * 16))) :
    turboReferenceOneRegionLaw22 R S =
      (Spin.Structured.Routing.shuffleLaw S).bind
        (fun T => regionalSelectedGiven22 T) := by
  rw [turboReferenceOneRegionLaw22_eq_perm_mixture]
  simp_rw [turboReferenceGiven_shuffle]
  rw [← FinPMF.bind_map
    (FinPMF.uniform (Equiv.Perm (Fin (R * 16))))
    (fun σ => Spin.Structured.Routing.shuffleSupport σ S)
    regionalSelectedGiven22]
  rw [Spin.Structured.Routing.uniform_permutation_shuffleLaw]

#print axioms regionalSelectedGiven22_apply
#print axioms turboReferenceOneRegionLaw22_eq_selected_mixture

end Spin.Turbo
