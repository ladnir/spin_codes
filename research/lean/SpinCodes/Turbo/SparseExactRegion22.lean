import SpinCodes.Turbo.RS8SelectedValues22
import SpinCodes.Turbo.PhysicalPacketMoment22

/-!
# Exact nonzero values on the actual sparse regional support

Audit: `PacketShuffleMixture.uniformLabeledLaw_apply` keeps complete values
on an exact support; `RS8SelectedValues22` supplies its variable-length
zero-inclusive counterpart. We use the exact support here, with fifteen
labels at each active slot. The regional permutation acts on the whole
selected set, including both distinct slots in the two-group case.
-/
noncomputable section
namespace Spin.Turbo
open Finset

/-- Nonzero labels fill exactly the selected positions. -/
def regionalNonzeroLabelWord22 {n : ℕ} (T : Finset (Fin n))
    (a : T → {p : Packet4 // p ≠ 0}) : Fin n → Packet4 :=
  fun i => if h : i ∈ T then (a ⟨i, h⟩).val else 0

def regionalNonzeroGiven22 {n : ℕ} (T : Finset (Fin n)) :
    FinPMF (Fin n → Packet4) :=
  (FinPMF.uniform (T → {p : Packet4 // p ≠ 0})).map
    (regionalNonzeroLabelWord22 T)

theorem regionalNonzeroLabelWord22_support {n : ℕ}
    (T : Finset (Fin n)) (a : T → {p : Packet4 // p ≠ 0}) :
    regionalPacketSupport22 (regionalNonzeroLabelWord22 T a) = T := by
  ext i
  simp only [regionalPacketSupport22, Finset.mem_filter, Finset.mem_univ, true_and]
  by_cases hi : i ∈ T
  · simp [regionalNonzeroLabelWord22, hi, (a ⟨i, hi⟩).property]
  · simp [regionalNonzeroLabelWord22, hi]

theorem regionalNonzeroLabelWord22_injective {n : ℕ}
    (T : Finset (Fin n)) : Function.Injective (regionalNonzeroLabelWord22 T) := by
  intro a b hab
  funext i
  apply Subtype.ext
  have h := congrFun hab i.val
  simpa [regionalNonzeroLabelWord22, i.property] using h

theorem regionalNonzeroLabelWord22_exists_iff {n : ℕ}
    (T : Finset (Fin n)) (x : Fin n → Packet4) :
    (∃ a : T → {p : Packet4 // p ≠ 0}, regionalNonzeroLabelWord22 T a = x) ↔
      regionalPacketSupport22 x = T := by
  constructor
  · rintro ⟨a, rfl⟩
    exact regionalNonzeroLabelWord22_support T a
  · intro hx
    have hnonzero (i : T) : x i.val ≠ 0 := by
      have hi : i.val ∈ regionalPacketSupport22 x := by rw [hx]; exact i.property
      exact (Finset.mem_filter.mp hi).2
    refine ⟨fun i => ⟨x i.val, hnonzero i⟩, ?_⟩
    funext i
    by_cases hi : i ∈ T
    · simp [regionalNonzeroLabelWord22, hi]
    · have hz : x i = 0 := by
        by_contra hn
        exact hi (hx ▸ Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩)
      simp [regionalNonzeroLabelWord22, hi, hz]

theorem regionalNonzeroGiven22_apply {n : ℕ}
    (T : Finset (Fin n)) (x : Fin n → Packet4) :
    (regionalNonzeroGiven22 T).p x =
      if regionalPacketSupport22 x = T then 1 / (15 : ℝ) ^ T.card else 0 := by
  rw [regionalNonzeroGiven22, FinPMF.map_p]
  by_cases hx : regionalPacketSupport22 x = T
  · rw [if_pos hx]
    obtain ⟨a, ha⟩ := (regionalNonzeroLabelWord22_exists_iff T x).mpr hx
    have he (b : T → {p : Packet4 // p ≠ 0}) :
        regionalNonzeroLabelWord22 T b = x ↔ b = a := by
      constructor
      · intro hb
        exact regionalNonzeroLabelWord22_injective T (hb.trans ha.symm)
      · rintro rfl
        exact ha
    simp only [he, Finset.sum_ite_eq', Finset.mem_univ, ite_true]
    change 1 / (Fintype.card (T → {p : Packet4 // p ≠ 0}) : ℝ) = _
    rw [Fintype.card_fun, packet4_nonzero_card, Fintype.card_coe]
    norm_cast
  · rw [if_neg hx]
    apply Finset.sum_eq_zero
    intro a _
    have ha : regionalNonzeroLabelWord22 T a ≠ x := by
      intro he
      exact hx ((regionalNonzeroLabelWord22_exists_iff T x).mp ⟨a, he⟩)
    simp [ha]

/-- The same regional permutation transports the entire exact support. -/
theorem regionalNonzeroSupport22_shuffle {R : ℕ}
    (σ : Equiv.Perm (Fin (R * 16))) (x : Fin (R * 16) → Packet4) :
    regionalPacketSupport22 (regionalValueShuffleEquiv22 σ x) =
      Spin.Structured.Routing.shuffleSupport σ (regionalPacketSupport22 x) := by
  ext i
  simp only [regionalPacketSupport22, Finset.mem_filter, Finset.mem_univ, true_and,
    Spin.Structured.Routing.mem_shuffleSupport]
  rfl

/-- Exact nonzero values stay uniform under fixed wiring. Both selected
slots, when present, are transported by this single permutation. -/
theorem regionalNonzeroGiven22_shuffle {R : ℕ}
    (T : Finset (Fin (R * 16))) (σ : Equiv.Perm (Fin (R * 16))) :
    (regionalNonzeroGiven22 T).map (regionalValueShuffleEquiv22 σ) =
      regionalNonzeroGiven22 (Spin.Structured.Routing.shuffleSupport σ T) := by
  apply FinPMF.ext
  funext y
  rw [FinPMF.map_equiv_apply, regionalNonzeroGiven22_apply, regionalNonzeroGiven22_apply]
  have hs : regionalPacketSupport22 ((regionalValueShuffleEquiv22 σ).symm y) = T ↔
      regionalPacketSupport22 y = Spin.Structured.Routing.shuffleSupport σ T := by
    have h := regionalNonzeroSupport22_shuffle σ ((regionalValueShuffleEquiv22 σ).symm y)
    rw [Equiv.apply_symm_apply] at h
    constructor
    · intro he
      exact h.trans (congrArg (Spin.Structured.Routing.shuffleSupport σ) he)
    · intro he
      apply Finset.map_injective σ.symm.toEmbedding
      exact h.symm.trans he
  simp only [hs, Spin.Structured.Routing.shuffleSupport_card]
/-- Average one common regional permutation of an exact labelled support. -/
def sparseExactRegionLaw22 (R : ℕ) (T : Finset (Fin (R * 16))) :
    FinPMF (Fin (R * 16) → Packet4) :=
  (FinPMF.uniform (Equiv.Perm (Fin (R * 16)))).bind
    (fun σ => (regionalNonzeroGiven22 T).map (regionalValueShuffleEquiv22 σ))

theorem sparseExactRegionLaw22_eq_support_mixture
    (R : ℕ) (T : Finset (Fin (R * 16))) :
    sparseExactRegionLaw22 R T =
      (Spin.Structured.Routing.shuffleLaw T).bind regionalNonzeroGiven22 := by
  rw [sparseExactRegionLaw22]
  simp_rw [regionalNonzeroGiven22_shuffle]
  rw [← FinPMF.bind_map
    (FinPMF.uniform (Equiv.Perm (Fin (R * 16))))
    (fun σ => Spin.Structured.Routing.shuffleSupport σ T)
    regionalNonzeroGiven22]
  rw [Spin.Structured.Routing.uniform_permutation_shuffleLaw]

/-- Only occupancy matters after the shared permutation has been averaged. -/
theorem sparseExactRegionLaw22_eq_of_card
    (R : ℕ) (T U : Finset (Fin (R * 16))) (h : T.card = U.card) :
    sparseExactRegionLaw22 R T = sparseExactRegionLaw22 R U := by
  rw [sparseExactRegionLaw22_eq_support_mixture,
    sparseExactRegionLaw22_eq_support_mixture]
  have he : Spin.Structured.Routing.shuffleLaw T =
      Spin.Structured.Routing.shuffleLaw U := by
    apply FinPMF.ext
    funext S
    simp only [Spin.Structured.Routing.shuffleLaw_apply, h]
  rw [he]
#print axioms regionalNonzeroGiven22_apply
#print axioms sparseExactRegionLaw22_eq_of_card
#print axioms regionalNonzeroGiven22_shuffle
end Spin.Turbo


