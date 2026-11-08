import SpinCodes.Turbo.PhysicalStep22
import SpinCodes.Structured.ConcretePlacementLaw

/-!
# Zero-inclusive selected slots at a Turbo physical step

The numerical producer first selects `j` packet slots, then gives each
selected slot a uniform four-bit value, including zero. This differs from
conditioning on `j` nonzero slots. The definitions below keep the selected
set and its labels in one finite probability experiment. The point-mass
identity is the local input for the producer's binomial thinning formula.
-/

noncomputable section

namespace Spin.Turbo

open Finset
open Spin.UniformMixer

abbrev SelectedSlotSet22 (j : Fin 17) :=
  {S : Finset (Fin 16) // S.card = j.val}

noncomputable instance (j : Fin 17) : Fintype (SelectedSlotSet22 j) :=
  Fintype.ofFinite _

noncomputable instance (j : Fin 17) : Nonempty (SelectedSlotSet22 j) := by
  classical
  have hj : j.val ≤ (Finset.univ : Finset (Fin 16)).card := by
    simp only [Finset.card_univ, Fintype.card_fin]
    omega
  obtain ⟨S, -, hS⟩ :=
    (Finset.univ : Finset (Fin 16)).exists_subset_card_eq hj
  exact ⟨⟨S, hS⟩⟩

def selectedLabelWord22 (S : Finset (Fin 16)) (a : S → Packet4) :
    PacketStep16 :=
  fun i => if h : i ∈ S then a ⟨i, h⟩ else 0

def selectedGiven22 (S : Finset (Fin 16)) : FinPMF PacketStep16 :=
  (FinPMF.uniform (S → Packet4)).map (selectedLabelWord22 S)

def selectedStepLaw22 (j : Fin 17) : FinPMF PacketStep16 :=
  (FinPMF.uniform (SelectedSlotSet22 j)).bind
    (fun S => selectedGiven22 S.val)

theorem selectedLabelWord22_support_subset
    (S : Finset (Fin 16)) (a : S → Packet4) :
    packetStepSupport22 (selectedLabelWord22 S a) ⊆ S := by
  intro i hi
  by_contra hnot
  have hnonzero := (Finset.mem_filter.mp hi).2
  apply hnonzero
  simp [selectedLabelWord22, hnot]

theorem selectedLabelWord22_injective (S : Finset (Fin 16)) :
    Function.Injective (selectedLabelWord22 S) := by
  intro a b hab
  funext i
  have h := congrFun hab i.val
  simpa [selectedLabelWord22, i.property] using h

theorem selectedLabelWord22_exists_iff
    (S : Finset (Fin 16)) (x : PacketStep16) :
    (∃ a : S → Packet4, selectedLabelWord22 S a = x) ↔
      packetStepSupport22 x ⊆ S := by
  constructor
  · rintro ⟨a, rfl⟩
    exact selectedLabelWord22_support_subset S a
  · intro hx
    refine ⟨fun i => x i.val, ?_⟩
    funext i
    by_cases hi : i ∈ S
    · simp [selectedLabelWord22, hi]
    · have hz : x i = 0 := by
        by_contra hn
        exact hi (hx (Finset.mem_filter.mpr ⟨Finset.mem_univ _, hn⟩))
      simp [selectedLabelWord22, hi, hz]

theorem selectedGiven22_apply (S : Finset (Fin 16))
    (x : PacketStep16) :
    (selectedGiven22 S).p x =
      if packetStepSupport22 x ⊆ S then
        1 / (16 : ℝ) ^ S.card else 0 := by
  classical
  rw [selectedGiven22, FinPMF.map_p]
  by_cases hx : packetStepSupport22 x ⊆ S
  · rw [if_pos hx]
    obtain ⟨a, ha⟩ := (selectedLabelWord22_exists_iff S x).mpr hx
    have he (b : S → Packet4) :
        selectedLabelWord22 S b = x ↔ b = a := by
      constructor
      · intro hb
        exact selectedLabelWord22_injective S (hb.trans ha.symm)
      · rintro rfl
        exact ha
    simp only [he, Finset.sum_ite_eq', Finset.mem_univ, if_pos]
    change 1 / (Fintype.card (S → Packet4) : ℝ) = _
    simp [Fintype.card_fun, Fintype.card_fin]
  · rw [if_neg hx]
    apply Finset.sum_eq_zero
    intro a _
    have ha : selectedLabelWord22 S a ≠ x := by
      intro he
      exact hx ((selectedLabelWord22_exists_iff S x).mp ⟨a, he⟩)
    simp [ha]

theorem selectedSlotSet22_card (j : Fin 17) :
    Fintype.card (SelectedSlotSet22 j) = Nat.choose 16 j.val := by
  classical
  rw [Fintype.card_subtype]
  have he : (Finset.univ.filter (fun S : Finset (Fin 16) =>
      S.card = j.val)) =
      (Finset.univ : Finset (Fin 16)).powersetCard j.val := by
    ext S
    simp [Finset.mem_powersetCard]
  rw [he, Finset.card_powersetCard, Finset.card_univ, Fintype.card_fin]

theorem selectedSlotSet22_supersets (j : Fin 17)
    (U : Finset (Fin 16)) (hU : U.card ≤ j.val) :
    (Finset.univ.filter (fun S : SelectedSlotSet22 j => U ⊆ S.val)).card =
      Nat.choose (16 - U.card) (j.val - U.card) := by
  classical
  let E : {S : SelectedSlotSet22 j // U ⊆ S.val} ≃
      {S : Finset (Fin 16) // S ∈
        ((Finset.univ : Finset (Fin 16)).powersetCard j.val).filter
          (fun S => U ⊆ S)} := {
    toFun := fun S => ⟨S.val.val, by
      simp [Finset.mem_powersetCard, S.val.property, S.property]⟩
    invFun := fun S => ⟨⟨S.val, by
      have h := (Finset.mem_filter.mp S.property).1
      exact (Finset.mem_powersetCard.mp h).2⟩,
        (Finset.mem_filter.mp S.property).2⟩
    left_inv := by intro S; apply Subtype.ext; apply Subtype.ext; rfl
    right_inv := by intro S; apply Subtype.ext; rfl
  }
  rw [← Fintype.card_subtype, Fintype.card_congr E, Fintype.card_subtype]
  have he : (Finset.univ.filter (fun S : Finset (Fin 16) =>
      S ∈ ((Finset.univ : Finset (Fin 16)).powersetCard j.val).filter
        (fun S => U ⊆ S))) =
      ((Finset.univ : Finset (Fin 16)).powersetCard j.val).filter
        (fun S => U ⊆ S) := by
    ext S
    simp
  rw [he]
  exact Finset.card_filter_powersetCard_subset U univ j.val
    (Finset.subset_univ U) hU

/-- The exact point mass of the producer's zero-inclusive slot experiment.
The numerator counts selected sets that contain every nonzero slot of `x`. -/
theorem selectedStepLaw22_apply (j : Fin 17) (x : PacketStep16) :
    (selectedStepLaw22 j).p x =
      if (packetStepSupport22 x).card ≤ j.val then
        (Nat.choose (16 - (packetStepSupport22 x).card)
          (j.val - (packetStepSupport22 x).card) : ℝ) /
          ((Nat.choose 16 j.val : ℝ) * (16 : ℝ) ^ j.val)
      else 0 := by
  classical
  rw [selectedStepLaw22, FinPMF.bind_p]
  simp only [selectedGiven22_apply]
  have hterm (S : SelectedSlotSet22 j) :
      (FinPMF.uniform (SelectedSlotSet22 j)).p S *
          (if packetStepSupport22 x ⊆ S.val then
            1 / (16 : ℝ) ^ S.val.card else 0) =
        if packetStepSupport22 x ⊆ S.val then
          1 / ((Nat.choose 16 j.val : ℝ) * (16 : ℝ) ^ j.val) else 0 := by
    rw [show S.val.card = j.val from S.property]
    change (1 / (Fintype.card (SelectedSlotSet22 j) : ℝ)) * _ = _
    rw [selectedSlotSet22_card]
    split_ifs <;> simp <;> ring
  simp only [hterm]
  rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul]
  by_cases hU : (packetStepSupport22 x).card ≤ j.val
  · rw [if_pos hU, selectedSlotSet22_supersets j _ hU]
    push_cast
    ring
  · rw [if_neg hU]
    have he : (Finset.univ.filter (fun S : SelectedSlotSet22 j =>
        packetStepSupport22 x ⊆ S.val)) = ∅ := by
      apply Finset.filter_eq_empty_iff.mpr
      intro S _ hs
      exact hU (by have hh := Finset.card_le_card hs; omega)
    rw [he]
    simp

/-- The occupancy of a complete 16-slot packet word. -/
def packetStepOccupancy22 (x : PacketStep16) : Fin 17 :=
  ⟨(packetStepSupport22 x).card, by
    have h : (packetStepSupport22 x).card ≤ 16 := by
      exact le_trans (Finset.card_le_card
        (Finset.filter_subset (fun i : Fin 16 => x i ≠ 0) Finset.univ))
        (by simp)
    omega⟩

/-- Probability that `j` independent uniform four-bit values contain
exactly `k` nonzero values. -/
def thinningWeight22 (j k : Fin 17) : ℝ :=
  (Nat.choose j.val k.val : ℝ) * (15 : ℝ) ^ k.val /
    (16 : ℝ) ^ j.val

/-- Pointwise version of the producer's `15/16` thinning rule. The complete
packet word, including its nonzero labels, is retained on both sides. -/
theorem selectedStepLaw22_apply_thinned (j : Fin 17)
    (x : PacketStep16) :
    (selectedStepLaw22 j).p x =
      (thinningWeight22 j (packetStepOccupancy22 x)) *
        (packetStepLaw22 (packetStepOccupancy22 x)).p x := by
  classical
  let k := packetStepOccupancy22 x
  change (selectedStepLaw22 j).p x =
    thinningWeight22 j k * (packetStepLaw22 k).p x
  have hk : (packetStepSupport22 x).card = k.val := rfl
  rw [selectedStepLaw22_apply, packetStepLaw22_apply]
  simp only [hk, if_pos rfl]
  by_cases hkj : k.val ≤ j.val
  · rw [if_pos hkj, thinningWeight22]
    have hj : j.val ≤ 16 := by omega
    have hc := Spin.Structured.Placement.choose_ratio
      (n := 16) hkj hj
    have hn : (Nat.choose 16 k.val : ℝ) ≠ 0 := by
      exact_mod_cast Nat.choose_ne_zero (by omega : k.val ≤ 16)
    have hp : (15 : ℝ) ^ k.val ≠ 0 := by positivity
    have hq : (16 : ℝ) ^ j.val ≠ 0 := by positivity
    calc
      _ = ((Nat.choose (16 - k.val) (j.val - k.val) : ℝ) /
          (Nat.choose 16 j.val : ℝ)) / (16 : ℝ) ^ j.val := by
            simp [div_eq_mul_inv, mul_inv_rev]
            ring
      _ = ((Nat.choose j.val k.val : ℝ) /
          (Nat.choose 16 k.val : ℝ)) / (16 : ℝ) ^ j.val := by
            rw [hc]
      _ = (Nat.choose j.val k.val : ℝ) * (15 : ℝ) ^ k.val /
          (16 : ℝ) ^ j.val *
          (1 / ((Nat.choose 16 k.val : ℝ) * (15 : ℝ) ^ k.val)) := by
            field_simp
  · rw [if_neg hkj, thinningWeight22]
    have hz : Nat.choose j.val k.val = 0 :=
      Nat.choose_eq_zero_of_lt (by omega)
    simp [hz]

/-- Exact finite mixture used by the physical-operator producer. It is a
pointwise identity of complete packet-word laws, not only of occupancies. -/
theorem selectedStepLaw22_mixture (j : Fin 17)
    (x : PacketStep16) :
    (selectedStepLaw22 j).p x =
      ∑ k : Fin 17, thinningWeight22 j k * (packetStepLaw22 k).p x := by
  classical
  rw [selectedStepLaw22_apply_thinned]
  let k₀ := packetStepOccupancy22 x
  rw [Finset.sum_eq_single k₀]
  · intro k _ hne
    have hk : (packetStepSupport22 x).card ≠ k.val := by
      intro h
      apply hne
      exact Fin.ext (by simpa [k₀, packetStepOccupancy22] using h.symm)
    rw [packetStepLaw22_apply, if_neg hk]
    simp
  · intro h
    exact (h (Finset.mem_univ k₀)).elim

/-- Full-state weighted kernel for zero-inclusive selected slots. -/
def turboSelectedKernel22 (j : Fin 17) (z : ℝ) :
    State22 → State22 → ℝ :=
  fun q r => (selectedStepLaw22 j).expect (fun x =>
    stepKernel z (packetStepInput22 x) q r)

/-- The producer's thinning is exact for the full `2^22`-state kernel.
The later 11-coordinate compression still needs a separate domination
theorem. -/
theorem turboSelectedKernel22_mixture (j : Fin 17) (z : ℝ)
    (q r : State22) :
    turboSelectedKernel22 j z q r =
      ∑ k : Fin 17, thinningWeight22 j k *
        turboPhysicalKernel22 k z q r := by
  classical
  unfold turboSelectedKernel22 turboPhysicalKernel22 FinPMF.expect
  simp_rw [selectedStepLaw22_mixture]
  simp only [Finset.sum_mul]
  rw [Finset.sum_comm]
  apply Finset.sum_congr rfl
  intro k _
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro x _
  ring

#print axioms selectedGiven22_apply
#print axioms selectedStepLaw22_apply
#print axioms selectedStepLaw22_apply_thinned
#print axioms selectedStepLaw22_mixture
#print axioms turboSelectedKernel22_mixture

end Spin.Turbo
