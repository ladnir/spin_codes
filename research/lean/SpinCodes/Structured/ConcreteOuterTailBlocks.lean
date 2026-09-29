import SpinCodes.Structured.ConcreteOuterEnvelopeShiftSupport
import SpinCodes.Structured.ConcreteOuterTailSparseMoment

set_option maxHeartbeats 1000000
noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset Spin.Numeric

def golayActive {k : ℕ} (x : LocalMessage k) : Finset (Fin k) := univ.filter (fun i => x i ≠ 0)
def golayOccupation {k : ℕ} (x : LocalMessage k) : ℕ := (golayActive x).card

lemma golayActive_zero_iff {k : ℕ} (x : LocalMessage k) : golayActive x = ∅ ↔ x=0 := by
  simp only [golayActive, filter_eq_empty_iff, mem_univ, true_implies, not_not]
  exact ⟨fun h => funext h, fun h i => by rw [h]; rfl⟩

lemma golayOccupation_pos {k : ℕ} {x : LocalMessage k} (hx : x≠0) : 0<golayOccupation x := by
  rw [golayOccupation, card_pos, Finset.nonempty_iff_ne_empty]
  exact fun h => hx ((golayActive_zero_iff x).mp h)

lemma golay_active_message_count {k : ℕ} (S : Finset (Fin k)) :
    (univ.filter (fun x : LocalMessage k => golayActive x=S)).card = 4095^S.card := by
  have he : univ.filter (fun x : LocalMessage k => golayActive x=S) =
      Fintype.piFinset (fun i : Fin k => if i∈S then univ.erase (0:Fin 12→Bool) else {0}) := by
    ext x
    simp only [mem_filter, mem_univ, true_and, Fintype.mem_piFinset]
    constructor
    · intro hx i
      have hh := Finset.ext_iff.mp hx i
      simp only [golayActive, mem_filter, mem_univ, true_and] at hh
      by_cases hi : i∈S
      · simp [hi, hh.mpr hi]
      · have hz : x i=0 := by simpa only [hi, iff_false, not_not] using hh
        simp [hi, hz]
    · intro hx
      ext i
      have hh := hx i
      by_cases hi : i∈S <;> simp_all [golayActive]
  rw [he, Fintype.card_piFinset]
  have hc : (univ.erase (0:Fin 12→Bool)).card = 4095 := by simp
  have hcard (i : Fin k) : (if i∈S then univ.erase (0:Fin 12→Bool) else {0}).card =
      if i∈S then 4095 else 1 := by split_ifs <;> simp [hc]
  simp_rw [hcard]
  rw [Finset.prod_ite]
  simp

/-- Exact count of messages with j active Golay blocks. -/
theorem golay_occupation_count (k j : ℕ) :
    (univ.filter (fun x : LocalMessage k => golayOccupation x=j)).card = k.choose j*4095^j := by
  have hh := sum_card_fiberwise_eq_card_filter (univ : Finset (LocalMessage k))
    ((univ : Finset (Fin k)).powersetCard j) golayActive
  have hfilter : univ.filter (fun x : LocalMessage k => golayActive x∈(univ : Finset (Fin k)).powersetCard j) =
      univ.filter (fun x : LocalMessage k => golayOccupation x=j) := by
    ext x
    simp only [mem_filter, mem_univ, true_and, mem_powersetCard, subset_univ]
    rfl
  rw [hfilter] at hh
  rw [←hh]
  have he : ∀ S ∈ (univ : Finset (Fin k)).powersetCard j,
      (univ.filter (fun x : LocalMessage k => golayActive x=S)).card = 4095^j := by
    intro S hS
    rw [golay_active_message_count, (mem_powersetCard.mp hS).2]
  rw [sum_const_nat he, card_powersetCard]
  simp

lemma golay_inputWeight_bounds {k : ℕ} (x : LocalMessage k) :
    8*golayOccupation x ≤ wtF (outerWord Golay.enc x) ∧
      wtF (outerWord Golay.enc x) ≤ 24*golayOccupation x := by
  have he : wtF (outerWord Golay.enc x)=∑ i ∈ golayActive x, Golay.W (x i) := by
    rw [wtF_outerWord]
    change (∑ i, Golay.W (x i))=_
    symm
    apply sum_subset (filter_subset _ _)
    intro i _ hi
    have hx : x i=0 := by simpa [golayActive] using hi
    simp [hx, Golay.W_zero]
  rw [he]
  constructor
  · have hh : (∑ _i ∈ golayActive x, 8) ≤ ∑ i ∈ golayActive x, Golay.W (x i) := by
      apply sum_le_sum
      intro i hi
      apply Golay.min_distance
      intro hz
      exact (mem_filter.mp hi).2 ((Golay.W_eq_zero_iff _).mp hz)
    simpa [golayOccupation, mul_comm] using hh
  · have hh : (∑ i ∈ golayActive x, Golay.W (x i)) ≤ ∑ _i ∈ golayActive x, 24 :=
      sum_le_sum (fun i _ => Golay.W_le (x i))
    simpa [golayOccupation, mul_comm] using hh

lemma golay_inputWeight_even {k : ℕ} (x : LocalMessage k) : Even (wtF (outerWord Golay.enc x)) := by
  rw [wtF_outerWord, even_iff_two_dvd]
  exact Finset.dvd_sum (fun i _ => even_iff_two_dvd.mp (golay_weight_even (x i)))

/-- The actual first-accumulator moment bound for a nonzero sparse constituent message. -/
theorem golay_sparse_moment {k : ℕ} (x : LocalMessage k) (hx : x≠0)
    (hsp : 1000*golayOccupation x ≤ k*24) :
    (∑ h ∈ range (k*24+1), Pt (k*24) (wtF (outerWord Golay.enc x)) h*zeta^h) ≤
      (4*kappa*(golayOccupation x:ℝ)/((k*24:ℕ)*(1-8*eta)))^(4*golayOccupation x) := by
  obtain ⟨ℓ,hℓ⟩ := golay_inputWeight_even x
  obtain ⟨hlo,hhi⟩ := golay_inputWeight_bounds x
  have he : wtF (outerWord Golay.enc x)=2*ℓ := by omega
  rw [he]
  exact sparse_accumulator_moment_le (golayOccupation_pos hx) (by omega) (by omega) hsp

end Spin.Structured.ConcreteOuter



