import SpinCodes.Structured.ConcreteNativeFirstMoment
import SpinCodes.Selection

noncomputable section
attribute [local instance] Classical.propDecidable
namespace Spin.Structured.ConcreteNativeFamily
open Finset ConcreteOuter ConcreteRoute ConcreteBinaryEncoding

/-- The exact paper window [13/125,112/125], expressed with integer inequalities. -/
def weightWindow (b : ℕ) : Finset ℕ :=
  (range (b+1)).filter (fun w => 13*b ≤ 125*w ∧ 125*w ≤ 112*b)

def constituentGood (k : ℕ) (seed : ConcreteOuter.Seed k) : Prop :=
  Good (ConcreteOuter.seedLaw k) ConcreteOuter.spectrum (k*24) (weightWindow (k*24)) seed

def nativeGood (m : ℕ) : NativeSeed m → Prop := constituentGood (nativeBlocks m)

theorem good_weight_mem {k : ℕ} {seed : ConcreteOuter.Seed k} (hg : constituentGood k seed)
    {x : LocalMessage k} (hx : x ≠ 0) : wtF (encode seed x) ∈ weightWindow (k*24) := by
  have hw : wtF (encode seed x) ≤ k*24 := by
    rw [← support_card]
    exact Routing.card_le_width _
  have hp : 0 < spectrum seed (wtF (encode seed x)) := by
    apply Finset.card_pos.mpr
    exact ⟨x, by simp [hx]⟩
  by_contra hn
  have hz := hg.1 (wtF (encode seed x)) (by
    simp only [badWeights, mem_sdiff, mem_range]
    exact ⟨Nat.lt_succ_of_le hw, hn⟩)
  omega

theorem good_weight_upper {k : ℕ} {seed : ConcreteOuter.Seed k} (hg : constituentGood k seed)
    (x : LocalMessage k) : 125 * wtF (encode seed x) ≤ 112 * (k*24) := by
  by_cases hx : x = 0
  · subst x
    rw [encode_zero, (wtF_eq_zero _).mpr rfl]
    omega
  · exact (mem_filter.mp (good_weight_mem hg hx)).2.2

theorem native_good_row_upper (m : ℕ) {seed : NativeSeed m} (hg : nativeGood m seed)
    (message : Fin (Nsched m / 2) → ZMod 2) (i : Fin (Lsched m)) :
    125 * (rows m seed message i).card ≤ 112 * bsched m := by
  have h := good_weight_upper hg
    (nativeMessageEquiv m ((wordEquiv (Nsched m / 2)).symm message) i)
  simpa only [rows, nativeRows, Finset.card_map, rowSupports, support_card, native_width] using h

theorem totalWeight_upper_of_rows {L b : ℕ} (x : Fin L → Finset (Fin b))
    (h : ∀ i, 125 * (x i).card ≤ 112 * b) :
    125 * totalWeight x ≤ L * (112 * b) := by
  have hs := Finset.sum_le_sum (fun i (_ : i ∈ (univ : Finset (Fin L))) => h i)
  simp only [← mul_sum, sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul] at hs
  simpa only [totalWeight, Nat.cast_id, mul_assoc, mul_left_comm, mul_comm] using hs

theorem native_good_total_weight (m : ℕ) {seed : NativeSeed m} (hg : nativeGood m seed)
    {message : Fin (Nsched m / 2) → ZMod 2} (hm : message ≠ 0) :
    0 < totalWeight (rows m seed message) ∧
      totalWeight (rows m seed message) < Lsched m * bsched m := by
  have ho : 0 < activeRows (rows m seed message) := (mem_Ico.mp (occ_mem m seed message hm)).1
  unfold activeRows at ho
  obtain ⟨i, hi⟩ := Finset.card_pos.mp ho
  have hn : rows m seed message i ≠ ∅ := (mem_filter.mp hi).2
  have hp : 0 < (rows m seed message i).card := Finset.card_pos.mpr (Finset.nonempty_iff_ne_empty.mpr hn)
  have hl : (rows m seed message i).card ≤ totalWeight (rows m seed message) :=
    Finset.single_le_sum (f := fun j => (rows m seed message j).card)
      (fun _ _ => Nat.zero_le _) (mem_univ i)
  have hs' := totalWeight_upper_of_rows (rows m seed message) (native_good_row_upper m hg message)
  have hd := Nat.mul_pos (Lsched_pos m) (bsched_pos m)
  constructor <;> nlinarith

theorem native_good_density (m : ℕ) {seed : NativeSeed m} (hg : nativeGood m seed)
    {message : Fin (Nsched m / 2) → ZMod 2} (hm : message ≠ 0) :
    0 < density (rows m seed message) ∧ density (rows m seed message) < 1 := by
  have hw := native_good_total_weight m hg hm
  have hd : (0 : ℝ) < (Lsched m : ℝ) * bsched m := by
    exact_mod_cast Nat.mul_pos (Lsched_pos m) (bsched_pos m)
  unfold density
  constructor
  · exact div_pos (by exact_mod_cast hw.1) hd
  · apply (div_lt_one hd).mpr
    exact_mod_cast hw.2

theorem selected_EZ_le_matrix_sum
    (hG : ∀ m, 0 < (nativeSeedLaw m).prob (nativeGood m)) (m Q : ℕ)
    {z : ℝ} (hz0 : 0 < z) (hz1 : z ≤ 1) :
    (family nativeGood hG).EZ m Q ≤
      ((nativeSeedLaw m).condition (nativeGood m) (hG m)).expect
        (fun out => ∑ message ∈ (nonzeroMsgs (Fin (Nsched m / 2) → ZMod 2)).filter
          (fun message => occ m out message = Q),
          ConcreteRoutedFirstMoment.matrixBound (rounds m) (threshold m) z (rows m out message)) :=
  EZ_le_matrix_sum nativeGood hG m Q (fun out hg message hm => native_good_density m hg hm) hz0 hz1

end Spin.Structured.ConcreteNativeFamily


