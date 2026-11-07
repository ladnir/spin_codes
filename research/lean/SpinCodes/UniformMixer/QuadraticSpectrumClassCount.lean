import SpinCodes.UniformMixer.QuadraticProductionTriple
import SpinCodes.UniformMixer.QuadraticRankWalshClasses
import SpinCodes.UniformMixer.QuadraticRadicalPartition

/-!
# Counting weights within one quadratic-radical class

For each fixed quadratic part and linear coefficient, the constant bit
complements the output word. A nonzero Walsh coefficient therefore gives
one word at each of the two complementary off-center weights.
-/

namespace Spin.UniformMixer

def quadraticWordWeight (b : ZMod 2) (a : State6)
    (c : Fin 15 → ZMod 2) : Nat :=
  boolWeight6 (fun x => b + dot6Linear x a + quadraticPart c x)

theorem walsh_quadraticWordWeight (b : ZMod 2) (a : State6)
    (c : Fin 15 → ZMod 2) :
    binarySign b * walsh6 (quadraticPart c) a =
      64 - 2 * (quadraticWordWeight b a c : ℤ) := by
  rw [← walsh6_constant_shift (quadraticPart c) b a]
  rw [walsh6_weight]
  have hfunc :
      (fun x => quadraticPart c x + b + dot6Linear x a) =
        (fun x => b + dot6Linear x a + quadraticPart c x) := by
    funext x
    ring
  rw [hfunc]
  rfl

theorem quadraticWordWeight_complement (a : State6)
    (c : Fin 15 → ZMod 2) :
    quadraticWordWeight 0 a c + quadraticWordWeight 1 a c = 64 := by
  have h0 := walsh_quadraticWordWeight 0 a c
  have h1 := walsh_quadraticWordWeight 1 a c
  rw [binarySign_zero, one_mul] at h0
  rw [binarySign_one] at h1
  omega

theorem quadraticWordWeight_active_pair (a : State6)
    (c : Fin 15 → ZMod 2) (amplitude low high : Nat)
    (hvalue : walsh6 (quadraticPart c) a = (amplitude : ℤ) ∨
      walsh6 (quadraticPart c) a = -(amplitude : ℤ))
    (hlow : 2 * low + amplitude = 64)
    (hhigh : 2 * high = 64 + amplitude) :
    (quadraticWordWeight 0 a c = low ∧
      quadraticWordWeight 1 a c = high) ∨
    (quadraticWordWeight 0 a c = high ∧
      quadraticWordWeight 1 a c = low) := by
  have h0 := walsh_quadraticWordWeight 0 a c
  have h1 := walsh_quadraticWordWeight 1 a c
  rw [binarySign_zero, one_mul] at h0
  rw [binarySign_one] at h1
  rcases hvalue with hp | hn
  · left
    constructor <;> omega
  · right
    constructor <;> omega

theorem quadraticWordWeight_inactive (a : State6)
    (c : Fin 15 → ZMod 2) (ha : a ∉ activeWalsh c) :
    quadraticWordWeight 0 a c = 32 ∧
      quadraticWordWeight 1 a c = 32 := by
  have hz : walsh6 (quadraticPart c) a = 0 := by
    simpa [activeWalsh] using ha
  have h0 := walsh_quadraticWordWeight 0 a c
  have h1 := walsh_quadraticWordWeight 1 a c
  rw [binarySign_zero, one_mul, hz] at h0
  rw [binarySign_one, hz] at h1
  constructor <;> omega

theorem quadraticWordWeight_low_pair_count (a : State6)
    (c : Fin 15 → ZMod 2) (amplitude low high : Nat)
    (hpos : 0 < amplitude)
    (hlow : 2 * low + amplitude = 64)
    (hhigh : 2 * high = 64 + amplitude)
    (hvalue : a ∈ activeWalsh c →
      walsh6 (quadraticPart c) a = (amplitude : ℤ) ∨
        walsh6 (quadraticPart c) a = -(amplitude : ℤ)) :
    (if quadraticWordWeight 0 a c = low then 1 else 0) +
      (if quadraticWordWeight 1 a c = low then 1 else 0) =
        if a ∈ activeWalsh c then 1 else 0 := by
  have hne : high ≠ low := by omega
  have hlo32 : (32 : Nat) ≠ low := by omega
  by_cases ha : a ∈ activeWalsh c
  · rcases quadraticWordWeight_active_pair a c amplitude low high
      (hvalue ha) hlow hhigh with ⟨h0, h1⟩ | ⟨h0, h1⟩
    · simp [ha, h0, h1, hne]
    · simp [ha, h0, h1, hne]
  · obtain ⟨h0, h1⟩ := quadraticWordWeight_inactive a c ha
    simp [ha, h0, h1, hlo32]

noncomputable def quadraticPairLowCount (c : Fin 15 → ZMod 2)
    (low : Nat) : Nat :=
  ∑ a : State6,
    ((if quadraticWordWeight 0 a c = low then 1 else 0) +
      (if quadraticWordWeight 1 a c = low then 1 else 0))

theorem quadraticPairLowCount_eq_active (c : Fin 15 → ZMod 2)
    (amplitude low high : Nat)
    (hpos : 0 < amplitude)
    (hlow : 2 * low + amplitude = 64)
    (hhigh : 2 * high = 64 + amplitude)
    (hvalue : ∀ a ∈ activeWalsh c,
      walsh6 (quadraticPart c) a = (amplitude : ℤ) ∨
        walsh6 (quadraticPart c) a = -(amplitude : ℤ)) :
    quadraticPairLowCount c low = (activeWalsh c).card := by
  unfold quadraticPairLowCount
  calc
    _ = ∑ a : State6, if a ∈ activeWalsh c then 1 else 0 := by
      apply Finset.sum_congr rfl
      intro a _
      exact quadraticWordWeight_low_pair_count a c amplitude low high
        hpos hlow hhigh (hvalue a)
    _ = (activeWalsh c).card := by
      exact (Finset.card_eq_sum_ite (Finset.subset_univ (activeWalsh c))).symm

theorem quadraticPairLowCount_radical_64 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 64) :
    quadraticPairLowCount c 0 = 1 := by
  rw [quadraticPairLowCount_eq_active c 64 0 64 (by omega)
    (by omega) (by omega) (fun a ha => walsh_value_radical_64 c hr a ha)]
  exact activeWalsh_card_radical_64 c hr

theorem quadraticPairLowCount_radical_16 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 16) :
    quadraticPairLowCount c 16 = 4 := by
  rw [quadraticPairLowCount_eq_active c 32 16 48 (by omega)
    (by omega) (by omega) (fun a ha => walsh_value_radical_16 c hr a ha)]
  exact activeWalsh_card_radical_16 c hr

theorem quadraticPairLowCount_radical_4 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 4) :
    quadraticPairLowCount c 24 = 16 := by
  rw [quadraticPairLowCount_eq_active c 16 24 40 (by omega)
    (by omega) (by omega) (fun a ha => walsh_value_radical_4 c hr a ha)]
  exact activeWalsh_card_radical_4 c hr

theorem quadraticPairLowCount_radical_1 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 1) :
    quadraticPairLowCount c 28 = 64 := by
  rw [quadraticPairLowCount_eq_active c 8 28 36 (by omega)
    (by omega) (by omega) (fun a ha => walsh_value_radical_1 c hr a ha)]
  exact activeWalsh_card_radical_1 c hr

theorem quadraticWordWeight_class (b : ZMod 2) (a : State6)
    (c : Fin 15 → ZMod 2) (amplitude low high : Nat)
    (hlow : 2 * low + amplitude = 64)
    (hhigh : 2 * high = 64 + amplitude)
    (hvalue : a ∈ activeWalsh c →
      walsh6 (quadraticPart c) a = (amplitude : ℤ) ∨
        walsh6 (quadraticPart c) a = -(amplitude : ℤ)) :
    quadraticWordWeight b a c = low ∨
      quadraticWordWeight b a c = 32 ∨
      quadraticWordWeight b a c = high := by
  fin_cases b
  · by_cases ha : a ∈ activeWalsh c
    · rcases quadraticWordWeight_active_pair a c amplitude low high
        (hvalue ha) hlow hhigh with h | h
      · exact Or.inl h.1
      · exact Or.inr (Or.inr h.1)
    · exact Or.inr (Or.inl (quadraticWordWeight_inactive a c ha).1)
  · by_cases ha : a ∈ activeWalsh c
    · rcases quadraticWordWeight_active_pair a c amplitude low high
        (hvalue ha) hlow hhigh with h | h
      · exact Or.inr (Or.inr h.2)
      · exact Or.inl h.2
    · exact Or.inr (Or.inl (quadraticWordWeight_inactive a c ha).2)

theorem quadraticPairLowCount_outside (c : Fin 15 → ZMod 2)
    (amplitude low high w : Nat)
    (hlow : 2 * low + amplitude = 64)
    (hhigh : 2 * high = 64 + amplitude)
    (hvalue : ∀ a ∈ activeWalsh c,
      walsh6 (quadraticPart c) a = (amplitude : ℤ) ∨
        walsh6 (quadraticPart c) a = -(amplitude : ℤ))
    (hwlow : w ≠ low) (hwcenter : w ≠ 32) (hwhigh : w ≠ high) :
    quadraticPairLowCount c w = 0 := by
  unfold quadraticPairLowCount
  apply Finset.sum_eq_zero
  intro a _
  have h0 : quadraticWordWeight 0 a c ≠ w := by
    rcases quadraticWordWeight_class 0 a c amplitude low high hlow hhigh
      (hvalue a) with h | h | h <;> omega
  have h1 : quadraticWordWeight 1 a c ≠ w := by
    rcases quadraticWordWeight_class 1 a c amplitude low high hlow hhigh
      (hvalue a) with h | h | h <;> omega
  simp [h0, h1]

theorem quadraticPairLowCount_complement (c : Fin 15 → ZMod 2)
    (w : Nat) (hw : w ≤ 64) :
    quadraticPairLowCount c w = quadraticPairLowCount c (64 - w) := by
  unfold quadraticPairLowCount
  apply Finset.sum_congr rfl
  intro a _
  have hc := quadraticWordWeight_complement a c
  have h0 : quadraticWordWeight 0 a c = w ↔
      quadraticWordWeight 1 a c = 64 - w := by omega
  have h1 : quadraticWordWeight 1 a c = w ↔
      quadraticWordWeight 0 a c = 64 - w := by omega
  have hi0 : (if quadraticWordWeight 0 a c = w then 1 else 0) =
      (if quadraticWordWeight 1 a c = 64 - w then 1 else 0) := by
    by_cases h : quadraticWordWeight 0 a c = w
    · rw [if_pos h, if_pos (h0.mp h)]
    · rw [if_neg h, if_neg (fun hx => h (h0.mpr hx))]
  have hi1 : (if quadraticWordWeight 1 a c = w then 1 else 0) =
      (if quadraticWordWeight 0 a c = 64 - w then 1 else 0) := by
    by_cases h : quadraticWordWeight 1 a c = w
    · rw [if_pos h, if_pos (h1.mp h)]
    · rw [if_neg h, if_neg (fun hx => h (h1.mpr hx))]
  rw [hi0, hi1]
  exact Nat.add_comm _ _

theorem quadraticPairCenterCount_eq_inactive (c : Fin 15 → ZMod 2)
    (amplitude low high : Nat)
    (hpos : 0 < amplitude)
    (hlow : 2 * low + amplitude = 64)
    (hhigh : 2 * high = 64 + amplitude)
    (hvalue : ∀ a ∈ activeWalsh c,
      walsh6 (quadraticPart c) a = (amplitude : ℤ) ∨
        walsh6 (quadraticPart c) a = -(amplitude : ℤ)) :
    quadraticPairLowCount c 32 =
      2 * (64 - (activeWalsh c).card) := by
  have hlow32 : low ≠ 32 := by omega
  have hhigh32 : high ≠ 32 := by omega
  have hpoint (a : State6) :
      (if quadraticWordWeight 0 a c = 32 then 1 else 0) +
        (if quadraticWordWeight 1 a c = 32 then 1 else 0) =
          if a ∈ activeWalsh c then 0 else 2 := by
    by_cases ha : a ∈ activeWalsh c
    · rcases quadraticWordWeight_active_pair a c amplitude low high
        (hvalue a ha) hlow hhigh with ⟨h0, h1⟩ | ⟨h0, h1⟩
      · simp [ha, h0, h1, hlow32, hhigh32]
      · simp [ha, h0, h1, hlow32, hhigh32]
    · obtain ⟨h0, h1⟩ := quadraticWordWeight_inactive a c ha
      simp [ha, h0, h1]
  unfold quadraticPairLowCount
  simp_rw [hpoint]
  rw [Finset.sum_ite]
  simp only [Finset.sum_const_zero, zero_add, Finset.sum_const,
    smul_eq_mul]
  rw [Finset.filter_notMem_eq_sdiff,
    Finset.card_sdiff_of_subset (Finset.subset_univ (activeWalsh c))]
  have hcard : (Finset.univ : Finset State6).card = 64 := by
    norm_num [Fintype.card_fun]
  rw [hcard]
  omega

theorem quadraticPairCenterCount_radical_64 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 64) :
    quadraticPairLowCount c 32 = 126 := by
  rw [quadraticPairCenterCount_eq_inactive c 64 0 64
    (by omega) (by omega) (by omega)
    (fun a ha => walsh_value_radical_64 c hr a ha),
    activeWalsh_card_radical_64 c hr]

theorem quadraticPairCenterCount_radical_16 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 16) :
    quadraticPairLowCount c 32 = 120 := by
  rw [quadraticPairCenterCount_eq_inactive c 32 16 48
    (by omega) (by omega) (by omega)
    (fun a ha => walsh_value_radical_16 c hr a ha),
    activeWalsh_card_radical_16 c hr]

theorem quadraticPairCenterCount_radical_4 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 4) :
    quadraticPairLowCount c 32 = 96 := by
  rw [quadraticPairCenterCount_eq_inactive c 16 24 40
    (by omega) (by omega) (by omega)
    (fun a ha => walsh_value_radical_4 c hr a ha),
    activeWalsh_card_radical_4 c hr]

theorem quadraticPairCenterCount_radical_1 (c : Fin 15 → ZMod 2)
    (hr : (quadraticRadical c).card = 1) :
    quadraticPairLowCount c 32 = 0 := by
  rw [quadraticPairCenterCount_eq_inactive c 8 28 36
    (by omega) (by omega) (by omega)
    (fun a ha => walsh_value_radical_1 c hr a ha),
    activeWalsh_card_radical_1 c hr]

theorem quadraticPairCenterCount_eq_class_expression
    (c : Fin 15 → ZMod 2) :
    quadraticPairLowCount c 32 =
      (if (quadraticRadical c).card = 64 then 126 else 0) +
      (if (quadraticRadical c).card = 16 then 120 else 0) +
      (if (quadraticRadical c).card = 4 then 96 else 0) := by
  rcases quadraticRadical_four_classes c with h64 | h16 | h4 | h1
  · rw [quadraticPairCenterCount_radical_64 c h64]
    have hn16 : (quadraticRadical c).card ≠ 16 := by omega
    have hn4 : (quadraticRadical c).card ≠ 4 := by omega
    simp only [if_pos h64, if_neg hn16, if_neg hn4, add_zero]
  · rw [quadraticPairCenterCount_radical_16 c h16]
    have hn64 : (quadraticRadical c).card ≠ 64 := by omega
    have hn4 : (quadraticRadical c).card ≠ 4 := by omega
    simp only [if_neg hn64, if_pos h16, if_neg hn4, zero_add, add_zero]
  · rw [quadraticPairCenterCount_radical_4 c h4]
    have hn64 : (quadraticRadical c).card ≠ 64 := by omega
    have hn16 : (quadraticRadical c).card ≠ 16 := by omega
    simp only [if_neg hn64, if_neg hn16, if_pos h4, zero_add]
  · rw [quadraticPairCenterCount_radical_1 c h1]
    have hn64 : (quadraticRadical c).card ≠ 64 := by omega
    have hn16 : (quadraticRadical c).card ≠ 16 := by omega
    have hn4 : (quadraticRadical c).card ≠ 4 := by omega
    simp only [if_neg hn64, if_neg hn16, if_neg hn4, zero_add]

end Spin.UniformMixer
