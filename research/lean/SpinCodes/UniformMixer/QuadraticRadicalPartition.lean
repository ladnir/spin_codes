import SpinCodes.UniformMixer.QuadraticRankHistogramSum

/-!
# Exhaustion of the quadratic radical classes

The checked histogram has 32768 entries. Since the four reported radical
classes are disjoint and their counts sum to 32768, every quadratic part
belongs to exactly one of them.
-/

namespace Spin.UniformMixer

private noncomputable def radicalClass (r : Nat) : Finset (Fin 15 → ZMod 2) :=
  Finset.univ.filter (fun c => (quadraticRadical c).card = r)

private theorem radicalClass_disjoint {r s : Nat} (hrs : r ≠ s) :
    Disjoint (radicalClass r) (radicalClass s) := by
  apply Finset.disjoint_left.mpr
  intro c hc hd
  simp only [radicalClass, Finset.mem_filter, Finset.mem_univ, true_and] at hc hd
  exact hrs (hc.symm.trans hd)

theorem quadraticRadical_four_classes (c : Fin 15 → ZMod 2) :
    (quadraticRadical c).card = 64 ∨
    (quadraticRadical c).card = 16 ∨
    (quadraticRadical c).card = 4 ∨
    (quadraticRadical c).card = 1 := by
  let s64 := radicalClass 64
  let s16 := radicalClass 16
  let s4 := radicalClass 4
  let s1 := radicalClass 1
  have hd64_16 : Disjoint s64 s16 := radicalClass_disjoint (by omega)
  have hd64_4 : Disjoint s64 s4 := radicalClass_disjoint (by omega)
  have hd16_4 : Disjoint s16 s4 := radicalClass_disjoint (by omega)
  have hd64_1 : Disjoint s64 s1 := radicalClass_disjoint (by omega)
  have hd16_1 : Disjoint s16 s1 := radicalClass_disjoint (by omega)
  have hd4_1 : Disjoint s4 s1 := radicalClass_disjoint (by omega)
  have hdA : Disjoint (s64 ∪ s16) s4 :=
    Finset.disjoint_union_left.mpr ⟨hd64_4, hd16_4⟩
  have hdB : Disjoint (s64 ∪ s16 ∪ s4) s1 :=
    Finset.disjoint_union_left.mpr
      ⟨Finset.disjoint_union_left.mpr ⟨hd64_1, hd16_1⟩, hd4_1⟩
  have hcard : (s64 ∪ s16 ∪ s4 ∪ s1).card = 32768 := by
    rw [Finset.card_union_of_disjoint hdB,
      Finset.card_union_of_disjoint hdA,
      Finset.card_union_of_disjoint hd64_16]
    change (radicalClass 64).card + (radicalClass 16).card +
      (radicalClass 4).card + (radicalClass 1).card = 32768
    simp only [radicalClass]
    rw [actualQuadraticRadicalCount_64, actualQuadraticRadicalCount_16,
      actualQuadraticRadicalCount_4, actualQuadraticRadicalCount_1]
  have hall : s64 ∪ s16 ∪ s4 ∪ s1 = Finset.univ := by
    apply Finset.eq_of_subset_of_card_le (Finset.subset_univ _)
    rw [hcard]
    norm_num [Fintype.card_fun]
  have hc : c ∈ s64 ∪ s16 ∪ s4 ∪ s1 := by
    rw [hall]
    exact Finset.mem_univ c
  simp only [Finset.mem_union, s64, s16, s4, s1,
    radicalClass, Finset.mem_filter, Finset.mem_univ, true_and] at hc
  simpa only [or_assoc] using hc

end Spin.UniformMixer
