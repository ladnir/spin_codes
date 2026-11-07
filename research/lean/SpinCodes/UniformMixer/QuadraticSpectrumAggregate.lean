import SpinCodes.UniformMixer.QuadraticSpectrumClassCount
import SpinCodes.UniformMixer.QuadraticRadicalPartition

/-!
# Exact off-center spectrum of the production map

The radical histogram counts quadratic parts. The Walsh-class lemmas count
affine shifts of each part. This file multiplies those two certified counts
without enumerating production words.
-/

namespace Spin.UniformMixer

noncomputable def quadraticTotalWeightCount (w : Nat) : Nat :=
  ∑ c : Fin 15 → ZMod 2, quadraticPairLowCount c w

theorem quadraticTotalWeightCount_unique_class (w radicalSize perForm : Nat)
    (hlocal : ∀ c : Fin 15 → ZMod 2,
      quadraticPairLowCount c w =
        if (quadraticRadical c).card = radicalSize then perForm else 0) :
    quadraticTotalWeightCount w =
      (Finset.univ.filter (fun c : Fin 15 → ZMod 2 =>
        (quadraticRadical c).card = radicalSize)).card * perForm := by
  unfold quadraticTotalWeightCount
  simp_rw [hlocal]
  rw [Finset.sum_ite]
  simp only [Finset.sum_const_zero, add_zero, Finset.sum_const,
    smul_eq_mul]

theorem quadraticTotalWeightCount_16 : quadraticTotalWeightCount 16 = 2604 := by
  have hlocal (c : Fin 15 → ZMod 2) :
      quadraticPairLowCount c 16 =
        if (quadraticRadical c).card = 16 then 4 else 0 := by
    rcases quadraticRadical_four_classes c with h64 | h16 | h4 | h1
    · have hz := quadraticPairLowCount_outside c 64 0 64 16
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_64 c h64 a ha)
        (by omega) (by omega) (by omega)
      simp [h64, hz]
    · simp [h16, quadraticPairLowCount_radical_16 c h16]
    · have hz := quadraticPairLowCount_outside c 16 24 40 16
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_4 c h4 a ha)
        (by omega) (by omega) (by omega)
      simp [h4, hz]
    · have hz := quadraticPairLowCount_outside c 8 28 36 16
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_1 c h1 a ha)
        (by omega) (by omega) (by omega)
      simp [h1, hz]
  unfold quadraticTotalWeightCount
  simp_rw [hlocal]
  rw [Finset.sum_ite]
  simp only [Finset.sum_const_zero, add_zero, Finset.sum_const,
    smul_eq_mul]
  rw [actualQuadraticRadicalCount_16]

theorem quadraticTotalWeightCount_0 : quadraticTotalWeightCount 0 = 1 := by
  have hlocal (c : Fin 15 → ZMod 2) :
      quadraticPairLowCount c 0 =
        if (quadraticRadical c).card = 64 then 1 else 0 := by
    rcases quadraticRadical_four_classes c with h64 | h16 | h4 | h1
    · simp [h64, quadraticPairLowCount_radical_64 c h64]
    · have hz := quadraticPairLowCount_outside c 32 16 48 0
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_16 c h16 a ha)
        (by omega) (by omega) (by omega)
      simp [h16, hz]
    · have hz := quadraticPairLowCount_outside c 16 24 40 0
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_4 c h4 a ha)
        (by omega) (by omega) (by omega)
      simp [h4, hz]
    · have hz := quadraticPairLowCount_outside c 8 28 36 0
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_1 c h1 a ha)
        (by omega) (by omega) (by omega)
      simp [h1, hz]
  rw [quadraticTotalWeightCount_unique_class 0 64 1 hlocal,
    actualQuadraticRadicalCount_64]

theorem quadraticTotalWeightCount_24 : quadraticTotalWeightCount 24 = 291648 := by
  have hlocal (c : Fin 15 → ZMod 2) :
      quadraticPairLowCount c 24 =
        if (quadraticRadical c).card = 4 then 16 else 0 := by
    rcases quadraticRadical_four_classes c with h64 | h16 | h4 | h1
    · have hz := quadraticPairLowCount_outside c 64 0 64 24
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_64 c h64 a ha)
        (by omega) (by omega) (by omega)
      simp [h64, hz]
    · have hz := quadraticPairLowCount_outside c 32 16 48 24
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_16 c h16 a ha)
        (by omega) (by omega) (by omega)
      simp [h16, hz]
    · simp [h4, quadraticPairLowCount_radical_4 c h4]
    · have hz := quadraticPairLowCount_outside c 8 28 36 24
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_1 c h1 a ha)
        (by omega) (by omega) (by omega)
      simp [h1, hz]
  rw [quadraticTotalWeightCount_unique_class 24 4 16 hlocal,
    actualQuadraticRadicalCount_4]

theorem quadraticTotalWeightCount_28 : quadraticTotalWeightCount 28 = 888832 := by
  have hlocal (c : Fin 15 → ZMod 2) :
      quadraticPairLowCount c 28 =
        if (quadraticRadical c).card = 1 then 64 else 0 := by
    rcases quadraticRadical_four_classes c with h64 | h16 | h4 | h1
    · have hz := quadraticPairLowCount_outside c 64 0 64 28
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_64 c h64 a ha)
        (by omega) (by omega) (by omega)
      simp [h64, hz]
    · have hz := quadraticPairLowCount_outside c 32 16 48 28
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_16 c h16 a ha)
        (by omega) (by omega) (by omega)
      simp [h16, hz]
    · have hz := quadraticPairLowCount_outside c 16 24 40 28
        (by omega) (by omega)
        (fun a ha => walsh_value_radical_4 c h4 a ha)
        (by omega) (by omega) (by omega)
      simp [h4, hz]
    · simp [h1, quadraticPairLowCount_radical_1 c h1]
  rw [quadraticTotalWeightCount_unique_class 28 1 64 hlocal,
    actualQuadraticRadicalCount_1]

theorem quadraticTotalWeightCount_complement (w : Nat) (hw : w ≤ 64) :
    quadraticTotalWeightCount w = quadraticTotalWeightCount (64 - w) := by
  unfold quadraticTotalWeightCount
  apply Finset.sum_congr rfl
  intro c _
  exact quadraticPairLowCount_complement c w hw

theorem quadraticTotalWeightCount_64 : quadraticTotalWeightCount 64 = 1 := by
  have h := quadraticTotalWeightCount_complement 0 (by omega)
  simpa using h.symm.trans quadraticTotalWeightCount_0

theorem quadraticTotalWeightCount_48 : quadraticTotalWeightCount 48 = 2604 := by
  have h := quadraticTotalWeightCount_complement 16 (by omega)
  simpa using h.symm.trans quadraticTotalWeightCount_16

theorem quadraticTotalWeightCount_40 : quadraticTotalWeightCount 40 = 291648 := by
  have h := quadraticTotalWeightCount_complement 24 (by omega)
  simpa using h.symm.trans quadraticTotalWeightCount_24

theorem quadraticTotalWeightCount_36 : quadraticTotalWeightCount 36 = 888832 := by
  have h := quadraticTotalWeightCount_complement 28 (by omega)
  simpa using h.symm.trans quadraticTotalWeightCount_28

end Spin.UniformMixer
