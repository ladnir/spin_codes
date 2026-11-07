import SpinCodes.UniformMixer.Map22SubmaskParity

/-!
# Recover a 64-bit truth table from its Boolean ANF

This symbolic step turns the submask parity into the exact Möbius
inversion identity for any 64-bit mask. It connects the production
row-degree check to the bits of the actual expansion map.
-/

noncomputable section

namespace Spin.UniformMixer

def truthBit (mask : ℕ) (j : Fin 64) : ZMod 2 :=
  if mask.testBit j.val then 1 else 0

theorem anfCoefficient_as_sum (mask : ℕ) (m : Fin 64) :
    anfCoefficient mask m =
      ∑ k : Fin 64,
        if isSubmask k m then truthBit mask k else 0 := by
  unfold anfCoefficient
  apply Finset.sum_congr rfl
  intro k _
  by_cases hsub : isSubmask k m <;>
    by_cases hbit : mask.testBit k.val <;>
    simp [isSubmask, truthBit, hbit]

theorem anf_reconstruct (mask : ℕ) (j : Fin 64) :
    (∑ m : Fin 64,
      if isSubmask m j then anfCoefficient mask m else 0) =
    truthBit mask j := by
  calc
    _ = ∑ m : Fin 64, ∑ k : Fin 64,
        if isSubmask k m ∧ isSubmask m j then truthBit mask k else 0 := by
      apply Finset.sum_congr rfl
      intro m _
      by_cases hm : isSubmask m j
      · calc
          (if isSubmask m j then anfCoefficient mask m else 0) =
              anfCoefficient mask m := ite_eq_left hm
          _ = ∑ k : Fin 64,
              if isSubmask k m then truthBit mask k else 0 :=
                anfCoefficient_as_sum mask m
          _ = ∑ k : Fin 64,
              if isSubmask k m ∧ isSubmask m j then
                truthBit mask k else 0 := by
                  apply Finset.sum_congr rfl
                  intro k _
                  by_cases hk : isSubmask k m
                  · rw [ite_eq_left hk, ite_eq_left ⟨hk, hm⟩]
                  · rw [ite_eq_right hk,
                      ite_eq_right (fun h => hk h.1)]
      · calc
          (if isSubmask m j then anfCoefficient mask m else 0) = 0 :=
            ite_eq_right hm
          _ = ∑ k : Fin 64,
              if isSubmask k m ∧ isSubmask m j then
                truthBit mask k else 0 := by
                  symm
                  apply Finset.sum_eq_zero
                  intro k _
                  exact ite_eq_right (fun h => hm h.2)
    _ = ∑ k : Fin 64, ∑ m : Fin 64,
        if isSubmask k m ∧ isSubmask m j then truthBit mask k else 0 := by
      rw [Finset.sum_comm]
    _ = ∑ k : Fin 64,
        (∑ m : Fin 64,
          if isSubmask k m ∧ isSubmask m j then (1 : ZMod 2) else 0) *
          truthBit mask k := by
      apply Finset.sum_congr rfl
      intro k _
      rw [Finset.sum_mul]
      apply Finset.sum_congr rfl
      intro m _
      split_ifs <;> simp
    _ = ∑ k : Fin 64,
        (if k = j then 1 else 0) * truthBit mask k := by
      apply Finset.sum_congr rfl
      intro k _
      rw [six_bit_submask_parity j k]
    _ = truthBit mask j := by simp

end Spin.UniformMixer
