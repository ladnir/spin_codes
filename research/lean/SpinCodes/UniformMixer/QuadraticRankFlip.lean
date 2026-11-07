import SpinCodes.UniformMixer.QuadraticRankClear

/-!
# Coordinate flips for polar solution counting

Flipping one input coordinate toggles precisely the row equations whose
mask has that bit set. After elimination the remaining rows have a zero
in the pivot column, so this flip pairs their solutions according to
the pivot equation.
-/

namespace Spin.UniformMixer

def flipCoordinate (bit : Fin 6) (x : State6) : State6 :=
  fun j => if j = bit then x j + 1 else x j

theorem flipCoordinate_involutive (bit : Fin 6) (x : State6) :
    flipCoordinate bit (flipCoordinate bit x) = x := by
  funext j
  by_cases hj : j = bit
  · simp [flipCoordinate, hj, add_assoc,
      show (1 : ZMod 2) + 1 = 0 by decide]
  · simp [flipCoordinate, hj]

theorem rowEval_flipCoordinate (row : Nat) (bit : Fin 6)
    (x : State6) :
    rowEval row (flipCoordinate bit x) =
      rowEval row x + bitScalar (row.testBit bit.val) := by
  unfold rowEval
  calc
    (∑ j : Fin 6, bitScalar (row.testBit j.val) *
        flipCoordinate bit x j) =
      ∑ j : Fin 6,
        (bitScalar (row.testBit j.val) * x j +
          if j = bit then bitScalar (row.testBit j.val) else 0) := by
            apply Finset.sum_congr rfl
            intro j _
            by_cases hj : j = bit <;>
              simp [flipCoordinate, hj, mul_add]
    _ = (∑ j : Fin 6, bitScalar (row.testBit j.val) * x j) +
        bitScalar (row.testBit bit.val) := by
          rw [Finset.sum_add_distrib]
          simp

theorem rowEval_flip_of_bit_false (row : Nat) (bit : Fin 6)
    (x : State6) (hbit : row.testBit bit.val = false) :
    rowEval row (flipCoordinate bit x) = rowEval row x := by
  rw [rowEval_flipCoordinate, hbit]
  simp [bitScalar]

theorem rowEval_flip_of_bit_true (row : Nat) (bit : Fin 6)
    (x : State6) (hbit : row.testBit bit.val = true) :
    rowEval row (flipCoordinate bit x) = rowEval row x + 1 := by
  rw [rowEval_flipCoordinate, hbit]
  simp [bitScalar]

theorem rowsVanish_flip_of_clear (rows : List Nat) (bit : Fin 6)
    (hclear : ∀ row ∈ rows, row.testBit bit.val = false)
    (x : State6) :
    rowsVanish rows (flipCoordinate bit x) ↔ rowsVanish rows x := by
  constructor <;> intro h row hr
  · have hrow := h row hr
    rwa [rowEval_flip_of_bit_false row bit x (hclear row hr)] at hrow
  · rw [rowEval_flip_of_bit_false row bit x (hclear row hr)]
    exact h row hr

end Spin.UniformMixer
