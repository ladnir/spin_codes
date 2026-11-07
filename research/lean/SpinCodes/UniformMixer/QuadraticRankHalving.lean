import SpinCodes.UniformMixer.QuadraticRankFlip

/-!
# One pivot halves the polar solution count

For residual rows with a cleared pivot column, flipping that coordinate
pairs solutions of the pivot equation with violations. The paired sets
have equal cardinality.
-/

namespace Spin.UniformMixer

def flipCoordinateEquiv (bit : Fin 6) : State6 ≃ State6 where
  toFun := flipCoordinate bit
  invFun := flipCoordinate bit
  left_inv := flipCoordinate_involutive bit
  right_inv := flipCoordinate_involutive bit

theorem pivot_flip_zero_iff_nonzero (pivot : Nat) (bit : Fin 6)
    (x : State6) (hbit : pivot.testBit bit.val = true) :
    rowEval pivot (flipCoordinate bit x) = 0 ↔
      rowEval pivot x ≠ 0 := by
  rw [rowEval_flip_of_bit_true pivot bit x hbit]
  rcases zmod2_zero_or_one (rowEval pivot x) with h0 | h1
  · simp [h0]
  · simp [h1, show (1 : ZMod 2) + 1 = 0 by decide]

noncomputable def solutionCount (rows : List Nat) : Nat := by
  classical
  exact (Finset.univ.filter (rowsVanish rows)).card

noncomputable def pivotZeroCount (rows : List Nat) (pivot : Nat) : Nat := by
  classical
  exact (Finset.univ.filter (fun x : State6 =>
    rowsVanish rows x ∧ rowEval pivot x = 0)).card

noncomputable def pivotNonzeroCount (rows : List Nat) (pivot : Nat) : Nat := by
  classical
  exact (Finset.univ.filter (fun x : State6 =>
    rowsVanish rows x ∧ rowEval pivot x ≠ 0)).card

theorem pivotZeroCount_eq_nonzero (rows : List Nat) (pivot : Nat)
    (bit : Fin 6)
    (hclear : ∀ row ∈ rows, row.testBit bit.val = false)
    (hbit : pivot.testBit bit.val = true) :
    pivotZeroCount rows pivot = pivotNonzeroCount rows pivot := by
  classical
  unfold pivotZeroCount pivotNonzeroCount
  rw [Finset.card_filter, Finset.card_filter]
  apply Fintype.sum_equiv (flipCoordinateEquiv bit)
  intro x
  change (if rowsVanish rows x ∧ rowEval pivot x = 0 then 1 else 0) =
    if rowsVanish rows (flipCoordinate bit x) ∧
      rowEval pivot (flipCoordinate bit x) ≠ 0 then 1 else 0
  rw [← (rowsVanish_flip_of_clear rows bit hclear x)]
  have hflip := pivot_flip_zero_iff_nonzero pivot bit x hbit
  have htwo := pivot_flip_zero_iff_nonzero pivot bit
    (flipCoordinate bit x) hbit
  rw [flipCoordinate_involutive] at htwo
  simp only [htwo]

theorem solutionCount_partition (rows : List Nat) (pivot : Nat) :
    solutionCount rows = pivotZeroCount rows pivot +
      pivotNonzeroCount rows pivot := by
  classical
  unfold solutionCount pivotZeroCount pivotNonzeroCount
  rw [Finset.card_filter, Finset.card_filter, Finset.card_filter,
    ← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro x _
  by_cases hr : rowsVanish rows x <;>
    by_cases hp : rowEval pivot x = 0 <;>
      simp [hr, hp]

theorem solutionCount_twice_pivotZero (rows : List Nat)
    (pivot : Nat) (bit : Fin 6)
    (hclear : ∀ row ∈ rows, row.testBit bit.val = false)
    (hbit : pivot.testBit bit.val = true) :
    solutionCount rows = 2 * pivotZeroCount rows pivot := by
  rw [solutionCount_partition,
    ← pivotZeroCount_eq_nonzero rows pivot bit hclear hbit]
  omega

theorem solutionCount_eliminate (rows : List Nat) (bit pivot : Nat)
    (hfind : rows.find? (fun row => row.testBit bit) = some pivot) :
    solutionCount rows =
      pivotZeroCount (eliminatePolarColumn bit rows) pivot := by
  classical
  unfold solutionCount pivotZeroCount
  congr 1
  ext x
  simp only [Finset.mem_filter, Finset.mem_univ, true_and]
  rw [rowsVanish_eliminate_iff rows bit pivot hfind x]
  exact and_comm

theorem solutionCount_elimination_doubles (rows : List Nat)
    (bit : Fin 6) (pivot : Nat)
    (hfind : rows.find? (fun row => row.testBit bit.val) = some pivot) :
    solutionCount (eliminatePolarColumn bit.val rows) =
      2 * solutionCount rows := by
  have hpbit : pivot.testBit bit.val = true :=
    List.find?_some (p := fun row : Nat => row.testBit bit.val) hfind
  have hclear : ∀ row ∈ eliminatePolarColumn bit.val rows,
      row.testBit bit.val = false :=
    fun row hr => eliminatePolarColumn_bit_false rows bit.val row hr
  calc
    solutionCount (eliminatePolarColumn bit.val rows) =
        2 * pivotZeroCount (eliminatePolarColumn bit.val rows) pivot :=
      solutionCount_twice_pivotZero _ pivot bit hclear hpbit
    _ = 2 * solutionCount rows := by
      rw [solutionCount_eliminate rows bit.val pivot hfind]

end Spin.UniformMixer
