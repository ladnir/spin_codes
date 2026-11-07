import SpinCodes.UniformMixer.QuadraticRadicalRank
import SpinCodes.UniformMixer.QuadraticRadicalBoolBridge

/-!
# Binary row operations for the quadratic polar matrix

The compact rank evaluator uses XOR row elimination. These lemmas give
its semantic invariant: XORing a pivot into a row preserves that row's
zero equation whenever the pivot equation holds. The full elimination
rank theorem remains to be proved.
-/

namespace Spin.UniformMixer

def rowEval (row : Nat) (x : State6) : ZMod 2 :=
  ∑ j : Fin 6, bitScalar (row.testBit j.val) * x j

theorem rowEval_zero (x : State6) : rowEval 0 x = 0 := by
  simp [rowEval, bitScalar]

theorem rowEval_xor (a b : Nat) (x : State6) :
    rowEval (Nat.xor a b) x = rowEval a x + rowEval b x := by
  unfold rowEval
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro j _
  change bitScalar ((a ^^^ b).testBit j.val) * x j =
    bitScalar (a.testBit j.val) * x j +
      bitScalar (b.testBit j.val) * x j
  rw [Nat.testBit_xor, bitScalar_xor, add_mul]

theorem rowEval_xor_iff (a pivot : Nat) (x : State6)
    (hp : rowEval pivot x = 0) :
    rowEval (Nat.xor a pivot) x = 0 ↔ rowEval a x = 0 := by
  rw [rowEval_xor, hp]
  simp

def rowsVanish (rows : List Nat) (x : State6) : Prop :=
  ∀ row ∈ rows, rowEval row x = 0

theorem rowsVanish_map_pivot (rows : List Nat) (pivot bit : Nat)
    (x : State6) (hp : rowEval pivot x = 0) :
    rowsVanish
      (rows.map (fun row => if row.testBit bit then Nat.xor row pivot else row)) x
      ↔ rowsVanish rows x := by
  simp only [rowsVanish, List.forall_mem_map]
  constructor
  · intro h row hr
    have hrow := h row hr
    by_cases hbit : row.testBit bit
    · exact (rowEval_xor_iff row pivot x hp).mp (by simpa [hbit] using hrow)
    · simpa [hbit] using hrow
  · intro h row hr
    have hrow := h row hr
    by_cases hbit : row.testBit bit
    · exact (by simpa [hbit] using
        (rowEval_xor_iff row pivot x hp).mpr hrow)
    · simpa [hbit] using hrow

theorem mem_erase_of_ne (rows : List Nat) (row pivot : Nat)
    (hne : row ≠ pivot) (hmem : row ∈ rows) :
    row ∈ rows.erase pivot := by
  induction rows with
  | nil => simp at hmem
  | cons head tail ih =>
      by_cases hhead : head = pivot
      · subst head
        rw [List.erase_cons_head]
        rcases List.mem_cons.mp hmem with heq | htail
        · exact False.elim (hne heq)
        · exact htail
      · rw [List.erase_cons_tail (by simpa using hhead)]
        rcases List.mem_cons.mp hmem with heq | htail
        · exact List.mem_cons.mpr (Or.inl heq)
        · exact List.mem_cons.mpr (Or.inr (ih htail))

theorem rowsVanish_erase_iff (rows : List Nat) (pivot : Nat)
    (hp : pivot ∈ rows) (x : State6) :
    rowsVanish rows x ↔
      rowEval pivot x = 0 ∧ rowsVanish (rows.erase pivot) x := by
  constructor
  · intro h
    refine ⟨h pivot hp, ?_⟩
    intro row hrow
    exact h row (List.mem_of_mem_erase hrow)
  · rintro ⟨hp0, hrest⟩ row hrow
    by_cases heq : row = pivot
    · simpa [heq] using hp0
    · exact hrest row (mem_erase_of_ne rows row pivot heq hrow)

theorem rowsVanish_eliminate_iff (rows : List Nat) (bit pivot : Nat)
    (hfind : rows.find? (fun row => row.testBit bit) = some pivot)
    (x : State6) :
    rowsVanish rows x ↔
      rowEval pivot x = 0 ∧
        rowsVanish (eliminatePolarColumn bit rows) x := by
  have hp : pivot ∈ rows := List.mem_of_find?_eq_some hfind
  have herase := rowsVanish_erase_iff rows pivot hp x
  have hmap (hp0 : rowEval pivot x = 0) :
      rowsVanish (eliminatePolarColumn bit rows) x ↔
        rowsVanish (rows.erase pivot) x := by
    simp only [eliminatePolarColumn, hfind]
    exact rowsVanish_map_pivot (rows.erase pivot) pivot bit x hp0
  constructor
  · intro h
    obtain ⟨hp0, hrest⟩ := herase.mp h
    exact ⟨hp0, (hmap hp0).mpr hrest⟩
  · rintro ⟨hp0, hrest⟩
    exact herase.mpr ⟨hp0, (hmap hp0).mp hrest⟩

theorem eliminatePolarColumn_length (rows : List Nat) (bit pivot : Nat)
    (hfind : rows.find? (fun row => row.testBit bit) = some pivot) :
    (eliminatePolarColumn bit rows).length = rows.length - 1 := by
  have hp : pivot ∈ rows := List.mem_of_find?_eq_some hfind
  simp [eliminatePolarColumn, hfind, List.length_erase_of_mem hp]

end Spin.UniformMixer
