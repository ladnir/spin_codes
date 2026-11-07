import SpinCodes.UniformMixer.QuadraticRankIterate

/-!
# Triangular invariant for descending elimination

Processed columns remain zero, including when later pivots are XORed into
remaining rows. At column zero every residual row is therefore zero.
-/

namespace Spin.UniformMixer

def rowsHighClear (start : Nat) (rows : List Nat) : Prop :=
  ∀ row ∈ rows, ∀ bit, start ≤ bit → row.testBit bit = false

theorem eliminatePolarColumn_preserves_bit (rows : List Nat)
    (pivotBit bit : Nat)
    (hclear : ∀ row ∈ rows, row.testBit bit = false) :
    ∀ row ∈ eliminatePolarColumn pivotBit rows,
      row.testBit bit = false := by
  intro row hrow
  unfold eliminatePolarColumn at hrow
  cases hfind : rows.find? (fun r => r.testBit pivotBit) with
  | none =>
      simp only [hfind] at hrow
      exact hclear row hrow
  | some pivot =>
      have hp : pivot ∈ rows := List.mem_of_find?_eq_some hfind
      simp only [hfind] at hrow
      obtain ⟨source, hs, rfl⟩ := List.mem_map.mp hrow
      have hsource := hclear source (List.mem_of_mem_erase hs)
      have hpivot := hclear pivot hp
      by_cases hb : source.testBit pivotBit
      · simp [hb, Nat.testBit_xor, hsource, hpivot]
      · simpa [hb] using hsource

theorem eliminatePolarColumn_highClear (n : Nat) (rows : List Nat)
    (hclear : rowsHighClear (n + 1) rows) :
    rowsHighClear n (eliminatePolarColumn n rows) := by
  intro row hr bit hbit
  rcases Nat.eq_or_lt_of_le hbit with rfl | hgt
  · exact eliminatePolarColumn_bit_false rows n row hr
  · exact eliminatePolarColumn_preserves_bit rows n bit
      (fun source hs => hclear source hs bit (by omega)) row hr

theorem eliminateDown_highClear (n : Nat) (rows : List Nat)
    (hclear : rowsHighClear n rows) :
    rowsHighClear 0 (eliminateDown n rows) := by
  induction n generalizing rows with
  | zero => simpa [eliminateDown] using hclear
  | succ n ih =>
      exact ih (eliminatePolarColumn n rows)
        (eliminatePolarColumn_highClear n rows hclear)

theorem rowsHighClear_zero_row (rows : List Nat)
    (hclear : rowsHighClear 0 rows) (row : Nat) (hr : row ∈ rows) :
    row = 0 := by
  apply Nat.eq_of_testBit_eq
  intro bit
  simpa using hclear row hr bit (Nat.zero_le bit)

end Spin.UniformMixer
