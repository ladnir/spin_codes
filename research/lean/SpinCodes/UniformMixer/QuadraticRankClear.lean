import SpinCodes.UniformMixer.QuadraticRankMatrixBridge

/-!
# Eliminated polar columns vanish

After a pivot is removed and XORed into the other rows, every remaining
row has a zero in the processed column. This supplies the triangular
invariant needed to count the solutions of the six polar equations.
-/

namespace Spin.UniformMixer

theorem eliminatePolarColumn_bit_false (rows : List Nat) (bit : Nat)
    (row : Nat) (hrow : row ∈ eliminatePolarColumn bit rows) :
    row.testBit bit = false := by
  unfold eliminatePolarColumn at hrow
  cases hfind : rows.find? (fun r => r.testBit bit) with
  | none =>
      simp only [hfind] at hrow
      have hnone := (List.find?_eq_none.mp hfind) row hrow
      cases hbit : row.testBit bit <;> simp_all
  | some pivot =>
      simp only [hfind] at hrow
      have hpbit : pivot.testBit bit = true :=
        List.find?_some (p := fun r : Nat => r.testBit bit) hfind
      obtain ⟨source, _, rfl⟩ := List.mem_map.mp hrow
      by_cases hs : source.testBit bit
      · simp [hs, Nat.testBit_xor, hpbit]
      · simp [hs]

end Spin.UniformMixer
