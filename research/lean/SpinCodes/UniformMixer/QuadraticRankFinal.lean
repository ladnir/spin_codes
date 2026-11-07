import SpinCodes.UniformMixer.QuadraticRankTriangular

/-!
# Six-column rank and the polar solution count

The packed polar rows have no coordinates above the six state bits.
Descending elimination therefore leaves only zero rows. The number of
successful pivots determines the radical size exactly.
-/

namespace Spin.UniformMixer

theorem polarRowContribution_high_bit (q : Nat) (i : Fin 6)
    (p : Fin 15) (bit : Nat) (hbit : 6 ≤ bit) :
    (polarRowContribution q i p).testBit bit = false := by
  have hleft : (pairLeft p).val ≠ bit := by
    have := (pairLeft p).isLt
    omega
  have hright : (pairRight p).val ≠ bit := by
    have := (pairRight p).isLt
    omega
  unfold polarRowContribution
  split_ifs <;>
    simp [Nat.testBit_xor, Nat.testBit_two_pow_of_ne,
      hleft, hright]

theorem fold_xor_high_bit (bit : Nat) (xs : List Nat)
    (hxs : ∀ x ∈ xs, x.testBit bit = false) :
    (xs.foldl Nat.xor 0).testBit bit = false := by
  have aux (acc : Nat) (hacc : acc.testBit bit = false)
      (ys : List Nat) (hys : ∀ x ∈ ys, x.testBit bit = false) :
      (ys.foldl Nat.xor acc).testBit bit = false := by
    induction ys generalizing acc with
    | nil => simpa using hacc
    | cons head tail ih =>
        simp only [List.foldl_cons]
        apply ih
        · simp [Nat.testBit_xor, hacc, hys head (by simp)]
        · intro x hx
          exact hys x (by simp [hx])
  exact aux 0 (by simp) xs hxs

theorem polarRows_highClear (q : Nat) :
    rowsHighClear 6 (polarRows q) := by
  intro row hr bit hbit
  obtain ⟨i, _, rfl⟩ := List.mem_ofFn.mp hr
  unfold polarRowMask
  apply fold_xor_high_bit
  intro x hx
  obtain ⟨p, _, rfl⟩ := List.mem_ofFn.mp hx
  exact polarRowContribution_high_bit q i p bit hbit

theorem solutionCount_of_highClear_zero (rows : List Nat)
    (hclear : rowsHighClear 0 rows) : solutionCount rows = 64 := by
  classical
  have hall (x : State6) : rowsVanish rows x := by
    intro row hr
    rw [rowsHighClear_zero_row rows hclear row hr]
    exact rowEval_zero x
  unfold solutionCount
  have hfilter : Finset.univ.filter (rowsVanish rows) = Finset.univ := by
    ext x
    simp [hall x]
  rw [hfilter]
  decide

theorem solutionCount_rank (rows : List Nat)
    (hlen : rows.length = 6)
    (hclear : rowsHighClear 6 rows) :
    solutionCount rows =
      2 ^ (6 - (6 - (eliminateDown 6 rows).length)) := by
  let p := pivotCountDown 6 rows
  have hp := eliminateDown_length_add_pivots 6 rows
  have hp_le : p ≤ 6 := by omega
  have hp_eq : p = 6 - (eliminateDown 6 rows).length := by omega
  have hfinal := solutionCount_of_highClear_zero
    (eliminateDown 6 rows) (eliminateDown_highClear 6 rows hclear)
  have hcount := solutionCount_eliminateDown 6 (by omega) rows
  rw [hfinal] at hcount
  have hpow : 64 = 2 ^ p * 2 ^ (6 - p) := by
    calc
      64 = 2 ^ 6 := by decide
      _ = 2 ^ (p + (6 - p)) := by congr 1; omega
      _ = 2 ^ p * 2 ^ (6 - p) := pow_add 2 p (6 - p)
  have heq : 2 ^ p * solutionCount rows = 2 ^ p * 2 ^ (6 - p) := by
    exact hcount.symm.trans hpow
  have hres : solutionCount rows = 2 ^ (6 - p) :=
    mul_left_cancel₀ (by positivity : 2 ^ p ≠ 0) heq
  rw [hp_eq] at hres
  exact hres

theorem polar_solutionCount_rank (q : Nat) :
    solutionCount (polarRows q) = fastRadicalSize q := by
  have hlen : (polarRows q).length = 6 := by simp [polarRows]
  have h := solutionCount_rank (polarRows q) hlen (polarRows_highClear q)
  unfold fastRadicalSize polarRank
  rw [← eliminateDown_eq_fold]
  exact h

end Spin.UniformMixer
