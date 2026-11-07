import SpinCodes.UniformMixer.QuadraticRankHalving

/-!
# Iterating polar elimination

The pivot count records successful eliminations in descending column order.
Each success doubles the solution count of the residual equations and removes
one row.
-/

namespace Spin.UniformMixer

def eliminateDown : Nat → List Nat → List Nat
  | 0, rows => rows
  | n + 1, rows => eliminateDown n (eliminatePolarColumn n rows)

def pivotCountDown : Nat → List Nat → Nat
  | 0, _ => 0
  | n + 1, rows =>
      (if (rows.find? (fun row => row.testBit n)).isSome then 1 else 0) +
        pivotCountDown n (eliminatePolarColumn n rows)

theorem eliminateDown_eq_fold (n : Nat) (rows : List Nat) :
    eliminateDown n rows =
      (List.range n).reverse.foldl
        (fun rs bit => eliminatePolarColumn bit rs) rows := by
  induction n generalizing rows with
  | zero => simp [eliminateDown]
  | succ n ih =>
      simp only [eliminateDown, List.range_succ, List.reverse_append,
        List.reverse_cons, List.reverse_nil, List.nil_append,
        ]
      exact ih _

theorem eliminateDown_length_add_pivots (n : Nat) (rows : List Nat) :
    (eliminateDown n rows).length + pivotCountDown n rows = rows.length := by
  induction n generalizing rows with
  | zero => simp [eliminateDown, pivotCountDown]
  | succ n ih =>
      cases hfind : rows.find? (fun row => row.testBit n) with
      | none => simpa [eliminateDown, pivotCountDown, hfind,
          eliminatePolarColumn] using ih rows
      | some pivot =>
          have hp : pivot ∈ rows := List.mem_of_find?_eq_some hfind
          have hlen := eliminatePolarColumn_length rows n pivot hfind
          have hrec := ih (eliminatePolarColumn n rows)
          simp only [eliminateDown, pivotCountDown, hfind, Option.isSome_some,
            ite_true] at *
          have hpos : 0 < rows.length := List.length_pos_of_mem hp
          omega

theorem solutionCount_eliminateDown (n : Nat) (hn : n ≤ 6)
    (rows : List Nat) :
    solutionCount (eliminateDown n rows) =
      2 ^ (pivotCountDown n rows) * solutionCount rows := by
  induction n generalizing rows with
  | zero => simp [eliminateDown, pivotCountDown]
  | succ n ih =>
      have hbit : n < 6 := by omega
      have hrec := ih (by omega) (eliminatePolarColumn n rows)
      cases hfind : rows.find? (fun row => row.testBit n) with
      | none =>
          simpa [eliminateDown, pivotCountDown, hfind,
            eliminatePolarColumn] using hrec
      | some pivot =>
          have hdouble := solutionCount_elimination_doubles rows
            ⟨n, hbit⟩ pivot hfind
          simp only [eliminateDown, pivotCountDown, hfind,
            Option.isSome_some, ite_true]
          rw [hrec, hdouble]
          simp only [pow_add, pow_one]
          ring

end Spin.UniformMixer
