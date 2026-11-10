import SpinCodes.Turbo.SparseListDefs22

/-! Core-only definitions of the packed (Kronecker) polynomial-matrix checks;
see `SparseListPacked22` for the design and soundness. Kernel-check modules
import only this file (no Mathlib), which keeps their resident memory small. -/
namespace Spin.Turbo

/-! ## Core-only definitions -/

/-- Horner packing of entry `(i,j)` of a polynomial-matrix list, base `2^w`. -/
def pkEntryL22 (w i j : Nat) : List (List (List Int)) → Nat
  | [] => 0
  | M :: P => ((M.getD i []).getD j 0).toNat + (pkEntryL22 w i j P <<< w)

def pkMatL22 (w : Nat) (P : List (List (List Int))) : List (List Nat) :=
  (List.range 11).map fun i => (List.range 11).map fun j => pkEntryL22 w i j P

/-- `Σ_{k < n} A[i][k] · B[k][j]`. -/
def pkDotL22 (A B : List (List Nat)) (i j : Nat) : Nat → Nat
  | 0 => 0
  | k + 1 => pkDotL22 A B i j k + (A.getD i []).getD k 0 * (B.getD k []).getD j 0

/-- Every entry is in `[0, 2^v)`. -/
def pkBoundL22 (v : Nat) (P : List (List (List Int))) : Bool :=
  P.all fun M => M.all fun r => r.all fun x => decide (0 ≤ x) && decide (x.toNat < 2 ^ v)

/-- Walk the base-`2^w` digits of `N` against `s · u_e[i][j]`. -/
def pkDigitsLeL22 (w : Nat) (s : Int) (i j : Nat) : Nat → List (List (List Int)) → Bool
  | _, [] => true
  | N, M :: us =>
    decide (((N % 2 ^ w : Nat) : Int) ≤ s * (M.getD i []).getD j 0) && pkDigitsLeL22 w s i j (N >>> w) us

/-- Row `i` of the product check. The packed product is truncated to its
`D+1` low digits before the walk (kernel memory: the walk shifts a number of
`w (D+1)` bits instead of `w (2D+1)`). -/
def pkRowCheckL22 (D w : Nat) (s : Int) (A B : List (List Nat))
    (u : List (List (List Int))) (i : Nat) : Bool :=
  (List.range 11).all fun j =>
    pkDigitsLeL22 w s i j (pkDotL22 A B i j 11 % 2 ^ (w * (D + 1))) u

/-- Entry `(i,j)` of the row check (finer kernel declarations: only the 22
packed entries `A[i][·]`, `B[·][j]` are built). -/
def pkEntryCheckL22 (D w : Nat) (s : Int) (A B : List (List Nat))
    (u : List (List (List Int))) (i j : Nat) : Bool :=
  pkDigitsLeL22 w s i j (pkDotL22 A B i j 11 % 2 ^ (w * (D + 1))) u

def pkProdCheckL22 (D w : Nat) (s : Int) (A B : List (List Nat))
    (u : List (List (List Int))) : Bool :=
  (List.range 11).all fun i => pkRowCheckL22 D w s A B u i

/-- Shapes, lengths, entry bounds and the no-carry width condition. -/
def pkHdrCheckL22 (D v w : Nat) (a b u : List (List (List Int))) : Bool :=
  ((decide (a.length = D + 1) && decide (b.length = D + 1)) &&
    (decide (u.length = D + 1) && u.all sparseMatShapeL22)) &&
  ((pkBoundL22 v a && pkBoundL22 v b) && decide (11 * (2 * D + 1) * (2 ^ v * 2 ^ v) < 2 ^ w))

/-- Packed analog of `sparsePolyMulCheckDL22`: `u` has `D+1` shaped matrices and
`Σ_{r ≤ d} a_r · b_{d-r} ≤ s · u_d` entrywise for every `d ≤ D`. -/
def pkPolyMulCheckL22 (D v w : Nat) (s : Int) (a b u : List (List (List Int))) : Bool :=
  pkHdrCheckL22 D v w a b u && pkProdCheckL22 D w s (pkMatL22 w a) (pkMatL22 w b) u

def pkRegionalPowerCheckL22 (D v w : Nat) (s : Int) (base prev : List (List (List Int))) :
    List (Bool × List (List (List Int)) × List (List (List Int))) → Bool
  | [] => true
  | (bit, sq, next) :: tail =>
    (pkPolyMulCheckL22 D v w s prev prev sq &&
      (if bit then pkPolyMulCheckL22 D v w s sq base next else true)) &&
      pkRegionalPowerCheckL22 D v w s base (if bit then next else sq) tail

namespace SparseListPacked22

theorem pkRowCheck_of_entries (D w : Nat) (s : Int) (A B : List (List Nat))
    (u : List (List (List Int))) (i : Nat)
    (h : ∀ j, j < 11 → pkEntryCheckL22 D w s A B u i j = true) :
    pkRowCheckL22 D w s A B u i = true := by
  simp only [pkRowCheckL22, List.all_eq_true, List.mem_range]
  exact h


theorem pkPolyMulCheck_of_rows (D v w : Nat) (s : Int) (a b u : List (List (List Int)))
    (hh : pkHdrCheckL22 D v w a b u = true)
    (hr : ∀ i, i < 11 → pkRowCheckL22 D w s (pkMatL22 w a) (pkMatL22 w b) u i = true) :
    pkPolyMulCheckL22 D v w s a b u = true := by
  simp only [pkPolyMulCheckL22, hh, Bool.true_and, pkProdCheckL22, List.all_eq_true,
    List.mem_range]
  exact hr

end SparseListPacked22
end Spin.Turbo
