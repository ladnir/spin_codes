/-
The extended binary Golay code `[24,12,8]` and its weight enumerator.

    G(u) = 1 + 759 u^8 + 2576 u^12 + 759 u^16 + u^24

The paper takes this as known.  It has only `2^12 = 4096` codewords, so it can
be *computed* rather than assumed, the same way the sparse Collatz certificate
was: by the Lean kernel, with no appeal to compiled evaluation.

Two things make the computation feasible.  The count is done over the list
enumeration `allWords 12` rather than over `Fin 12 -> Bool`, because
`Fintype.piFinset` reduction is far slower in the kernel (a *trivial* predicate
over `Fin 12 -> Bool` already takes a minute).  And the whole weight
distribution is accumulated in a single pass, so the 4096 encodings are
performed once instead of once per weight class.  `card_filter_ofFn_eq_countP`
then carries the result back to the tuple side, where the enumerator lives.
-/
import SpinCodes.Structured.AccTuple
import SpinCodes.Structured.Enumerator

set_option linter.unusedSectionVars false
set_option maxRecDepth 1000000

namespace Spin.Structured.Golay

open Finset Polynomial Spin.Structured

/-! ## The code -/

/-- Pointwise XOR of two words. -/
def xorW (x y : List Bool) : List Bool := List.zipWith xor x y

@[simp] lemma length_xorW (x y : List Bool) :
    (xorW x y).length = min x.length y.length := by
  simp [xorW]

/-- The 12 generator rows, in systematic form `[I | B]`. -/
def rows : List (List Bool) :=
   [[true, false, false, false, false, false, false, false, false, false, false, false, false, true, true, true, true, true, true, true, true, true, true, true],
    [false, true, false, false, false, false, false, false, false, false, false, false, true, true, true, false, true, true, true, false, false, false, true, false],
    [false, false, true, false, false, false, false, false, false, false, false, false, true, true, false, true, true, true, false, false, false, true, false, true],
    [false, false, false, true, false, false, false, false, false, false, false, false, true, false, true, true, true, false, false, false, true, false, true, true],
    [false, false, false, false, true, false, false, false, false, false, false, false, true, true, true, true, false, false, false, true, false, true, true, false],
    [false, false, false, false, false, true, false, false, false, false, false, false, true, true, true, false, false, false, true, false, true, true, false, true],
    [false, false, false, false, false, false, true, false, false, false, false, false, true, true, false, false, false, true, false, true, true, false, true, true],
    [false, false, false, false, false, false, false, true, false, false, false, false, true, false, false, false, true, false, true, true, false, true, true, true],
    [false, false, false, false, false, false, false, false, true, false, false, false, true, false, false, true, false, true, true, false, true, true, true, false],
    [false, false, false, false, false, false, false, false, false, true, false, false, true, false, true, false, true, true, false, true, true, true, false, false],
    [false, false, false, false, false, false, false, false, false, false, true, false, true, true, false, true, true, false, true, true, true, false, false, false],
    [false, false, false, false, false, false, false, false, false, false, false, true, true, false, true, true, false, true, true, true, false, false, false, true]]

/-- Encode a 12-bit message: XOR the selected generator rows. -/
def encList (m : List Bool) : List Bool :=
  (m.zip rows).foldl (fun acc p => if p.1 then xorW acc p.2 else acc)
    (List.replicate 24 false)

lemma length_rows_mem : ∀ r ∈ rows, r.length = 24 := by decide

lemma length_encList (m : List Bool) : (encList m).length = 24 := by
  unfold encList
  have key : ∀ (l : List (Bool × List Bool)) (acc : List Bool),
      acc.length = 24 → (∀ p ∈ l, p.2.length = 24) →
      (l.foldl (fun acc p => if p.1 then xorW acc p.2 else acc) acc).length = 24 := by
    intro l
    induction l with
    | nil => intro acc h _; exact h
    | cons p t ih =>
        intro acc h hl
        refine ih _ ?_ (fun q hq => hl q (List.mem_cons_of_mem _ hq))
        by_cases hp : p.1
        · simp [hp, h, hl p (List.mem_cons_self ..)]
        · simpa [hp] using h
  refine key _ _ (by simp) ?_
  intro p hp
  exact length_rows_mem p.2 (List.of_mem_zip hp).2

/-! ## The weight histogram, computed in one pass -/

/-- Increment bucket `k` of a histogram. -/
def bumpAt : Nat -> List Nat -> List Nat
  | _, [] => []
  | 0, x :: t => (x + 1) :: t
  | (n + 1), x :: t => x :: bumpAt n t

@[simp] lemma length_bumpAt (k : Nat) (h : List Nat) : (bumpAt k h).length = h.length := by
  induction h generalizing k with
  | nil => cases k <;> rfl
  | cons x t ih => cases k <;> simp [bumpAt, ih]

lemma getD_bumpAt (k i : Nat) (h : List Nat) (hi : i < h.length) :
    (bumpAt k h).getD i 0 = h.getD i 0 + (if k = i then 1 else 0) := by
  induction h generalizing k i with
  | nil => simp at hi
  | cons x t ih =>
      cases k with
      | zero =>
          cases i with
          | zero => simp [bumpAt]
          | succ j => simp [bumpAt, List.getD_cons_succ]
      | succ n =>
          cases i with
          | zero => simp [bumpAt]
          | succ j =>
              simp only [bumpAt, List.getD_cons_succ]
              rw [ih n j (by simpa using Nat.lt_of_succ_lt_succ hi)]
              simp

/-- A one-pass histogram fold agrees with counting. -/
lemma getD_foldl_bumpAt {α : Type*} (f : α → Nat) :
    ∀ (l : List α) (h : List Nat) (i : Nat), i < h.length →
      (l.foldl (fun h x => bumpAt (f x) h) h).getD i 0
        = h.getD i 0 + l.countP (fun x => f x == i) := by
  intro l
  induction l with
  | nil => intro h i _; simp
  | cons x t ih =>
      intro h i hi
      simp only [List.foldl_cons, List.countP_cons]
      rw [ih _ i (by simpa using hi), getD_bumpAt _ _ _ hi]
      by_cases hx : f x = i <;> simp [hx] <;> omega

/-- The weight distribution of the code, accumulated in one pass. -/
def hist : List Nat :=
  (allWords 12).foldl (fun h m => bumpAt (wtL (encList m)) h) (List.replicate 25 0)

/-- **The extended Golay weight distribution, verified by the kernel.** -/
theorem hist_eq :
    hist = [1, 0, 0, 0, 0, 0, 0, 0, 759, 0, 0, 0, 2576,
            0, 0, 0, 759, 0, 0, 0, 0, 0, 0, 0, 1] := by
  decide

end Spin.Structured.Golay
