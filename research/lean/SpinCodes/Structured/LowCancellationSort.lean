import Mathlib.Data.List.Sort

/-! A structurally recursive sort for kernel reduction. Only preservation of
the multiset is used; sortedness is independently checked on each certificate. -/
namespace Spin.Structured.LowCancellation

def mergeNat (xs ys : List Nat) : List Nat :=
  List.rec (motive := fun _ => List Nat → List Nat) (fun ys => ys)
    (fun x xs mergeTail ys =>
      List.rec (motive := fun _ => List Nat) (x :: xs)
        (fun y ys merged => if x ≤ y then x :: mergeTail (y :: ys) else y :: merged) ys) xs ys

theorem mergeNat_eq (xs ys : List Nat) : mergeNat xs ys = List.merge xs ys (· ≤ ·) := by
  induction xs generalizing ys with
  | nil => simp [mergeNat]
  | cons x xs ih =>
    induction ys with
    | nil => simp [mergeNat]
    | cons y ys ihy =>
      simp only [mergeNat, List.merge]
      by_cases h : x ≤ y
      · simp only [h, ↓reduceIte, decide_true]
        exact congrArg (x :: ·) (ih (y :: ys))
      · simp only [h, ↓reduceIte, decide_false, Bool.false_eq_true]
        simpa only [mergeNat] using congrArg (y :: ·) ihy

def sortFuel : Nat → List Nat → List Nat
  | 0, xs => xs
  | n + 1, xs => match xs with
    | [] => []
    | [x] => [x]
    | x :: y :: zs =>
      let xs := x :: y :: zs
      let k := xs.length / 2
      mergeNat (sortFuel n (xs.take k)) (sortFuel n (xs.drop k))

theorem sortFuel_perm (n : Nat) (xs : List Nat) : (sortFuel n xs).Perm xs := by
  induction n generalizing xs with
  | zero => exact List.Perm.refl _
  | succ n ih =>
    rcases xs with _ | ⟨x, _ | ⟨y, zs⟩⟩
    · exact List.Perm.refl _
    · exact List.Perm.refl _
    · simp only [sortFuel, mergeNat_eq]
      exact (List.merge_perm_append _).trans
        ((List.Perm.append (ih _) (ih _)).trans (by rw [List.take_append_drop]))

def certSort (xs : List Nat) : List Nat := sortFuel 14 xs

theorem certSort_perm (xs : List Nat) : (certSort xs).Perm xs := sortFuel_perm 14 xs

end Spin.Structured.LowCancellation
