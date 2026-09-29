import SpinCodes.Structured.ConcreteOuterNative

namespace Spin.Structured.ConcreteBinaryEncoding
open Finset

/-- The binary support convention: false is zero and true is one. -/
def bitEquiv : Bool ≃ ZMod 2 where
  toFun b := if b then 1 else 0
  invFun x := decide (x ≠ 0)
  left_inv b := by cases b <;> decide
  right_inv := by decide +kernel

@[simp] theorem bitEquiv_false : bitEquiv false = 0 := rfl
@[simp] theorem bitEquiv_true : bitEquiv true = 1 := rfl
@[simp] theorem bitEquiv_symm_zero : bitEquiv.symm 0 = false := rfl

def wordEquiv (n : ℕ) : (Fin n → Bool) ≃ (Fin n → ZMod 2) :=
  Equiv.piCongrRight (fun _ => bitEquiv)

@[simp] theorem wordEquiv_zero (n : ℕ) : wordEquiv n 0 = 0 := rfl
@[simp] theorem wordEquiv_symm_zero (n : ℕ) : (wordEquiv n).symm 0 = 0 := rfl

/-- Flatten rows in row-major order into the outer-word type used by Family. -/
def rowsToWord {L b : ℕ} (rows : Fin L → Finset (Fin b)) : Fin (L*b) → ZMod 2 :=
  fun k => bitEquiv (decide ((finProdFinEquiv.symm k).2 ∈ rows (finProdFinEquiv.symm k).1))

def wordToRows {L b : ℕ} (word : Fin (L*b) → ZMod 2) : Fin L → Finset (Fin b) :=
  fun i => univ.filter (fun j => word (finProdFinEquiv (i,j)) ≠ 0)

@[simp] theorem wordToRows_rowsToWord {L b : ℕ} (rows : Fin L → Finset (Fin b)) :
    wordToRows (rowsToWord rows) = rows := by
  funext i
  ext j
  simp [wordToRows, rowsToWord, bitEquiv]

@[simp] theorem rowsToWord_wordToRows {L b : ℕ} (word : Fin (L*b) → ZMod 2) :
    rowsToWord (wordToRows word) = word := by
  funext k
  simp only [rowsToWord, wordToRows, mem_filter, mem_univ, true_and,
    Prod.mk.eta, Equiv.apply_symm_apply]
  exact bitEquiv.apply_symm_apply (word k)

def rowsEquiv (L b : ℕ) : (Fin L → Finset (Fin b)) ≃ (Fin (L*b) → ZMod 2) where
  toFun := rowsToWord
  invFun := wordToRows
  left_inv := wordToRows_rowsToWord
  right_inv := rowsToWord_wordToRows

end Spin.Structured.ConcreteBinaryEncoding
