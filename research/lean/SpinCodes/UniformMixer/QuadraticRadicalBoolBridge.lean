import SpinCodes.UniformMixer.QuadraticRadicalPacked
import SpinCodes.UniformMixer.QuadraticSpectrumBridge

/-!
# Boolean polar equations and the formal quadratic radical

This bridge expresses the six basis tests with Boolean XORs. It avoids
the expensive `ZMod 2` reductions in finite certificates while keeping
the checked predicate tied to the mathematical polar map.
-/

namespace Spin.UniformMixer

theorem bitScalar_and (a b : Bool) :
    bitScalar (a && b) = bitScalar a * bitScalar b := by
  cases a <;> cases b <;> decide

theorem bitScalar_decide_eq {α : Type*} [DecidableEq α] (a b : α) :
    bitScalar (decide (a = b)) = if a = b then 1 else 0 := by
  by_cases h : a = b <;> simp [bitScalar, h]

theorem bitScalar_foldl_xor (xs : List Bool) (initial : Bool) :
    bitScalar (xs.foldl Bool.xor initial) =
      bitScalar initial + (xs.map bitScalar).sum := by
  induction xs generalizing initial with
  | nil => simp
  | cons b bs ih =>
      rw [List.foldl_cons, ih, bitScalar_xor]
      simp [add_assoc]

def polarTermBool (q x : Nat) (i : Fin 6) (p : Fin 15) : Bool :=
  q.testBit p.val && Bool.xor
    (decide (pairLeft p = i) && x.testBit (pairRight p).val)
    (decide (pairRight p = i) && x.testBit (pairLeft p).val)

def polarBit (q x : Nat) (i : Fin 6) : Bool :=
  (List.ofFn (fun p : Fin 15 => polarTermBool q x i p)).foldl
    Bool.xor false

theorem polarTermBool_correct (q : Nat) (x : Fin 64)
    (i : Fin 6) (p : Fin 15) :
    bitScalar (polarTermBool q x.val i p) =
      bitScalar (q.testBit p.val) *
        (outputIndexState x (pairLeft p) * basis6 i (pairRight p) +
          basis6 i (pairLeft p) * outputIndexState x (pairRight p)) := by
  simp only [polarTermBool, bitScalar_and, bitScalar_xor,
    outputIndexState]
  by_cases hl : pairLeft p = i <;>
    by_cases hr : pairRight p = i <;>
      simp [basis6, bitScalar, hl, hr,
        mul_add, mul_comm]

theorem polarBit_correct (q : Nat) (x : Fin 64) (i : Fin 6) :
    bitScalar (polarBit q x.val i) =
      quadraticPolarLinear
        (fun p : Fin 15 => bitScalar (q.testBit p.val))
        (outputIndexState x) (basis6 i) := by
  rw [polarBit, bitScalar_foldl_xor]
  simp only [List.map_ofFn, List.sum_ofFn,
    quadraticPolarLinear, LinearMap.coe_mk, AddHom.coe_mk]
  rw [show bitScalar false = (0 : ZMod 2) by decide, zero_add]
  apply Finset.sum_congr rfl
  intro p _
  change bitScalar (polarTermBool q x.val i p) = _
  rw [polarTermBool_correct]
  ring

theorem bitScalar_eq_zero_iff (b : Bool) :
    bitScalar b = 0 ↔ b = false := by
  cases b <;> decide

theorem polarBit_radical_iff (q : Fin (2 ^ 15)) (x : Fin 64) :
    (∀ i : Fin 6, polarBit q.val x.val i = false) ↔
      outputIndexState x ∈ quadraticRadical (packedQuadratic q) := by
  rw [quadraticRadical_basis_mem]
  constructor
  · intro h i
    change quadraticPolarLinear
      (fun p : Fin 15 => bitScalar (q.val.testBit p.val))
      (outputIndexState x) (basis6 i) = 0
    rw [← polarBit_correct]
    exact (bitScalar_eq_zero_iff _).2 (h i)
  · intro h i
    have hi := h i
    change quadraticPolarLinear
      (fun p : Fin 15 => bitScalar (q.val.testBit p.val))
      (outputIndexState x) (basis6 i) = 0 at hi
    rw [← polarBit_correct, bitScalar_eq_zero_iff] at hi
    exact hi

def booleanRadicalSize (q : Fin (2 ^ 15)) : Nat :=
  (Finset.univ.filter (fun x : Fin 64 =>
    decide (∀ i : Fin 6, polarBit q.val x.val i = false) = true)).card

theorem booleanRadicalSize_eq (q : Fin (2 ^ 15)) :
    booleanRadicalSize q = packedRadicalSize q := by
  classical
  rw [packedRadicalSize, basisRadical_eq_quadraticRadical]
  unfold booleanRadicalSize
  rw [quadraticRadical]
  rw [Finset.card_filter, Finset.card_filter]
  apply Fintype.sum_equiv outputIndexEquiv
  intro x
  have h := polarBit_radical_iff q x
  simp only [quadraticRadical, Finset.mem_filter, Finset.mem_univ,
    true_and] at h
  simp [h, outputIndexEquiv]

end Spin.UniformMixer
