import SpinCodes.Structured.PackedMap
import SpinCodes.Structured.ConcreteMapData

/-! The paper's selected maps on packed binary words and on finite supports.
The basis data and finite inverse checks come from `ConcreteMapData`; the
inverse checks apply to every input by the proved packed-map algebra. -/

namespace Spin.Structured.ConcreteMaps
open PackedMap

theorem a_bound (q : ℕ) : eval aRows q < 2 ^ 128 :=
  eval_lt aRows (by simpa only [List.all_eq_true, decide_eq_true_eq] using a_rows_bounded) q

theorem cTranspose_bound (q : ℕ) : eval cTransposeRows q < 2 ^ 128 :=
  eval_lt cTransposeRows
    (by simpa only [List.all_eq_true, decide_eq_true_eq] using cTranspose_rows_bounded) q

theorem c_bound (q : ℕ) : eval cRows q < 2 ^ 19 :=
  eval_lt cRows (by simpa only [List.all_eq_true, decide_eq_true_eq] using c_rows_bounded) q

def A (q : Fin (2 ^ 19)) : Fin (2 ^ 128) := ⟨eval aRows q, a_bound q⟩
def Ctranspose (q : Fin (2 ^ 19)) : Fin (2 ^ 128) :=
  ⟨eval cTransposeRows q, cTranspose_bound q⟩
def C (x : Fin (2 ^ 128)) : Fin (2 ^ 19) := ⟨eval cRows x, c_bound x⟩

theorem A_injective : Function.Injective A := by
  intro q p h
  exact eval_injective aRows aLeftInverse a_inverse_checked (congrArg Fin.val h)

theorem Ctranspose_injective : Function.Injective Ctranspose := by
  intro q p h
  exact eval_injective cTransposeRows cTransposeLeftInverse cTranspose_inverse_checked
    (congrArg Fin.val h)

def CrightInverse (q : Fin (2 ^ 19)) : Fin (2 ^ 128) :=
  ⟨eval cRightInverse q, eval_lt cRightInverse
    (by simpa only [List.all_eq_true, decide_eq_true_eq] using cRightInverse_bounded) q⟩

theorem C_rightInverse (q : Fin (2 ^ 19)) : C (CrightInverse q) = q := by
  apply Fin.ext
  exact eval_inverse cRightInverse cRows c_inverse_checked q.isLt

theorem C_surjective : Function.Surjective C := fun q => ⟨CrightInverse q, C_rightInverse q⟩

def inputBasis (i : Fin 128) : Fin (2 ^ 128) :=
  ⟨2 ^ i.val, Nat.pow_lt_pow_right (by decide) i.isLt⟩

theorem C_inputBasis (i : Fin 128) : (C (inputBasis i)).val = cRows.getD i 0 :=
  eval_onehot cRows i

theorem C_inputBasis_weight (i : Fin 128) : weight 19 (C (inputBasis i)) = 5 := by
  have hi : (i : ℕ) < cRows.length := by simpa only [c_length] using i.isLt
  have hm : cRows.getD i 0 ∈ cRows := by
    rw [List.getD_eq_getElem _ _ hi]
    exact List.getElem_mem hi
  have hh : ∀ r ∈ cRows, weight 19 r = 5 := by
    simpa only [List.all_eq_true, beq_iff_eq] using c_columns_weight_checked
  rw [C_inputBasis]
  exact hh _ hm

theorem C_inputBasis_nonzero (i : Fin 128) : (C (inputBasis i)).val ≠ 0 := by
  intro h
  have hh := C_inputBasis_weight i
  rw [h, weight_zero] at hh
  contradiction

theorem C_inputBasis_distinct : Function.Injective (fun i => C (inputBasis i)) := by
  intro i j h
  have hi : (i : ℕ) < cRows.length := by simpa only [c_length] using i.isLt
  have hj : (j : ℕ) < cRows.length := by simpa only [c_length] using j.isLt
  have hh := congrArg Fin.val h
  rw [C_inputBasis, C_inputBasis, List.getD_eq_getElem _ _ hi,
    List.getD_eq_getElem _ _ hj] at hh
  exact Fin.ext (c_columns_distinct_checked.getElem_inj.mp hh)

def xorWord {n : ℕ} (q p : Fin (2 ^ n)) : Fin (2 ^ n) :=
  ⟨q.val ^^^ p.val, Nat.xor_lt_two_pow q.isLt p.isLt⟩

theorem A_xor (q p : Fin (2 ^ 19)) : A (xorWord q p) = xorWord (A q) (A p) := by
  apply Fin.ext
  exact eval_xor aRows q p

theorem C_xor (q p : Fin (2 ^ 128)) : C (xorWord q p) = xorWord (C q) (C p) := by
  apply Fin.ext
  exact eval_xor cRows q p

theorem Ctranspose_xor (q p : Fin (2 ^ 19)) :
    Ctranspose (xorWord q p) = xorWord (Ctranspose q) (Ctranspose p) := by
  apply Fin.ext
  exact eval_xor cTransposeRows q p

open scoped symmDiff

theorem supportEquiv_xor {n : ℕ} (q p : Fin (2 ^ n)) :
    supportEquiv n (xorWord q p) = supportEquiv n q ∆ supportEquiv n p :=
  support_xor n q p

theorem supportEquiv_symm_xor {n : ℕ} (q p : Finset (Fin n)) :
    (supportEquiv n).symm (q ∆ p) =
      xorWord ((supportEquiv n).symm q) ((supportEquiv n).symm p) := by
  apply (supportEquiv n).injective
  simp only [Equiv.apply_symm_apply, supportEquiv_xor]

noncomputable def Aset (q : Finset (Fin 19)) : Finset (Fin 128) :=
  supportEquiv 128 (A ((supportEquiv 19).symm q))

noncomputable def Cset (x : Finset (Fin 128)) : Finset (Fin 19) :=
  supportEquiv 19 (C ((supportEquiv 128).symm x))

noncomputable def CtransposeSet (q : Finset (Fin 19)) : Finset (Fin 128) :=
  supportEquiv 128 (Ctranspose ((supportEquiv 19).symm q))

theorem Aset_xor (q p : Finset (Fin 19)) : Aset (q ∆ p) = Aset q ∆ Aset p := by
  simp only [Aset, supportEquiv_symm_xor, A_xor, supportEquiv_xor]

theorem Cset_xor (q p : Finset (Fin 128)) : Cset (q ∆ p) = Cset q ∆ Cset p := by
  simp only [Cset, supportEquiv_symm_xor, C_xor, supportEquiv_xor]

theorem CtransposeSet_xor (q p : Finset (Fin 19)) :
    CtransposeSet (q ∆ p) = CtransposeSet q ∆ CtransposeSet p := by
  simp only [CtransposeSet, supportEquiv_symm_xor, Ctranspose_xor, supportEquiv_xor]

theorem Aset_empty : Aset ∅ = ∅ := by simpa using Aset_xor ∅ ∅
theorem Cset_empty : Cset ∅ = ∅ := by simpa using Cset_xor ∅ ∅
theorem CtransposeSet_empty : CtransposeSet ∅ = ∅ := by simpa using CtransposeSet_xor ∅ ∅

theorem Aset_injective : Function.Injective Aset :=
  (supportEquiv 128).injective.comp (A_injective.comp (supportEquiv 19).symm.injective)

theorem Cset_surjective : Function.Surjective Cset :=
  (supportEquiv 19).surjective.comp (C_surjective.comp (supportEquiv 128).symm.surjective)

theorem CtransposeSet_injective : Function.Injective CtransposeSet :=
  (supportEquiv 128).injective.comp
    (Ctranspose_injective.comp (supportEquiv 19).symm.injective)

end Spin.Structured.ConcreteMaps
