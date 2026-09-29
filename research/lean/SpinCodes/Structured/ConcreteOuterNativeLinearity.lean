import SpinCodes.Structured.ConcreteOuterGolayLinearity
import SpinCodes.Structured.ConcreteRoutedOutput
import SpinCodes.Structured.ConcreteNativeFamily

noncomputable section
namespace Spin.Structured.ConcreteOuter
open Finset
open scoped symmDiff

lemma support_xor {n : ℕ} (x y : Fin n→Bool) :
    support (fun i => xor (x i) (y i))=support x ∆ support y := by
  ext i
  simp only [support,mem_filter,mem_univ,true_and,mem_symmDiff]
  cases x i <;> cases y i <;> decide

theorem rowSupports_xor {L k : ℕ} (seed : Seed k) (x y : Fin L→LocalMessage k) :
    rowSupports seed (fun i j t => xor (x i j t) (y i j t))=
      fun i => rowSupports seed x i ∆ rowSupports seed y i := by
  funext i
  simp only [rowSupports,encode_xor,support_xor]

/-- A single shared outer seed yields a binary linear native row encoder. -/
theorem nativeRows_xor (m : ℕ) (seed : NativeSeed m) (x y : Fin (Nsched m/2)→Bool) :
    nativeRows m seed (fun i => xor (x i) (y i))=
      fun i => nativeRows m seed x i ∆ nativeRows m seed y i := by
  have he : nativeMessageEquiv m (fun i => xor (x i) (y i))=
      fun i j t => xor (nativeMessageEquiv m x i j t) (nativeMessageEquiv m y i j t) := rfl
  unfold nativeRows
  rw [he,rowSupports_xor]
  funext i
  ext j
  simp only [mem_map_equiv,mem_symmDiff]

end Spin.Structured.ConcreteOuter

namespace Spin.Structured.ConcreteBinaryEncoding

lemma bitEquiv_symm_add (x y : ZMod 2) :
    bitEquiv.symm (x+y)=xor (bitEquiv.symm x) (bitEquiv.symm y) := by
  apply bitEquiv.injective
  simp only [Equiv.apply_symm_apply,bitEquiv_xor]

lemma wordEquiv_symm_add (n : ℕ) (x y : Fin n→ZMod 2) :
    (wordEquiv n).symm (x+y)=fun i => xor ((wordEquiv n).symm x i) ((wordEquiv n).symm y i) := by
  funext i
  exact bitEquiv_symm_add _ _

end Spin.Structured.ConcreteBinaryEncoding

namespace Spin.Structured.ConcreteNativeFamily
open ConcreteOuter ConcreteBinaryEncoding
open scoped symmDiff

/-- Binary linearity, directly on the actual Family message carrier. -/
theorem rows_add (m : ℕ) (seed : NativeSeed m) (x y : Fin (Nsched m/2)→ZMod 2) :
    rows m seed (x+y)=fun i => rows m seed x i ∆ rows m seed y i := by
  unfold rows
  rw [wordEquiv_symm_add,nativeRows_xor]

/-- The concrete setup's outer map preserves addition without any selection assumption. -/
theorem nativeSetup_outer_add (m : ℕ) (seed : NativeSeed m) (x y : Fin (Nsched m/2)→ZMod 2) :
    (nativeSetup m).outer seed (x+y)=(nativeSetup m).outer seed x+(nativeSetup m).outer seed y := by
  change rowsToWord (rows m seed (x+y))=rowsToWord (rows m seed x)+rowsToWord (rows m seed y)
  rw [rows_add,rowsToWord_xor]

end Spin.Structured.ConcreteNativeFamily

