import SpinCodes.Structured.ConcreteNativeCodeword
import SpinCodes.Structured.ConcreteOuterNativeLinearity

noncomputable section
namespace Spin.Structured.ConcreteNativeFamily
open ConcreteOuter ConcreteRoute ConcreteEncoder ConcreteBinaryEncoding Finset

set_option backward.isDefEq.respectTransparency false in
/-- The complete emitted encoder is additive for every fixed realization of all seeds. -/
theorem codeword_add (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m)
    (x y : Fin (Nsched m/2)→ZMod 2) :
    codeword m out inner (x+y)=codeword m out inner x+codeword m out inner y := by
  have h := ConcreteRoutedEncoder.word_xor (R := rounds m)
    (streamWiring (native_round_divisibility m)) (rows m out x) (rows m out y) inner
  funext i
  simpa only [codeword,rows_add,Pi.add_apply] using congrFun h (finCongr (output_length m).symm i)

/-- The concrete encoder, as a linear map over the binary field. -/
def encoderLinearMap (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (Fin (Nsched m/2)→ZMod 2) →ₗ[ZMod 2] (Fin (Nsched m)→ZMod 2) where
  toFun := codeword m out inner
  map_add' := codeword_add m out inner
  map_smul' c x := by
    have hc : c=0 ∨ c=1 := (by decide : ∀ c:ZMod 2,c=0 ∨ c=1) c
    rcases hc with rfl|rfl <;> simp [codeword_zero]

/-- The realized native code is the range of its actual linear encoder. -/
def realizedCode (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    Submodule (ZMod 2) (Fin (Nsched m)→ZMod 2) := (encoderLinearMap m out inner).range

@[simp] theorem mem_realizedCode (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m)
    (c : Fin (Nsched m)→ZMod 2) :
    c∈realizedCode m out inner ↔ c∈codewords m out inner := by
  simp only [realizedCode,LinearMap.mem_range,codewords,mem_image,mem_univ,true_and]
  rfl

/-- The finite code has the exact intended binary dimension. -/
theorem realizedCode_finrank (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    Module.finrank (ZMod 2) (realizedCode m out inner)=Nsched m/2 := by
  rw [realizedCode,LinearMap.finrank_range_of_inj (codeword_injective m out inner)]
  simp

lemma native_dimension_twice (m : ℕ) : Nsched m/2*2=Nsched m := by
  apply Nat.div_mul_cancel
  rw [Nsched_eq]
  apply dvd_mul_of_dvd_left
  unfold Lsched
  omega

/-- Every realized code has exact rate one half, before taking a limit. -/
theorem realizedCode_rate (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (Module.finrank (ZMod 2) (realizedCode m out inner):ℝ)/(Nsched m:ℝ)=1/2 := by
  rw [realizedCode_finrank]
  have hN : (0:ℝ)<Nsched m := by exact_mod_cast Nsched_pos m
  have he : (Nsched m/2:ℕ)* (2:ℝ)=(Nsched m:ℝ) := by exact_mod_cast native_dimension_twice m
  apply (div_eq_iff hN.ne').mpr
  linarith

end Spin.Structured.ConcreteNativeFamily



