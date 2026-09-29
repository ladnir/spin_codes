import SpinCodes.Structured.ConcreteNativeLinearCodeDistance

open Spin.Structured Spin.Structured.ConcreteOuter Spin.Structured.ConcreteNativeFamily

example (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    Module.finrank (ZMod 2) (realizedCode m out inner)=Nsched m/2 := realizedCode_finrank m out inner
example (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) :
    (Module.finrank (ZMod 2) (realizedCode m out inner):ℝ)/(Nsched m:ℝ)=1/2 := realizedCode_rate m out inner
example (m : ℕ) (out : NativeSeed m) (inner : InnerSeed m) (d : ℕ) :
    d < minimumDistance m out inner ↔
      ∀ x∈codewords m out inner,∀ y∈codewords m out inner,x≠y → d < hammingDist x y :=
  minimumDistance_gt_iff m out inner d

#print axioms codeword_add
#print axioms encoderLinearMap
#print axioms realizedCode_finrank
#print axioms realizedCode_rate
#print axioms low_weight_iff_close_pair
#print axioms nonzeroCodewords_nonempty
#print axioms minimumDistance_le_iff_pair
#print axioms minimumDistance_gt_iff
#print axioms concrete_probBad_minimumDistance
