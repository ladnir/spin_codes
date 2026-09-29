import SpinCodes.Structured.ConcreteOuterNativeLinearity

open Spin.Structured Spin.Structured.ConcreteOuter Spin.Structured.ConcreteNativeFamily
open scoped symmDiff

example {k : ℕ} (seed : Seed k) (x y : LocalMessage k) :
    encode seed (fun i j => xor (x i j) (y i j))=fun j => xor (encode seed x j) (encode seed y j) :=
  encode_xor seed x y
example (m : ℕ) (seed : NativeSeed m) (x y : Fin (Nsched m/2)→ZMod 2) :
    rows m seed (x+y)=fun i => rows m seed x i ∆ rows m seed y i := rows_add m seed x y
example (m : ℕ) (seed : NativeSeed m) (x y : Fin (Nsched m/2)→ZMod 2) :
    (nativeSetup m).outer seed (x+y)=(nativeSetup m).outer seed x+(nativeSetup m).outer seed y :=
  nativeSetup_outer_add m seed x y

#print axioms accF_xor
#print axioms golay_enc_xor
#print axioms encode_xor
#print axioms rowSupports_xor
#print axioms nativeRows_xor
#print axioms rows_add
#print axioms nativeSetup_outer_add
