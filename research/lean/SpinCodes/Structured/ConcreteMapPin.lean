import SpinCodes.Structured.ConcreteMaps

/-! Supplemental pins for the concrete maps; the original theorem pins
remain unchanged. -/

namespace Spin.Structured.ConcreteMapPin
open ConcreteMaps
open scoped symmDiff

example : aRows.length = 19 := a_length
example : cTransposeRows.length = 19 := cTranspose_length
example : cRows.length = 128 := c_length
example : aRows.getD 0 0 = 0x55555555555555555555555555555555 := by decide
example : aRows.getD 18 0 = 0x4e824e82e428e4284e824e82e428e428 := by decide
example : cTransposeRows.getD 0 0 = 0x800003d868100423210234004341a004 := by decide
example : cTransposeRows.getD 18 0 = 0x128a40312091018038456bcd40011874 := by decide
example : Function.Injective Aset := Aset_injective
example : Function.Surjective Cset := Cset_surjective
example : Function.Injective CtransposeSet := CtransposeSet_injective
example (q p : Finset (Fin 19)) : Aset (q ∆ p) = Aset q ∆ Aset p := Aset_xor q p
example (x y : Finset (Fin 128)) : Cset (x ∆ y) = Cset x ∆ Cset y := Cset_xor x y
example : cRows.Nodup := c_columns_distinct_checked
example : cRows.all (fun r => PackedMap.weight 19 r == 5) = true := c_columns_weight_checked
example (i : Fin 128) : PackedMap.weight 19 (C (inputBasis i)) = 5 := C_inputBasis_weight i
example (i : Fin 128) : (C (inputBasis i)).val ≠ 0 := C_inputBasis_nonzero i
example : Function.Injective (fun i => C (inputBasis i)) := C_inputBasis_distinct

end Spin.Structured.ConcreteMapPin
