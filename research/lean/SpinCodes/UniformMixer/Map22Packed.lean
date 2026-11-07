import SpinCodes.UniformMixer.Map22
import SpinCodes.Structured.MapSpectrum

/-!
# Packed certificate interface for the production 22-to-64 map

The existing Structured SPIN spectrum checker accepts a list of binary
row masks. This module pins that list to the production map before using
the checker on small blocks of input states.
-/

namespace Spin.UniformMixer

def productionPackedRows : List Nat := [
  0xffffffffffffffff, 0xaaaaaaaaaaaaaaaa, 0xcccccccccccccccc,
  0xf0f0f0f0f0f0f0f0, 0xff00ff00ff00ff00, 0xffff0000ffff0000,
  0xffffffff00000000, 0x3693af0a9cc6faa0, 0xd882d88272287228,
  0x666669960ff00000, 0xb874edded1e28448, 0xa5ffcc96cc965a00,
  0x8e2b4de88ed4b2e8, 0x2de14b78ddee4488, 0xcca5cc5a69009600,
  0x9556599a300cfcc0, 0x8888888888888888, 0xa0a0a0a0a0a0a0a0,
  0xaa00aa00aa00aa00, 0xaaaa0000aaaa0000, 0xaaaaaaaa00000000,
  0xcccc0000cccc0000]

theorem productionPackedRows_eq : productionPackedRows = List.ofFn rowMask := by
  decide

theorem productionPacked_basis (i : Fin 22) :
    Spin.Structured.PackedMap.eval productionPackedRows (2 ^ i.val) =
      rowMask i := by
  rw [Spin.Structured.PackedMap.eval_onehot, productionPackedRows_eq]
  rw [List.getD_eq_getElem _ _ (by simp)]
  change (List.ofFn rowMask)[i.val] = rowMask i
  rw [List.getElem_ofFn]

theorem productionPacked_basis_coordinate (i : Fin 22) (j : Fin 64) :
    (if (Spin.Structured.PackedMap.eval productionPackedRows
      (2 ^ i.val)).testBit j.val then (1 : ZMod 2) else 0) =
        coordMap j (basisVector i) := by
  rw [productionPacked_basis, coordMap_basis]
  rfl

def firstPackedBlockHistogram : List Nat :=
  (List.range 129).map fun w =>
    if w = 0 then 1 else if w = 24 then 16 else
    if w = 32 then 222 else if w = 40 then 16 else
    if w = 64 then 1 else 0

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem firstPackedBlock_checked :
    Spin.Structured.MapSpectrum.block productionPackedRows 0 256 =
      firstPackedBlockHistogram := by
  decide +kernel

def highPackedBlockHistogram : List Nat :=
  (List.range 129).map fun w =>
    if w = 16 then 4 else if w = 24 then 16 else
    if w = 32 then 216 else if w = 40 then 16 else
    if w = 48 then 4 else 0

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem highPackedBlock_checked :
    Spin.Structured.MapSpectrum.block productionPackedRows (2 ^ 21) 256 =
      highPackedBlockHistogram := by
  decide +kernel

end Spin.UniformMixer
