import SpinCodes.Turbo.PhysicalBirth22
import SpinCodes.UniformMixer.Map22PacketProfileSmallDefs

/-!
# Packed actual physical inputs for the small birth census

Historical audit: Map22PackedBridge proves coordinatewise agreement with the
production expansion; the spectrum block checker separates packed evaluation
from fast weight checking. The same bridges here cover A^T x and the zero-state
output weight. The exact occupancy law remains PhysicalStep22's uniform layer.
-/
namespace Spin.Turbo
open Spin.UniformMixer
open Spin.Structured

def smallBirthFeedbackRows22 : List Nat :=
  PackedMap.transpose 64 productionPackedRows

set_option maxRecDepth 100000 in
set_option maxHeartbeats 1000000 in
theorem smallBirthFeedbackRows22_coordinate :
    ∀ (j : Fin 64) (i : Fin 22),
      bitScalar ((smallBirthFeedbackRows22.getD j.val 0).testBit i.val) = coordCoeff j i := by
  decide +kernel

def smallBirthPackedInput22 (x : Nat) : Input64 :=
  fun j => bitScalar (x.testBit j.val)

theorem smallBirthPackedFeedback22 (x : Nat) :
    packedState (PackedMap.eval smallBirthFeedbackRows22 x) =
      stepFeedback (smallBirthPackedInput22 x) := by
  funext i
  rw [packedState, packed_eval_coordinate]
  have hlen : smallBirthFeedbackRows22.length = 64 := by decide
  change (∑ j : Fin smallBirthFeedbackRows22.length,
    bitScalar ((smallBirthFeedbackRows22.getD j.val 0).testBit i.val) *
      bitScalar (x.testBit j.val)) = _
  rw [show smallBirthFeedbackRows22.length = 64 from hlen]
  unfold stepFeedback smallBirthPackedInput22
  apply Finset.sum_congr rfl
  intro j _
  rw [smallBirthFeedbackRows22_coordinate]
  exact mul_comm _ _

theorem smallBirthPackedWeight22 (x : Nat) :
    PackedMap.weight 64 x = stepWeight (smallBirthPackedInput22 x) 0 := by
  rw [PackedMap.weight_eq_sum_bits]
  unfold stepWeight smallBirthPackedInput22
  apply Finset.sum_congr rfl
  intro j _
  rw [packedBit_weight]
  simp only [map_zero, zero_add]

theorem smallBirthPackedExpansionWeight22 (x : Nat) :
    PackedMap.weight 64 (PackedMap.eval productionPackedRows
      (PackedMap.eval smallBirthFeedbackRows22 x)) =
      coordinateWeight coordMap (stepFeedback (smallBirthPackedInput22 x)) := by
  rw [productionPacked_weight, smallBirthPackedFeedback22]

theorem smallBirthPackedInput22_xor (x y : Nat) :
    smallBirthPackedInput22 (x ^^^ y) = smallBirthPackedInput22 x + smallBirthPackedInput22 y := by
  funext j
  change bitScalar ((x ^^^ y).testBit j.val) = _
  rw [Nat.testBit_xor, bitScalar_xor]
  rfl
end Spin.Turbo



