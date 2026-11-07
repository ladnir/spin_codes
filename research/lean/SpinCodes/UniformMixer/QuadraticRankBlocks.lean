import SpinCodes.UniformMixer.QuadraticRankActual

/-!
# Bounded blocks for the full quadratic-rank histogram

The kernel can check small rank blocks, but a direct 1,024-form reduction
exceeds its memory limit. This range form is equivalent to the original
checker and splits exactly into adjoining blocks.
-/

namespace Spin.UniformMixer

def fastRadicalRangeCount (start length radicalSize : Nat) : Nat :=
  ((List.range' start length).filter (fun q =>
    fastRadicalSize q = radicalSize)).length

theorem fastRadicalRangeCount_eq_block (start length radicalSize : Nat) :
    fastRadicalRangeCount start length radicalSize =
      fastRadicalBlockCount start length radicalSize := by
  unfold fastRadicalRangeCount fastRadicalBlockCount
  rw [List.range'_eq_map_range, List.filter_map, List.length_map]
  rfl

theorem fastRadicalRangeCount_add (start m n radicalSize : Nat) :
    fastRadicalRangeCount start (m + n) radicalSize =
      fastRadicalRangeCount start m radicalSize +
        fastRadicalRangeCount (start + m) n radicalSize := by
  unfold fastRadicalRangeCount
  rw [← List.range'_append_1, List.filter_append, List.length_append]

set_option maxRecDepth 100000 in
set_option maxHeartbeats 0 in
theorem first256_fast_radical_4 :
    fastRadicalRangeCount 0 256 4 = 196 := by
  decide +kernel

end Spin.UniformMixer
